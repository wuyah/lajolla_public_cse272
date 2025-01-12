#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    
    // Homework 1: implement this!
    // check incoming light's position
    Real eta = dot( vertex.geometric_normal, dir_in ) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Spectrum baseColor = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * eta);
    }
    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    Real n_dot_in = dot( dir_in, frame.n );
    Real h_dot_out = dot( dir_out, half_vector );
    Real h_dot_in = dot( dir_in, half_vector );

    Real F_g = fresnel_dielectric( fabs(h_dot_in), fabs(h_dot_out),eta );

    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real r2 = roughness * roughness;
    Real anisotropic = eval( bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool );
    Real aspect = sqrt( 1 - 0.9 * anisotropic );
    Real alpha_x = fmax( 0.0001, r2 / aspect );
    Real alpha_y = fmax( 0.0001, r2 * aspect );

    Vector3 h_l = to_local(frame, half_vector);
    Real D_g = 1 / (c_PI  * alpha_x * alpha_y *
                pow( (pow(h_l.x / alpha_x, 2) + pow(h_l.y / alpha_y, 2) + pow(h_l.z, 2)), 2) );

    Real G_g = smith_masking_disney(to_local(frame, dir_in), Vector2(alpha_x, alpha_y)) *
                smith_masking_disney(to_local(frame, dir_out), Vector2(alpha_x, alpha_y));

    Spectrum f_g = make_const_spectrum(Real(0));
    if(reflect) {
        f_g = baseColor * F_g * D_g * G_g / (4 * fabs(n_dot_in));
    } else {
        if (h_dot_in < 0) {
            h_dot_in = -h_dot_in;
            h_dot_out = -h_dot_out;
        }
        // exit(-1);
        f_g =  sqrt(baseColor) * (1-F_g) * D_g * G_g * fabs(h_dot_in * h_dot_out) /
            (fabs(n_dot_in) * pow(h_dot_in + eta * h_dot_out, Real(2)) );
    }

    return f_g;
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real r2 = roughness * roughness;

    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool );

    Real aspect = sqrt( 1 - 0.9 * anisotropic );
    Real alpha_x = fmax( 0.0001, r2 / aspect );
    Real alpha_y = fmax( 0.0001, r2 * aspect );

    // We sample the visible normals, also we use F to determine
    // whether to sample reflection or refraction
    // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
    Real h_dot_in = dot(half_vector, dir_in);
    Real n_dot_in = dot( dir_in, frame.n );
    Real h_dot_out = dot( dir_out, half_vector );

    Real F_g = fresnel_dielectric(fabs(h_dot_in), fabs(h_dot_out), eta);
    Vector3 h_l = to_local(frame, half_vector);
    Real D_g = 1 / (c_PI  * alpha_x * alpha_y *
            pow( (pow(h_l.x / alpha_x, 2) + pow(h_l.y / alpha_y, 2) + pow(h_l.z, 2)), 2)
            );

    Real G_in = smith_masking_disney(to_local(frame, dir_in), Vector2(alpha_x, alpha_y));
    if (reflect) {
        return (F_g * D_g * G_in) / (4 * fabs(n_dot_in));
    } else {
        Real h_dot_out = dot(half_vector, dir_out);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
        return (1 - F_g) * D_g * G_in * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    // Sample a micro normal and transform it to world space -- this is our half-vector.
    Real r2 = roughness * roughness;

    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool );

    Real aspect = sqrt( 1 - 0.9 * anisotropic );
    Real alpha_x = fmax( 0.0001, r2 / aspect );
    Real alpha_y = fmax( 0.0001, r2 * aspect );

    Vector3 local_dir_in = to_local(frame, dir_in);
    Vector3 local_micro_normal =
        sample_visible_normals_ani(local_dir_in, Vector2(alpha_x, alpha_y) , rnd_param_uv);

    Vector3 half_vector = normalize( to_world(frame, local_micro_normal) );
    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    // Now we need to decide whether to reflect or refract.
    // We do this using the Fresnel term.
    Real h_dot_in = dot(half_vector, dir_in);
    // Real h_dot_out = dot( dir_out, half_vector );
    // calculate h_dot_out by using Snell's law:
    if (h_dot_in < 0) {
        half_vector = -half_vector;
    }

    Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
    if (h_dot_out_sq <= 0) {
        // Handle total internal reflection
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        // set eta to 0 since we are not transmitting
        return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
    }
    Real h_dot_out = sqrt(h_dot_out_sq);

    Real F_g = fresnel_dielectric(fabs(h_dot_in), fabs(h_dot_out),eta);
    if (rnd_param_w <= F_g) {
        // Reflection
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        // set eta to 0 since we are not transmitting
        return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
    } else {
        // Refraction
        // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
        // (note that our eta is eta2 / eta1, and l = -dir_in)

        // flip half_vector if needed
        Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
        assert( dot(refracted, half_vector) * dot(dir_in, half_vector) < 0  );
        return BSDFSampleRecord{refracted, eta, roughness};
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
