#include "../microfacet.h"

// The only different to roughplastic is introducing the anistropic to extend the normal to different spaces.

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Vector3 half_vector = normalize(dir_in + dir_out);
    Vector3 h_l = to_local(frame, half_vector);
    
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real r2 = pow(roughness, 2);

    Real n_dot_out = fmax((dot(dir_out, frame.n)) ,Real(0));
    Real n_dot_in = fmax(dot(dir_in, frame.n), Real(0));
    Real out_dot_h = fmax(dot(dir_out, half_vector), Real(0));

    Spectrum base_color = eval(
        bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool );
    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool );
    anisotropic = std::clamp(roughness, Real(0.01), Real(1));

    Vector3 F_m = base_color + (1-base_color) * pow(1 - out_dot_h, Real(5));


    Real aspect = sqrt( 1 - 0.9 * anisotropic );
    Real alpha_x = fmax( 0.0001, r2 / aspect );
    Real alpha_y = fmax( 0.0001, r2 * aspect );

    Real D_m = 1 / (c_PI  * alpha_x * alpha_y *
                pow( (pow(h_l.x / alpha_x, 2) + pow(h_l.y / alpha_y, 2) + pow(h_l.z, 2)), 2)
                );
    Real G = smith_masking_disney(to_local(frame, dir_in), Vector2(alpha_x, alpha_y)) *
             smith_masking_disney(to_local(frame, dir_out), Vector2(alpha_x, alpha_y));

    Vector3 f_metal = F_m * D_m * G / (4 * n_dot_in);
    return f_metal;
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return 0;
    }

    // Spectrum S = eval(
    //     bsdf., vertex.uv, vertex.uv_screen_size, texture_pool);
    // Spectrum R = eval(
    //     bsdf.diffuse_reflectance, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Real lS = luminance(S), lR = luminance(R);
    // if (lS + lR <= 0) {
    //     return 0;
    // }
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval( bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool );
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.001), Real(1));
    anisotropic = std::clamp(anisotropic, Real(0.001), Real(1));


    Real aspect = sqrt( 1 - 0.9 * anisotropic );
    Real alpha_x = fmax( 0.0001, pow(roughness, 2) / aspect );
    Real alpha_y = fmax( 0.0001, pow(roughness, 2) * aspect );


    Real G_in = smith_masking_disney(to_local(frame, dir_in), Vector2(alpha_x, alpha_y));

    Vector3 h_l = to_local(frame, half_vector);
    Real D_m = 1 / (c_PI  * alpha_x * alpha_y *
                pow( (pow(h_l.x / alpha_x, 2) + pow(h_l.y / alpha_y, 2) + pow(h_l.z, 2)), 2 )
                );
    // (4 * cos_theta_v) is the Jacobian of the reflectiokn
    Real spec_prob = (G_in * D_m) / (4 * n_dot_in);
    // For the diffuse lobe, we importance sample cos_theta_out
    return spec_prob;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval( bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool );
    Real aspect = sqrt( 1 - 0.9 * anisotropic );
    Real alpha_x = fmax( 0.0001, pow(roughness, 2) / aspect );
    Real alpha_y = fmax( 0.0001, pow(roughness, 2) * aspect );
        
    // Transform the micro normal to world space
    Vector3 half_vector = to_world(frame, 
            sample_visible_normals_metal( to_local(frame, dir_in), Vector2(alpha_x, alpha_y), rnd_param_uv));

    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);

    return BSDFSampleRecord{ reflected,
        Real(0) /* eta */, roughness /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
