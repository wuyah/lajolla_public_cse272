#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneySheen &bsdf) const {
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
    Vector3 half_vector = normalize( dir_in + dir_out );

    Real h_dot_out = dot( half_vector, dir_out) ;
    Real n_dot_out = dot( frame.n, dir_out );

    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real L = luminance(base_color);
    Spectrum c_tint = L > 0 ? base_color / L : make_const_spectrum(1);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum c_sheen = (1 - sheen_tint) + sheen_tint * c_tint;
    Spectrum f_tint = c_sheen * pow( 1 - fabs(h_dot_out), Real(5) ) * fabs(n_dot_out);

    // Homework 1: implement this!
    return f_tint;
}

Real pdf_sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
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
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
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
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(1)};  /* roughness */
}

TextureSpectrum get_texture_op::operator()(const DisneySheen &bsdf) const {
    return bsdf.base_color;
}
