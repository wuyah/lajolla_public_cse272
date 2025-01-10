#include "../microfacet.h"

#define DISNEY_CC_ROUGHNESS Real(0.25)

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real out_dot_h = fmax(dot(dir_out, half_vector), Real(0));
    Real n_dot_in = fmax(dot( dir_in, vertex.geometric_normal), Real(0));
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool );
    clearcoat_gloss = std::clamp(clearcoat_gloss, Real(0.001), Real(1));

    Real alpha_g = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
    Real a2 = alpha_g * alpha_g;

    // hard code number for hardcode \eta = 1.5 and R_0 = (\eta - 1)^2 / (\eta + 1)^2
    Real R0_eta = Real(0.04);

    //  R0_eta + (1-R0_eta) * pow(1 - out_dot_h, Real(5));
    Real F_c = schlick_fresnel(R0_eta, out_dot_h);

    Real D_c = (a2 - 1) / (c_PI  * log(a2) * (1 + (a2 - 1) * h_l.z * h_l.z));
    Real G_c = smith_masking_gtr2(to_local(frame, dir_in), DISNEY_CC_ROUGHNESS) *
             smith_masking_gtr2(to_local(frame, dir_out), DISNEY_CC_ROUGHNESS);

    Real f = F_c * D_c * G_c / (4 * n_dot_in);
    
    return make_const_spectrum(fmax(f, 0));
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Vector3 h_l = to_local(frame, half_vector);
    Real n_dot_h = fmax(dot( half_vector, vertex.geometric_normal), Real(0));
    Real out_dot_h = fmax(dot( dir_out, half_vector ), Real(0));

    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool );
    // clearcoat_gloss = std::clamp(clearcoat_gloss, Real(0.001), Real(1));

    Real alpha_g = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
    Real a2 = alpha_g * alpha_g;
    Real D_c = (a2 - 1) / (c_PI  * log(a2) * (1 + (a2 - 1) * h_l.z * h_l.z));

    // return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
    return  D_c * n_dot_h / (4 * out_dot_h);
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Vector3 half_vector = to_world(frame, sample_visible_normals_cc( to_local(frame, dir_in), DISNEY_CC_ROUGHNESS, rnd_param_uv) );
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    // return BSDFSampleRecord{
    //     to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
    //     Real(0) /* eta */, DISNEY_CC_ROUGHNESS};  /* roughness */

    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, DISNEY_CC_ROUGHNESS};  /* roughness */
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
