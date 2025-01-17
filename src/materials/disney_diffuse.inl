Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real n_dot_out = fabs(dot(dir_out, frame.n));
    Real n_dot_in = fabs(dot(dir_in, frame.n));
    Real out_dot_h_2 = pow(dot(dir_out, half_vector), Real(2));
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool );
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool );

    Real F_D_90 = 0.5 + 2 * roughness * out_dot_h_2;

    Real F_D_in = (1 + (F_D_90 - 1) * pow( 1-n_dot_in, Real(5)));
    Real F_D_out = (1 + (F_D_90 - 1) * pow( 1-n_dot_out, Real(5)));

    Real f = F_D_out * F_D_in * n_dot_out;

    Spectrum f_diff = f * base_color / c_PI;

    Real F_SS_90 = roughness * out_dot_h_2;
    Real F_SS_in = (1 + (F_SS_90 - 1) * pow( 1-n_dot_in, Real(5)));
    Real F_SS_out = (1 + (F_SS_90 - 1) * pow( 1-n_dot_out, Real(5)));

    Spectrum f_SS = Real(1.25) * base_color / c_PI * n_dot_out *
                    ( F_SS_in * F_SS_out * (1/(n_dot_in + n_dot_out) - 0.5) + 0.5);
    // assert(subsurface == Real(0.5));


    return (1 - subsurface) * f_diff + subsurface * f_SS;
    // return f_SS;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    // Generate a cosine-weighted sample in local space

    // Create and return the sample record
    
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);

    roughness = std::clamp(roughness, Real(0.001), Real(1));

    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, roughness};  /* roughness */
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
