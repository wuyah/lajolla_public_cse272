#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!    
    DisneyDiffuse bsdf_diffuse {
        bsdf.base_color,  
        bsdf.roughness,    
        bsdf.subsurface    
    };
    DisneyClearcoat bsdf_clear_coat {
        bsdf.clearcoat_gloss
    };
    DisneyGlass bsdf_glass {
        bsdf.base_color,
        bsdf.roughness,
        bsdf.anisotropic,
        bsdf.eta
    };
    DisneySheen bsdf_sheen {
        bsdf.base_color,
        bsdf.sheen_tint
    };

    // Only f_g will be calculated ecery time.
    Spectrum f_diff = make_zero_spectrum();
    Spectrum f_metal = make_zero_spectrum();
    Spectrum f_cc = make_zero_spectrum();
    Spectrum f_g = eval(bsdf_glass, dir_in, dir_out, vertex, texture_pool, dir);
    Spectrum f_sh = make_zero_spectrum();

    Vector3 half_vector = normalize( dir_in + dir_out );
    Real n_dot_in = dot( dir_in, vertex.geometric_normal );
    Real n_dot_out = dot( dir_out, vertex.geometric_normal );
    Real h_dot_out = dot( dir_out, half_vector );
    Real h_dot_in = dot( dir_in, half_vector );
    Vector3 h_l = to_local(frame, half_vector);

    Real specular = eval(
        bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_trans = eval( 
        bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool );
    Real specular_tint = eval(
        bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool );
    Real metallic = eval(
        bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.001), Real(1));
    // metallic = std::clamp(metallic, Real(0.001), Real(1));
    // specular_tint = std::clamp(specular_tint, Real(0.01), Real(1));
    // specular_trans = std::clamp(specular_trans, Real(0.01), Real(1));
    // specular = std::clamp(specular, Real(0.01), Real(1));

    Real sheen = eval( 
        bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clear_coat = eval( 
        bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real out_dot_h = fabs(dot(dir_out, half_vector));

    Spectrum base_color = eval(
        bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool );
    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool );
    // anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Real r2 = pow(roughness, 2);
    Real R_0_eta = ( 1 - eta ) * (1 - eta) / ( (1+eta) * (1+eta) );
    Real L = luminance(base_color);
    Spectrum c_tint = L > 0 ? base_color / L : make_const_spectrum(1);
    Spectrum Ks = (1 - specular_tint) + specular_tint * c_tint;
    
    Spectrum C0 = specular * R_0_eta * (1 - metallic) * Ks + metallic * base_color;

    Vector3 F_m = C0 + (1-C0) * pow(1 - out_dot_h, Real(5));

    if(  n_dot_in <=0 )
    { 
        return
           (1-metallic) * specular_trans * f_g;
    } else {
        // Re-calculate the Metal to make it closer to the base color.
        Real aspect = sqrt( 1 - 0.9 * anisotropic );
        Real alpha_x = fmax( 0.0001, r2 / aspect );
        Real alpha_y = fmax( 0.0001, r2 * aspect );

        Real D_m = 1 / (c_PI  * alpha_x * alpha_y *
                    pow( (pow(h_l.x / alpha_x, 2) + pow(h_l.y / alpha_y, 2) + pow(h_l.z, 2)), 2)
                    );
        Real G = smith_masking_disney(to_local(frame, dir_in), Vector2(alpha_x, alpha_y)) *
                smith_masking_disney(to_local(frame, dir_out), Vector2(alpha_x, alpha_y));

        f_metal = F_m * D_m * G / (4 * n_dot_in );
        f_diff = eval(bsdf_diffuse, dir_in, dir_out, vertex, texture_pool, dir);
        f_cc = eval(bsdf_clear_coat, dir_in, dir_out, vertex, texture_pool, dir);
        f_sh = eval(bsdf_sheen, dir_in, dir_out, vertex, texture_pool, dir);
        Spectrum f = (1 - specular_trans) * (1 - metallic) * f_diff +
            (1 - metallic) * sheen * f_sh +
            (1-specular_trans*(1-metallic)) * f_metal +
            0.25 * clear_coat * f_cc +
            (1-metallic) * specular_trans * f_g;
        return f;
    }

}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!


    DisneyDiffuse bsdf_diffuse {
        bsdf.base_color,  
        bsdf.roughness,    
        bsdf.subsurface    
    };
    DisneyMetal bsdf_metal {
        bsdf.base_color,
        bsdf.roughness,
        bsdf.anisotropic
    };
    DisneyClearcoat bsdf_clear_coat {
        bsdf.clearcoat_gloss
    };
    DisneyGlass bsdf_glass {
        bsdf.base_color,
        bsdf.roughness,
        bsdf.anisotropic,
        bsdf.eta
    };

    // wrap weights in a vector4, diff-0, metal-1, glass-2, clearcoat-3.
    Real specular_trans = eval( 
        bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool );

    Real clear_coat = eval( 
        bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(
        bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    // metallic = std::clamp( metallic, Real(0.0001), Real(1) );
    // specular_trans = std::clamp( specular_trans, Real(0.001), Real(1) );
    // clear_coat = std::clamp( clear_coat, Real(0.001), Real(1) );

    Vector4 weights = {0., 0., 0., 0.};
    Real n_dot_in = dot(vertex.geometric_normal, dir_in);
    Real n_dot_out = dot( dir_out, vertex.geometric_normal );

    weights[3] = (1-metallic) * specular_trans;

    // metal weight should be zero when this happen 
    if( !reflect || n_dot_in <=0)
    { 
        return pdf_sample_bsdf( bsdf_glass, dir_in, dir_out, vertex, texture_pool, dir);
    } else {
        weights[0] = (1 - metallic) * (1-specular_trans);
        weights[1] = (1-specular_trans * (1-metallic));
        weights[2] = Real(0.25) * clear_coat;
    }
    Real total = 0;
    for(int i=0; i < 4; i++) {
        total += weights[i];
    }

    assert(total > 0.00001);
    for(int i=0; i < 4; i++) {
        weights[i] /= total;
    }
    // weights = weights / total;
    Real pdf = 0;

    pdf += weights[0] * this->operator()(bsdf_diffuse);
    pdf += weights[1] * this->operator()(bsdf_metal);
    pdf += weights[2] * this->operator()(bsdf_clear_coat);
    
    pdf += weights[3] * this->operator()(bsdf_glass);

    return pdf;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }

    // Evaluate texture parameters
    Real specular_trans = eval(
        bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    // specular_trans = std::clamp(specular_trans, Real(0.01), Real(1));

    Real metallic = eval(
        bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    // metallic = std::clamp(metallic, Real(0.0001), Real(1));

    Real clearcoat = eval(
        bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    // clearcoat = std::clamp(clearcoat, Real(0.001), Real(1));

    DisneyDiffuse bsdf_diffuse {
        bsdf.base_color,  
        bsdf.roughness,    
        bsdf.subsurface    
    };
    DisneyMetal bsdf_metal {
        bsdf.base_color,
        bsdf.roughness,
        bsdf.anisotropic
    };
    DisneyClearcoat bsdf_clear_coat {
        bsdf.clearcoat_gloss
    };
    DisneyGlass bsdf_glass {
        bsdf.base_color,
        bsdf.roughness,
        bsdf.anisotropic,
        bsdf.eta
    };

    Vector4 weights = {0., 0., 0., 0.};
    Real n_dot_in = dot(vertex.geometric_normal, dir_in);
    // Real n_dot_out = dot( dir_out, vertex.geometric_normal );

    if (n_dot_in > 0) {
        weights[0] = (1 - metallic) * (1 - specular_trans);
        weights[1] = (1 - specular_trans * (1 - metallic));
        weights[2] = 0.25 * clearcoat;
    } else {
        return sample_bsdf(bsdf_glass,
            dir_in, vertex, texture_pool, rnd_param_uv,
            rnd_param_w, dir);
    }

    weights[3] = (1 - metallic) * specular_trans;  // Glass weight

    // Compute total weight and check for near-zero scenario
    Real total = weights[0] + weights[1] + weights[2] + weights[3];

    for(int i=0; i < 4; i++) {
        weights[i] /= total;
    }

    assert( weights[0] + weights[1] + weights[2] + weights[3] == 1 );

    Real weightPrefix[5];
    weightPrefix[0] = 0.0;
    for (int i = 0; i < 4; ++i) {
        weightPrefix[i + 1] = weightPrefix[i] + weights[i];
    }

    assert(weightPrefix[4] == 1);
    // Real scaled = std::clamp( rnd_param_w, Real(0), Real(weightPrefix[4]) );
    // Select BSDF component based on random sample rnd_param_w
    int bsdf_selected = -1;
    for (int i = 1; i <= 4; ++i) {
        if (rnd_param_w <= weightPrefix[i]) {
            bsdf_selected = i - 1;
            break;
        }
    }

    assert(bsdf_selected < 4);

    // Rescale the random variable for the selected component

    // Real rescaleW = rnd_param_w;
    // Select appropriate BSDF based on bsdf_selected
    switch (bsdf_selected) {
        case 0:
            return sample_bsdf(bsdf_diffuse,
            dir_in, vertex, texture_pool, rnd_param_uv,
            rnd_param_w, dir);
        case 1:
            return sample_bsdf(bsdf_metal,
            dir_in, vertex, texture_pool, rnd_param_uv,
            rnd_param_w, dir);
        case 2:            
            return sample_bsdf(bsdf_clear_coat,
            dir_in, vertex, texture_pool, rnd_param_uv,
            rnd_param_w, dir);
        case 3:
            assert(weights[3] != 0);
            Real rescaleW = ( rnd_param_w - weightPrefix[3] ) / ( 1- weightPrefix[3] ) ;
            return sample_bsdf(bsdf_glass,
            dir_in, vertex, texture_pool, rnd_param_uv,
            rescaleW, dir);
    }
    // we wont go into this, if so, we return blanlk
    std::cout << "Not Selecting any lobe!" << std::endl;
    return sample_bsdf(bsdf_diffuse,
                        dir_in, vertex, texture_pool, rnd_param_uv,
                        rnd_param_w, dir);
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
