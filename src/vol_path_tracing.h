#pragma once

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos = {(x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h};
 
    Ray camera_ray = sample_primary(scene.camera, screen_pos);

    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    auto isect = intersect(scene, camera_ray, ray_diff);
    if (isect) {
        int medium_id = isect->exterior_medium_id;
        Spectrum sigma_a = get_sigma_a(scene.media[medium_id], 
                                isect->position);
        Spectrum sigma_e = get_sigma_s(scene.media[medium_id], 
                                isect->position);

        Real t = distance(camera_ray.org, isect->position);
        Spectrum transmittance = exp(-sigma_a * t);
        Spectrum Le = make_zero_spectrum();

        if (is_light(scene.shapes[isect->shape_id])) {
            Le = emission(*isect, -camera_ray.dir, scene);
        }

        return transmittance * Le;
    }
    return make_zero_spectrum();
}

// Test for adding emission
// Spectrum vol_path_tracing_1(const Scene &scene,
//                                         int x, int y,
//                                         pcg32_state &rng) {
//     int w = scene.camera.width, h = scene.camera.height;
//     Vector2 screen_pos = {(x + next_pcg32_real<Real>(rng)) / w,
//                        (y + next_pcg32_real<Real>(rng)) / h};
 
//     Ray camera_ray = sample_primary(scene.camera, screen_pos);

//     RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
//     auto isect = intersect(scene, camera_ray, ray_diff);

//     int medium_id = scene.camera.medium_id;
//     const Medium &medium = scene.media[medium_id];
    
//     Spectrum sigma_a = get_sigma_a(medium, camera_ray.org);
//     Spectrum Le_media = make_const_spectrum(0.2); // 假设介质有Le_media参数

//     Real t = isect ? distance(camera_ray.org, isect->position) : Real(1e10);

//     Spectrum transmittance = exp(-sigma_a * t);

//     Spectrum volume_contrib = Le_media * (make_const_spectrum(1.f) - transmittance);

//     Spectrum surface_contrib = make_zero_spectrum();
//     if (isect && is_light(scene.shapes[isect->shape_id])) {
//         surface_contrib = transmittance * emission(*isect, -camera_ray.dir, scene);
//     }

//     return volume_contrib;

// }

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos = {(x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h};
 
    Ray camera_ray = sample_primary(scene.camera, screen_pos);

    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    auto isect = intersect(scene, camera_ray, ray_diff);
    int medium_id = scene.camera.medium_id;

    // monochromatic sigma only need the first 
    Real sigma_s = get_sigma_s(scene.media[medium_id], 
                            isect->position)[0];
    Real sigma_a = get_sigma_a(scene.media[medium_id], 
                        isect->position)[0];

    Real sigma_t = sigma_s + sigma_a;
    Real u = next_pcg32_real<Real>(rng);
    Real t = -log(1 - u) / sigma_t;
    if(!isect || t < distance( isect->position, camera_ray.org ) ) {
        Real p_t = exp(-sigma_t * t);
        Real trans_pdf =  p_t* sigma_t;
        Real transmittance  = p_t;
        Vector3 p = camera_ray.org +  camera_ray.dir * t;

        // sample a light
        Real rnd_light = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, rnd_light);

        Real rnd_w = next_pcg32_real<Real>(rng);
        Vector2 rnd_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        auto light = scene.lights[light_id];

        auto point_on_light = sample_point_on_light(light, 
                                            p, rnd_uv, rnd_w, scene);

        Vector3 dir_light = normalize(point_on_light.position - p);
        Real L_s1_pdf = light_pmf(scene, light_id) * pdf_point_on_light(light, point_on_light, p, scene);

        Spectrum L_s1_estimate;

        Ray shadow_ray{p, dir_light,
                       get_shadow_epsilon(scene),
                       (1 - get_shadow_epsilon(scene)) *
                           distance(point_on_light.position, p)};
        if (occluded(scene, shadow_ray)) {
            L_s1_estimate = make_zero_spectrum();
        } else {
            Spectrum Le = emission(light, -dir_light, Real(0), point_on_light, scene);

            Spectrum phase_val = eval(get_phase_function(scene.media[medium_id]), 
                                    -camera_ray.dir, dir_light);

            Real dist = distance(p, point_on_light.position);
            Real cos_light = abs(-dot(dir_light, point_on_light.normal));
            L_s1_estimate = phase_val * Le * exp(-sigma_t * dist) * cos_light / (dist * dist);
        }

        return (transmittance / trans_pdf) * sigma_s * (L_s1_estimate / L_s1_pdf);
    } else {
         Real trans_pdf = exp(-sigma_t * distance(camera_ray.org, isect->position));
        Real transmittance = exp(-sigma_t * distance(camera_ray.org, isect->position));
        Spectrum Le = make_zero_spectrum();
        if (isect && is_light(scene.shapes[isect->shape_id])) {
            Le = emission(*isect, -camera_ray.dir, scene);
        }

        return (transmittance / trans_pdf) * Le;
    }

    return make_zero_spectrum();
}

int update_medium(const PathVertex &isect, Ray &ray, int medium) {
    int new_medium = medium;

    if (isect.interior_medium_id != isect.exterior_medium_id) {
        // At medium transition, update medium
        if (dot(ray.dir, isect.geometric_normal) > 0) {
            new_medium = isect.exterior_medium_id;
        } else {
            new_medium = isect.interior_medium_id;
        }
    }
    return new_medium;
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos = {(x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h};
 
    Ray ray = sample_primary(scene.camera, screen_pos);

    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput = make_const_spectrum(1.);
    Spectrum radiance = make_const_spectrum(0.f);
    u_int bounces = 0;
    while (true) {
        bool scatter = false;
        auto isect = intersect(scene, ray, ray_diff);
        Real transmittance = 1;
        Real trans_pdf = 1;
        if(medium_id != -1) {
            const Medium &medium = scene.media[medium_id];
             // Homogeneous
            Real sigma_s = get_sigma_s(medium, ray.org)[0]; 
            Real sigma_a = get_sigma_a(medium, ray.org)[0];  
            Real sigma_t = sigma_s + sigma_a;     

            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;

            Real p_t = distance(ray.org, isect->position);
            Real exp_p_t = exp(-sigma_t * t);
            Real exp_p_t_hit = exp(-sigma_t * p_t);;

            if (!isect || t < p_t) {
                scatter = true;
                transmittance = exp_p_t;
                trans_pdf = exp_p_t * sigma_t;
            } else {
                transmittance = exp_p_t_hit;
                trans_pdf = exp_p_t_hit;
            }
            ray.org = isect->position + ray.dir * get_intersection_epsilon(scene);

        }

        current_path_throughput *= (transmittance / trans_pdf);
        if(!scatter && isect && is_light(scene.shapes[isect->shape_id])) {
            
            radiance += current_path_throughput * emission(*isect, -ray.dir, scene);

        }
        if(bounces == scene.options.max_depth - 1 && 
                scene.options.max_depth != -1) {
            break;
        }

        if(!scatter && isect && isect->material_id == -1) {
            medium_id = update_medium( *isect, ray, medium_id );
    
            bounces++;
            continue;
        }

        if(scatter) {
            Real sigma_s = get_sigma_s(scene.media[medium_id] , ray.org)[0]; 

            Vector2 rnd_param{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            auto phase_func = get_phase_function( scene.media[medium_id] );
            auto next_dir = sample_phase_function(phase_func, -ray.dir, rnd_param );
            current_path_throughput *= eval(phase_func, -ray.dir, *next_dir) /
                                       pdf_sample_phase(phase_func, -ray.dir, *next_dir) * sigma_s;
            ray.dir = *next_dir;
        } else {
            break;
        }

        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max(current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }
    

    return radiance;
}

Spectrum next_event_estimation(const Scene& scene, Vector3 dir_view, int bounces,
                     Vector3 p, int current_medium, pcg32_state &rng ) {

    Real rnd_light = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, rnd_light);
    const Light &light = scene.lights[light_id];

    Real rnd_w = next_pcg32_real<Real>(rng);
    Vector2 rnd_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            
    auto light_sample = sample_point_on_light(light, p, rnd_uv, rnd_w, scene);
    Vector3 p_prime = light_sample.position;

    Vector3 dir_light = normalize(p_prime - p);
    Vector3 init_p = p;

    Real T_light = (1.f);
    int shadow_medium = current_medium;
    int shadow_bounces = 0;
    Real p_trans_dir = 1.f;  // for multiple importance sampling

    while (true) {
        Ray shadow_ray{p, dir_light,
                       get_shadow_epsilon(scene),
                       (1 - get_shadow_epsilon(scene)) *
                           distance(p_prime, p)};
        auto isect = intersect(scene, shadow_ray);
        Real next_t = distance(p, p_prime);
        if(isect) {
            next_t = distance(p, isect->position);
        }

        if(shadow_medium != -1) {
            auto medium = scene.media[shadow_medium];
            Real sigma_s = get_sigma_s(medium, shadow_ray.org)[0]; 
            Real sigma_a = get_sigma_a(medium, shadow_ray.org)[0];  
            Real sigma_t = sigma_s + sigma_a;     

            Real p_t = exp(-sigma_t * next_t);
            T_light *= p_t;
            p_trans_dir *= p_t; 
        }

        if(!isect) { break; } 
        else {
            if( isect->material_id >= 0 ) {
                return make_zero_spectrum();
            }
            shadow_bounces++;
            if (scene.options.max_depth != -1 && 
                bounces + shadow_bounces + 1 >= scene.options.max_depth) {
                return make_zero_spectrum();
            } 
            shadow_medium = update_medium(*isect, shadow_ray, shadow_medium);  
            p = p + next_t * dir_light;     
        }
        
    }
    
    if( T_light > 0) {      
        Real dist = distance_squared(p, p_prime);
  
        Real G = abs(dot(dir_light, light_sample.normal)) / (dist);
        auto phase = get_phase_function(scene.media[shadow_medium]);
        Spectrum f = eval(phase, dir_view, dir_light);
        Spectrum L = emission(light, -dir_light, Real(0), light_sample, scene);
        Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, light_sample, init_p, scene);
        Spectrum contrib = T_light * G * f * L / pdf_nee;

        Real pdf_phase = pdf_sample_phase(phase, dir_view, dir_light) * G * p_trans_dir;
        Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
        return w * contrib;
    }

    return make_zero_spectrum();
}


// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos = {(x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h};
 
    Ray ray = sample_primary(scene.camera, screen_pos);

    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput = make_const_spectrum(1.);
    Spectrum radiance = make_const_spectrum(0.f);
    u_int bounces = 0;

    // MIS
    bool never_scatter = true;
    Vector3 nee_p_cache = make_zero_spectrum();
    Real multi_trans_pdf = 1;
    Real dir_pdf = 0;

    while (true) {
        bool scatter = false;
        auto isect = intersect(scene, ray, ray_diff);
        Real transmittance = 1;
        Real trans_pdf = 1;
        if(medium_id != -1) {
            const Medium &medium = scene.media[medium_id];
             // Homogeneous
            Real sigma_s = get_sigma_s(medium, ray.org)[0]; 
            Real sigma_a = get_sigma_a(medium, ray.org)[0];  
            Real sigma_t = sigma_s + sigma_a;     

            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;

            Real p_t_hit = distance(ray.org, isect->position);
            Real exp_p_t = exp(-sigma_t * t);
            Real exp_p_t_hit = exp(-sigma_t * p_t_hit);
            if (!isect || t < p_t_hit) {
                scatter = true;
                transmittance = exp_p_t;
                trans_pdf = exp_p_t * sigma_t;
                ray.org += t * ray.dir;
            } else {
                transmittance = exp_p_t_hit;
                trans_pdf = exp_p_t_hit;
                ray.org = isect->position + ray.dir * get_intersection_epsilon(scene);
            }
        } else if (isect) {
            ray.org = isect->position + ray.dir * get_intersection_epsilon(scene);
        } else {
            break;
        }

        multi_trans_pdf *= trans_pdf;
        current_path_throughput *= (transmittance / trans_pdf);
        if(!scatter && isect && is_light(scene.shapes[isect->shape_id])) {
            if(never_scatter) {
                radiance += current_path_throughput * emission(*isect, -ray.dir, scene);
            } else {
                int light = get_area_light_id(scene.shapes[isect->shape_id]);

                PointAndNormal light_point{isect->position, isect->geometric_normal};
                Real pdf_nee =  light_pmf(scene, light) *
                            pdf_point_on_light(scene.lights[light], light_point, nee_p_cache, scene);

                Vector3 dir_light = normalize(isect->position - nee_p_cache);
                Real G = abs(-dot(dir_light, isect->geometric_normal)) / 
                            distance_squared(nee_p_cache, isect->position);
                Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);

                radiance += current_path_throughput * emission(*isect, -ray.dir, scene) * w;
                break;
            }
        }
        if(bounces == scene.options.max_depth - 1 && 
                scene.options.max_depth != -1) {
            break;
        }

        if(!scatter && isect && isect->material_id == -1) {
            medium_id = update_medium( *isect, ray, medium_id );    
            bounces++;
            continue;
        }

        if(scatter) {
            never_scatter = false;

            Spectrum nee = next_event_estimation( scene, -ray.dir, bounces, ray.org, medium_id, rng );

            Real sigma_s = get_sigma_s(scene.media[medium_id] , ray.org)[0]; 

            radiance += current_path_throughput * sigma_s * nee;

            Vector2 rnd_param{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            auto phase_func = get_phase_function( scene.media[medium_id] );
            auto next_dir = sample_phase_function(phase_func, -ray.dir, rnd_param );

            Real phase_pdf = pdf_sample_phase(phase_func, -ray.dir, *next_dir);

            current_path_throughput *= eval(phase_func, -ray.dir, *next_dir) /
                                       phase_pdf * sigma_s;
            
            dir_pdf = phase_pdf;
            nee_p_cache = ray.org;
            multi_trans_pdf = 1;

            ray.dir = *next_dir;
        } else {
            break;
        }

        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max(current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }
    
    return radiance;
}

Spectrum next_event_estimation_bsdf(
    const Scene& scene, Vector3 dir_view, int bounces, 
    const Material &mat, const PathVertex &vertex,
    Vector3 p, int current_medium, pcg32_state &rng ) {

    Real rnd_light = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, rnd_light);
    const Light &light = scene.lights[light_id];

    Real rnd_w = next_pcg32_real<Real>(rng);
    Vector2 rnd_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            
    auto light_sample = sample_point_on_light(light, p, rnd_uv, rnd_w, scene);
    Vector3 p_prime = light_sample.position;

    Vector3 dir_light = normalize(p_prime - p);
    Vector3 init_p = p;

    Real T_light = (1.f);
    int shadow_medium = current_medium;
    int shadow_bounces = 0;
    Real p_trans_dir = 1.f;  // for multiple importance sampling

    Vector3 p_temp  =p;
    while (true) {
        Ray shadow_ray{p_temp, dir_light,
                       get_shadow_epsilon(scene),
                       (1 - get_shadow_epsilon(scene)) *
                           distance(p_prime, p_temp)};
        auto isect = intersect(scene, shadow_ray);
        Real next_t = distance(p_temp, p_prime);
        if(isect) {
            next_t = distance(p_temp, isect->position);
        }

        if(shadow_medium != -1) {
            auto medium = scene.media[shadow_medium];
            Real sigma_s = get_sigma_s(medium, p)[0]; 
            Real sigma_a = get_sigma_a(medium, p)[0];  
            Real sigma_t = sigma_s + sigma_a;     

            Real p_t = exp(-sigma_t * next_t);
            T_light *= p_t;
            p_trans_dir *= p_t; 
        }

        if(!isect) { break; } 
        else {
            if( isect->material_id >= 0 ) {
                return make_zero_spectrum();
            }
            shadow_bounces++;
            if (scene.options.max_depth != -1 && 
                bounces + shadow_bounces + 1 >= scene.options.max_depth) {
                return make_zero_spectrum();
            } 
            shadow_medium = update_medium(*isect, shadow_ray, shadow_medium);  
            p_temp = p_temp + next_t * dir_light;     
        }
        
    }
    
    if( T_light > 0) {      
        Real dist = distance_squared(p_temp, p_prime);
  
        Real G = abs(-dot(dir_light, light_sample.normal)) / (dist);

        Spectrum f =  eval(mat, dir_view, dir_light, vertex, scene.texture_pool);
        Spectrum L = emission(light, -dir_light, Real(0), light_sample, scene);
        Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, light_sample, init_p, scene);
        Spectrum contrib = T_light * G * f * L / pdf_nee;

        Real pdf_bsdf = pdf_sample_bsdf(mat, dir_view, dir_light, vertex, scene.texture_pool) * G * p_trans_dir;
        Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_bsdf * pdf_bsdf);
        return w * contrib;
    }

    return make_zero_spectrum();
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos = {(x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h};
 
    Ray ray = sample_primary(scene.camera, screen_pos);

    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput = make_const_spectrum(1.);
    Spectrum radiance = make_const_spectrum(0.f);
    u_int bounces = 0;

    // MIS
    bool never_scatter = true;
    bool never_bsdf = true;
    Vector3 nee_p_cache = make_zero_spectrum();
    Real multi_trans_pdf = 1;
    Real dir_pdf = 0;
    Real eta_scale = Real(1);


    while (true) {
        bool scatter = false;
        auto isect = intersect(scene, ray, ray_diff);

        Real transmittance = 1;
        Real trans_pdf = 1;
        if(medium_id > -1) {
            const Medium &medium = scene.media[medium_id];
             // Homogeneous
            Real sigma_s = get_sigma_s(medium, ray.org)[0]; 
            Real sigma_a = get_sigma_a(medium, ray.org)[0];  
            Real sigma_t = sigma_s + sigma_a;     

            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;

            Real p_t_hit = isect ? distance(ray.org, isect->position) : INFINITY;
            Real exp_p_t = exp(-sigma_t * t);
            Real exp_p_t_hit = exp(-sigma_t * p_t_hit);
            if (!isect || t < p_t_hit) {
                scatter = true;
                transmittance = exp_p_t;
                trans_pdf = exp_p_t * sigma_t;
                ray.org += t * ray.dir;
                ray.tnear = get_intersection_epsilon(scene);

            } else {
                transmittance = exp_p_t_hit;
                trans_pdf = exp_p_t_hit;
                ray.org = isect->position ;
                ray.tnear = get_intersection_epsilon(scene);
            }
        } else if(isect) {
            ray.org = isect->position;
        }

        multi_trans_pdf *= trans_pdf;
        current_path_throughput *= (transmittance / trans_pdf);

        // reach a surface
        if(!scatter && isect && is_light(scene.shapes[isect->shape_id])) {
            // dont need to count nee contribution if we did not scatter or hit bsdf.
            if(never_scatter && never_bsdf) {
                radiance += current_path_throughput * emission(*isect, -ray.dir, scene);
            } else {
                // light pdf
                int light = get_area_light_id(scene.shapes[isect->shape_id]);

                PointAndNormal light_point{isect->position, isect->geometric_normal};
                Real pdf_nee =  light_pmf(scene, light) *
                            pdf_point_on_light(scene.lights[light], light_point, nee_p_cache, scene);

                Vector3 dir_light = normalize(isect->position - nee_p_cache);
                Real G = abs(-dot(dir_light, isect->geometric_normal)) / 
                            distance_squared(nee_p_cache, isect->position);
                Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);

                radiance += current_path_throughput * emission(*isect, -ray.dir, scene) * w;
                break;
            }
        }
        if(bounces == scene.options.max_depth - 1 && 
                scene.options.max_depth != -1) {
            break;
        }

        // index-matching medium update
        if(!scatter && isect && isect->material_id == -1) {
            medium_id = update_medium( *isect, ray, medium_id );    
            bounces++;
            // ray.org += ray.dir * get_intersection_epsilon(scene);
                ray.tnear = get_intersection_epsilon(scene);

            continue;
        }

        if(scatter) {
            // nee for scattering
            never_scatter = false;

            Spectrum nee = next_event_estimation( scene, -ray.dir, bounces, ray.org, medium_id, rng );

            Real sigma_s = get_sigma_s(scene.media[medium_id] , ray.org)[0]; 

            radiance += current_path_throughput * sigma_s * nee;

            // Sampling phase
            Vector2 rnd_param{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            auto phase_func = get_phase_function( scene.media[medium_id] );
            auto next_dir = sample_phase_function(phase_func, -ray.dir, rnd_param );

            Real phase_pdf = pdf_sample_phase(phase_func, -ray.dir, *next_dir);
            Spectrum f = eval(phase_func, -ray.dir, *next_dir);

            current_path_throughput *=  f / phase_pdf * sigma_s;
            
            
            dir_pdf = phase_pdf;
            nee_p_cache = ray.org;
            multi_trans_pdf = 1;

            ray = Ray{ ray.org, *next_dir,  get_intersection_epsilon(scene), infinity<Real>()};
        } else if(isect){ 
            // hit the brdf surface.
            never_bsdf = false;

            if( &ray == nullptr ) {
                std::cout << "dhusaodsa"<< std::endl;
            }

            const Material &mat = scene.materials[isect->material_id];
            const Vector3 dir_view = -ray.dir;
            Real sigma_s = get_sigma_s(scene.media[medium_id] , ray.org)[0]; 

            // nee for bsdf.
            // std::cout << "main nee running" << std::endl;
            Spectrum nee = next_event_estimation_bsdf(
                scene, -ray.dir, bounces, mat, *isect, ray.org, medium_id, rng);
            // std::cout << "main nee finished" << std::endl;

            radiance += current_path_throughput * nee;

            // Sampling BSDF 
            Vector2 rnd_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real rnd_w = next_pcg32_real<Real>(rng);

            // std::cout << "" << std::endl;

            // std::cout << "main sampling running" << std::endl;
            auto bsdf_sample = sample_bsdf(mat, dir_view, *isect, 
                                            scene.texture_pool, rnd_uv, rnd_w);
            // std::cout << "main sampling finished" << std::endl;

            if (!bsdf_sample) {
                break;
            }
            Vector3 dir_out = bsdf_sample->dir_out;

            // std::cout << "main pdf running" << std::endl;
            Real bsdf_pdf = pdf_sample_bsdf(mat, dir_view, dir_out, *isect, scene.texture_pool);
            // std::cout << "main pdf finished" << std::endl;

            if (bsdf_pdf <= 0) {
                break;
            }


            // std::cout << "main eval running" << std::endl;
            Spectrum f = eval(mat, dir_view, dir_out, *isect, scene.texture_pool);
            // std::cout << "main eval finished" << std::endl;

            current_path_throughput *= f / bsdf_pdf; 

            
            if (bsdf_sample->eta != 0) {
                // refracted, update material
                medium_id = update_medium(*isect, ray, medium_id);
                eta_scale /= (bsdf_sample->eta * bsdf_sample->eta);

            }

            dir_pdf = bsdf_pdf;
            nee_p_cache = ray.org;
            multi_trans_pdf = 1;

            // Update ray

            ray = Ray{ ray.org, dir_out,  get_intersection_epsilon(scene), infinity<Real>()};
        } 

        // Russian
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max((1 / eta_scale) * current_path_throughput), Real(0.95));
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }
    
    return radiance;
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos = {(x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h};
 
    Ray ray = sample_primary(scene.camera, screen_pos);

    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput = make_const_spectrum(1.);
    Spectrum radiance = make_const_spectrum(0.f);
    u_int bounces = 0;

    // MIS
    bool never_scatter = true;
    bool never_bsdf = true;
    Vector3 nee_p_cache = make_zero_spectrum();
    Real multi_trans_pdf = 1;
    Real dir_pdf = 0;
    Real eta_scale = Real(1);


    while (true) {
        bool scatter = false;
        auto isect = intersect(scene, ray, ray_diff);

        Real transmittance = 1;
        Real trans_pdf = 1;
        if(medium_id > -1) {
            const Medium &medium = scene.media[medium_id];
             // Homogeneous
            Real sigma_s = get_sigma_s(medium, ray.org)[0]; 
            Real sigma_a = get_sigma_a(medium, ray.org)[0];  
            Real sigma_t = sigma_s + sigma_a;     

            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;

            Real p_t_hit = isect ? distance(ray.org, isect->position) : INFINITY;
            Real exp_p_t = exp(-sigma_t * t);
            Real exp_p_t_hit = exp(-sigma_t * p_t_hit);
            if (!isect || t < p_t_hit) {
                scatter = true;
                transmittance = exp_p_t;
                trans_pdf = exp_p_t * sigma_t;
                ray.org += t * ray.dir;
                ray.tnear = get_intersection_epsilon(scene);

            } else {
                transmittance = exp_p_t_hit;
                trans_pdf = exp_p_t_hit;
                ray.org = isect->position ;
                ray.tnear = get_intersection_epsilon(scene);
            }
        } else if(isect) {
            ray.org = isect->position;
        }

        multi_trans_pdf *= trans_pdf;
        current_path_throughput *= (transmittance / trans_pdf);

        // reach a surface
        if(!scatter && isect && is_light(scene.shapes[isect->shape_id])) {
            // dont need to count nee contribution if we did not scatter or hit bsdf.
            if(never_scatter && never_bsdf) {
                radiance += current_path_throughput * emission(*isect, -ray.dir, scene);
            } else {
                // light pdf
                int light = get_area_light_id(scene.shapes[isect->shape_id]);

                PointAndNormal light_point{isect->position, isect->geometric_normal};
                Real pdf_nee =  light_pmf(scene, light) *
                            pdf_point_on_light(scene.lights[light], light_point, nee_p_cache, scene);

                Vector3 dir_light = normalize(isect->position - nee_p_cache);
                Real G = abs(-dot(dir_light, isect->geometric_normal)) / 
                            distance_squared(nee_p_cache, isect->position);
                Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);

                radiance += current_path_throughput * emission(*isect, -ray.dir, scene) * w;
                break;
            }
        }
        if(bounces == scene.options.max_depth - 1 && 
                scene.options.max_depth != -1) {
            break;
        }

        // index-matching medium update
        if(!scatter && isect && isect->material_id == -1) {
            medium_id = update_medium( *isect, ray, medium_id );    
            bounces++;
            // ray.org += ray.dir * get_intersection_epsilon(scene);
                ray.tnear = get_intersection_epsilon(scene);

            continue;
        }

        if(scatter) {
            // nee for scattering
            never_scatter = false;

            Spectrum nee = next_event_estimation( scene, -ray.dir, bounces, ray.org, medium_id, rng );

            Real sigma_s = get_sigma_s(scene.media[medium_id] , ray.org)[0]; 

            radiance += current_path_throughput * sigma_s * nee;

            // Sampling phase
            Vector2 rnd_param{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            auto phase_func = get_phase_function( scene.media[medium_id] );
            auto next_dir = sample_phase_function(phase_func, -ray.dir, rnd_param );

            Real phase_pdf = pdf_sample_phase(phase_func, -ray.dir, *next_dir);
            Spectrum f = eval(phase_func, -ray.dir, *next_dir);

            current_path_throughput *=  f / phase_pdf * sigma_s;
            
            
            dir_pdf = phase_pdf;
            nee_p_cache = ray.org;
            multi_trans_pdf = 1;

            ray = Ray{ ray.org, *next_dir,  get_intersection_epsilon(scene), infinity<Real>()};
        } else if(isect){ 
            // hit the brdf surface.
            never_bsdf = false;

            if( &ray == nullptr ) {
                std::cout << "dhusaodsa"<< std::endl;
            }

            const Material &mat = scene.materials[isect->material_id];
            const Vector3 dir_view = -ray.dir;
            Real sigma_s = get_sigma_s(scene.media[medium_id] , ray.org)[0]; 

            // nee for bsdf.
            // std::cout << "main nee running" << std::endl;
            Spectrum nee = next_event_estimation_bsdf(
                scene, -ray.dir, bounces, mat, *isect, ray.org, medium_id, rng);
            // std::cout << "main nee finished" << std::endl;

            radiance += current_path_throughput * nee;

            // Sampling BSDF 
            Vector2 rnd_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real rnd_w = next_pcg32_real<Real>(rng);

            // std::cout << "" << std::endl;

            // std::cout << "main sampling running" << std::endl;
            auto bsdf_sample = sample_bsdf(mat, dir_view, *isect, 
                                            scene.texture_pool, rnd_uv, rnd_w);
            // std::cout << "main sampling finished" << std::endl;

            if (!bsdf_sample) {
                break;
            }
            Vector3 dir_out = bsdf_sample->dir_out;

            // std::cout << "main pdf running" << std::endl;
            Real bsdf_pdf = pdf_sample_bsdf(mat, dir_view, dir_out, *isect, scene.texture_pool);
            // std::cout << "main pdf finished" << std::endl;

            if (bsdf_pdf <= 0) {
                break;
            }


            // std::cout << "main eval running" << std::endl;
            Spectrum f = eval(mat, dir_view, dir_out, *isect, scene.texture_pool);
            // std::cout << "main eval finished" << std::endl;

            current_path_throughput *= f / bsdf_pdf; 

            
            if (bsdf_sample->eta != 0) {
                // refracted, update material
                medium_id = update_medium(*isect, ray, medium_id);
                eta_scale /= (bsdf_sample->eta * bsdf_sample->eta);

            }

            dir_pdf = bsdf_pdf;
            nee_p_cache = ray.org;
            multi_trans_pdf = 1;

            // Update ray

            ray = Ray{ ray.org, dir_out,  get_intersection_epsilon(scene), infinity<Real>()};
        } 

        // Russian
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max((1 / eta_scale) * current_path_throughput), Real(0.95));
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }
    
    return radiance;
}
