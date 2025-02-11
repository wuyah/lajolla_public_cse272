#pragma once
#include <time.h>
#include <chrono>
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

Spectrum next_event_estimation_scatter(const Scene& scene, Vector3 dir_view, int bounces,
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

            Spectrum nee = next_event_estimation_scatter( scene, -ray.dir, bounces, ray.org, medium_id, rng );

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

Spectrum next_event_estimation_blend(
    const Scene& scene, Vector3 dir_view, int bounces, 
    std::optional<Material> mat, const PathVertex &isec,
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
        Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, light_sample, init_p, scene);
        Spectrum L = emission(light, -dir_light, Real(0), light_sample, scene);

        Real pdf_sample = 0;
        Spectrum f = make_const_spectrum(0);

        // If a material passed, means we hit bsdf surface, else, only a 
        if(mat) {
            f =  eval(*mat, dir_view, dir_light, isec, scene.texture_pool);

            pdf_sample = pdf_sample_bsdf(*mat, dir_view, dir_light, isec, scene.texture_pool) * G * p_trans_dir;
        } else {
            auto phase = get_phase_function(scene.media[current_medium]);
            f = eval(phase, dir_view, dir_light);
    
            pdf_sample = pdf_sample_phase(phase, dir_view, dir_light) * G * p_trans_dir;
        }

        Spectrum contrib = T_light * G * f * L / pdf_nee;
        Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_sample * pdf_sample);
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

            // Spectrum nee = 
            //     next_event_estimation_scatter( 
            //         scene, -ray.dir, bounces, 
            //         ray.org, medium_id, rng );

            Spectrum nee = 
                    next_event_estimation_blend(
                        scene, -ray.dir, bounces, std::nullopt, *isect, ray.org, medium_id, rng);
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

            Material mat = scene.materials[isect->material_id];
            const Vector3 dir_view = -ray.dir;
            Real sigma_s = get_sigma_s(scene.media[medium_id] , ray.org)[0]; 

            // nee for bsdf.
            Spectrum nee = next_event_estimation_blend(
                scene, -ray.dir, bounces, std::optional(mat), *isect, ray.org, medium_id, rng);

            radiance += current_path_throughput * nee;

            // Sampling BSDF 
            Vector2 rnd_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real rnd_w = next_pcg32_real<Real>(rng);


            auto bsdf_sample = sample_bsdf(mat, dir_view, *isect, 
                                            scene.texture_pool, rnd_uv, rnd_w);

            if (!bsdf_sample) {
                break;
            }
            Vector3 dir_out = bsdf_sample->dir_out;

            Real bsdf_pdf = pdf_sample_bsdf(mat, dir_view, dir_out, *isect, scene.texture_pool);

            if (bsdf_pdf <= 0) {
                break;
            }


            Spectrum f = eval(mat, dir_view, dir_out, *isect, scene.texture_pool);

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

const int OUTE_M_LOOP_THRESHOLD_MS = 1000;
const int INNER_M_LOOP_THRESHOLD_MS = 1000;
const int OUTER_LOOP_THRESHOLD_MS = 1000;
const int INNER_LOOP_THRESHOLD_MS = 1000;


// final nee with both bsdf and surface,
Spectrum next_event_estimation(
    const Scene& scene, Vector3 dir_view, int bounces, 
    std::optional<Material> mat, PathVertex isec_org,
    Vector3 p_org, int current_medium, pcg32_state &rng 
) {
    // Sampling lights
    Real rnd_light = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, rnd_light);
    const Light &light = scene.lights[light_id];

    Real rnd_w = next_pcg32_real<Real>(rng);
    Vector2 rnd_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            
    auto light_sample = sample_point_on_light(light, p_org, rnd_uv, rnd_w, scene);

    // Prepare
    Vector3 p_prime = light_sample.position;
    Vector3 dir_light = normalize(p_prime - p_org);
    Vector3 p = p_org;

    Spectrum T_light = make_const_spectrum(1.f);
    int shadow_medium = current_medium;
    int shadow_bounces = 0;

    // MIS
    Spectrum p_trans_dir = make_const_spectrum(1);
    Spectrum p_trans_nee = make_const_spectrum(1);

    auto outer_start = std::chrono::high_resolution_clock::now();
    // int i = 0;
    while (true) {
        // i++;
        // assert(i<6);
        // std::cout << p << std::endl;

        auto outer_now = std::chrono::high_resolution_clock::now();
        if (std::chrono::duration_cast<std::chrono::milliseconds>(outer_now - outer_start).count() > OUTER_LOOP_THRESHOLD_MS) {
            std::cout << "Warning: Outer while loop in NEE time is too long."<< std::endl;
            // std::cout << "Loop count" << i << std::endl;
            // break;
            // 可选：break; 或者仅仅输出警告继续循环
        }
        Ray shadow_ray{p, dir_light,
                       get_shadow_epsilon(scene),
                       (1 - get_shadow_epsilon(scene)) *
                           distance(p_prime, p)};
        auto isect = intersect(scene, shadow_ray);
        Real next_t = distance(p, p_prime);
        if(isect) {
            next_t = distance(p, isect->position);
        }

        if(shadow_medium >= 0) {
            assert(shadow_medium >= 0 && shadow_medium <scene.media.size());
            auto medium = scene.media[shadow_medium];
            Spectrum majorant = get_majorant(medium, shadow_ray);

            // sample channel
            Real u = next_pcg32_real<Real>(rng);
            int channel = std::clamp( int(u*3), 0, 2 );
            Real accum_t = 0;
            int iteration = 0;

            auto inner_start = std::chrono::high_resolution_clock::now();

            while (true) {
                auto inner_now = std::chrono::high_resolution_clock::now();
                if (std::chrono::duration_cast<std::chrono::milliseconds>(inner_now - inner_start).count() > INNER_LOOP_THRESHOLD_MS) {
                    std::cout << "Warning: INNER while loop in NEE time is too long."  << std::endl;
                }
        
                // violation or exceed maxium
                if (majorant[channel] <= 0 || 
                    iteration >= scene.options.max_null_collisions) 
                    break;
                
                Real t = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = next_t - accum_t;
                accum_t = min(accum_t + t, next_t);
                if (t < dt) {
                    // scatter
                    p = p + t * dir_light;
                    Spectrum sigma_t = get_sigma_a(medium, p) + get_sigma_s(medium, p);
                    T_light *= exp(-majorant * t) * (majorant - sigma_t) / max(majorant);
                    p_trans_nee *= exp(-majorant * t) * majorant / max(majorant);
                    Spectrum real_prob = sigma_t / majorant;
                    p_trans_dir *= exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant);
                    if (max(T_light) <= 0) {
                        break;
                    }
                } else {
                    // hit surface
                    p = p + dt * dir_light;
                    T_light *= exp(-majorant * dt);
                    p_trans_nee *= exp(-majorant * dt);
                    p_trans_dir *= exp(-majorant * dt);

                    // for fast rendering
                    break;
                }
                iteration++;
            }
        }        

        
        if(!isect) { break; } 
        // occulded
        if( isect->material_id >= 0 ) {
            return make_zero_spectrum();
        }
        shadow_bounces++;
        if (scene.options.max_depth != -1 && 
            bounces + shadow_bounces + 1 >= scene.options.max_depth) {
            return make_zero_spectrum();
        } 
        shadow_medium = update_medium(*isect, shadow_ray, shadow_medium);  
        p = isect->position;     
    }
    
    if( max(T_light) > 0) {    
  
        Real dist_sq = distance_squared(p_org, p_prime);
  
        Real G = abs(dot(dir_light, light_sample.normal)) / (dist_sq);
        Real pdf_nee = light_pmf(scene, light_id) * 
                        pdf_point_on_light(light, light_sample, p_org, scene) * 
                        avg(p_trans_nee);
        Spectrum L = emission(light, -dir_light, Real(0), light_sample, scene);

        Real pdf_sample = 0;
        Spectrum f = make_const_spectrum(0);

        // If a material passed, means we nee a bsdf surface
        if(mat != std::nullopt) {
            // bsdf
            f =  eval(*mat, dir_view, dir_light, isec_org, scene.texture_pool);

            pdf_sample = pdf_sample_bsdf(*mat, dir_view, dir_light, 
                                        isec_org, scene.texture_pool) * 
                        G * avg(p_trans_dir);
        } else {
            // phase
            auto phase = get_phase_function(scene.media[current_medium]);
            f = eval(phase, dir_view, dir_light);
            pdf_sample = pdf_sample_phase(phase, dir_view, dir_light) * G * avg(p_trans_dir);
        }

        Spectrum contrib = T_light * G * f * L / pdf_nee;
        Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_sample * pdf_sample);
        return w * contrib;    
    }

    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
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
    Spectrum multi_trans_pdf = make_const_spectrum(1);
    Spectrum multi_nee_pdf = make_const_spectrum(1);

    Real dir_pdf = 0;
    Real eta_scale = Real(1);

    auto outer_start = std::chrono::high_resolution_clock::now();
    while (true) {
        auto outer_now = std::chrono::high_resolution_clock::now();
        if (std::chrono::duration_cast<std::chrono::milliseconds>(abs(outer_now - outer_start)).count() > OUTE_M_LOOP_THRESHOLD_MS) {
            std::cout << "Warning: Outer while loop in Main time is too long." <<std::endl;
            std::cout << ray.org << std::endl;
            // << "bounces" << bounces << std::endl;
            // break;
            // 可选：break; 或者仅仅输出警告继续循环
        }

        bool scatter = false;
        auto isect = intersect(scene, ray, ray_diff);

        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_dir_pdf = make_const_spectrum(1);
        Spectrum trans_nee_pdf = make_const_spectrum(1);
        Real t_hit = infinity<Real>();
        if(isect) t_hit = distance(ray.org, isect->position);

        if(medium_id >= 0) {
            assert( medium_id < scene.media.size() && medium_id >=0 );
            const Medium &medium = scene.media[medium_id];

            // Heter sampling channel
            Spectrum majorant = get_majorant(medium, ray);
            Real u = next_pcg32_real<Real>(rng);
            int channel = std::clamp( int(u*3), 0, 2 );
            Real accum_t = 0.f;
            int iteration = 0;
            // std::cout << "start channel sampling" << std::endl;
            auto inner_start = std::chrono::high_resolution_clock::now();
            while(true) {
                auto inner_now = std::chrono::high_resolution_clock::now();
                if (std::chrono::duration_cast<std::chrono::milliseconds>(abs(inner_now - inner_start)).count() > INNER_M_LOOP_THRESHOLD_MS) {
                    std::cout << "Warning: INNER while loop in Main time is too long." << std::endl;
                    // break;
                    // 可选：break; 或者仅仅输出警告继续循环
                }
        
                if(majorant[channel] <= 0 || 
                    iteration >= scene.options.max_null_collisions) 
                    break;
                
                Real t = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel] ;
                Real dt = t_hit - accum_t;

                // if( dt < 0.1) break;
                accum_t = min(t_hit, accum_t + t);

                if(t < dt) {
                    ray.org = ray.org + t * ray.dir;

                    Spectrum sigma_t = get_sigma_a(medium, ray.org) + get_sigma_s(medium, ray.org);
                    Spectrum real_prob = sigma_t / majorant;

                    if( next_pcg32_real<Real>(rng) < real_prob[channel] ) {
                        // hit a real particle
                        scatter = true;
                        never_scatter = false;

                        Spectrum exp_t = exp(-majorant * t);

                        transmittance *= exp_t / max(majorant);
                        trans_dir_pdf *= exp_t * majorant * real_prob / max(majorant);

                        break;
                    } else {
                        //  hit fake particle
                        Spectrum sigma_n = majorant - sigma_t;
                        Spectrum exp_t = exp(-majorant * t);
                        transmittance *= exp_t * sigma_n / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant);
                        trans_nee_pdf *= exp(-majorant * t) * majorant / max(majorant);
                    }
                } else {
                    // reach the surface
                    Spectrum exp_dt = exp(-majorant * dt);

                    transmittance *= exp_dt;
                    trans_dir_pdf *= exp_dt;
                    trans_nee_pdf *= exp_dt;

                    ray.org += t_hit * ray.dir;
                    // for fast loading
                    break;
                }
                iteration++;
            }
            // ray.org = ray.org + accum_t * ray.dir;
            multi_trans_pdf *= trans_dir_pdf;
            multi_nee_pdf *= trans_nee_pdf;    
            if (!scatter && !isect) break;  // handle the case where we always hit fake particles until the last iteration

        } else if(isect) {
            ray.org = isect->position;
        } else {
            break;
        }

        current_path_throughput *= (transmittance / avg(trans_dir_pdf));
        
        if(isect) {
            assert( isect->shape_id < scene.shapes.size() && isect->shape_id >=0 );
        }

        if(!scatter && isect  && is_light(scene.shapes[isect->shape_id])) {
            // dont need to count nee contribution if we did not scatter or hit bsdf.
            if(never_scatter && never_bsdf) {
                radiance += current_path_throughput * emission(*isect, -ray.dir, scene);
            } else {
                // light pdf
                int light = get_area_light_id(scene.shapes[isect->shape_id]);

                PointAndNormal light_point{isect->position, isect->geometric_normal};
                Real pdf_nee =  light_pmf(scene, light) *
                            pdf_point_on_light(scene.lights[light], light_point, nee_p_cache, scene)  * avg(multi_nee_pdf);

                Vector3 dir_light = normalize(isect->position - nee_p_cache);
                Real G = abs(-dot(dir_light, isect->geometric_normal)) / 
                            distance_squared(nee_p_cache, isect->position);
                Real dir_pdf_ = dir_pdf * avg(multi_trans_pdf) * G ;
                Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);

                radiance += current_path_throughput * emission(*isect, -ray.dir, scene) * w;
            }
        }
        
        // maxium bounce, exit
        if(bounces == scene.options.max_depth - 1 && 
                scene.options.max_depth != -1) {
            break;
        }
        // index-matching medium update
        if(!scatter && isect && isect->material_id == -1) {
            ray.org += ray.dir * get_intersection_epsilon(scene);

            medium_id = update_medium( *isect, ray, medium_id );    
            if(scatter)
                assert(medium_id != -1 && medium_id < scene.media.size());

            bounces++;
            ray.tnear = get_intersection_epsilon(scene);
            continue;
        }

        if(scatter) { 
            // nee for scattering
            never_scatter = false;
            assert(medium_id != -1 && medium_id < scene.media.size());

            // nee for scatter, material should pass std::nullopt
            Spectrum nee = next_event_estimation(
                scene, -ray.dir, bounces, std::nullopt, *isect, ray.org, medium_id, rng);


            Spectrum sigma_s = get_sigma_s(scene.media[medium_id] , ray.org); 

            radiance += current_path_throughput * sigma_s * nee;

            // Sampling phase
            Vector2 rnd_param{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            auto phase_func = get_phase_function( scene.media[medium_id] );
            auto phase_sample = sample_phase_function(phase_func, -ray.dir, rnd_param );

            Real phase_pdf = pdf_sample_phase(phase_func, -ray.dir, *phase_sample);
            Spectrum f = eval(phase_func, -ray.dir, *phase_sample);

            current_path_throughput *=  f / phase_pdf * sigma_s;
            
            dir_pdf = phase_pdf;
            nee_p_cache = ray.org;
            multi_trans_pdf = make_const_spectrum(1);
            multi_nee_pdf = make_const_spectrum(1);
            ray.org = ray.org;
            ray.dir = *phase_sample;
            ray.tnear = get_intersection_epsilon(scene);

            // ray = Ray{ ray.org, *phase_sample,  get_intersection_epsilon(scene), infinity<Real>()};
        } else if(isect){ // bsdf_branch
            // hit the brdf surface.
            never_bsdf = false;

            assert(isect->material_id < scene.materials.size() && isect->material_id>=0);
            const Material &mat = scene.materials[isect->material_id];
            const Vector3 dir_view = -ray.dir;
            Spectrum sigma_s = get_sigma_s(scene.media[medium_id] , ray.org); 

            // nee for bsdf.

            Spectrum nee = next_event_estimation(
                scene, -ray.dir, bounces, mat, *isect, ray.org, medium_id, rng);

            radiance += current_path_throughput * nee;

            // Sampling BSDF 
            Vector2 rnd_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real rnd_w = next_pcg32_real<Real>(rng);

            auto bsdf_sample = sample_bsdf(mat, dir_view, *isect, 
                                            scene.texture_pool, rnd_uv, rnd_w);

            if (!bsdf_sample) { break; }
            Vector3 dir_out = bsdf_sample->dir_out;

            Real bsdf_pdf = pdf_sample_bsdf(mat, dir_view, dir_out, *isect, scene.texture_pool);

            if (bsdf_pdf <= 0) { break; }

            Spectrum f = eval(mat, dir_view, dir_out, *isect, scene.texture_pool);
            current_path_throughput *= f / bsdf_pdf; 

            
            if (bsdf_sample->eta != 0) {
                // refracted, update material
                medium_id = update_medium(*isect, ray, medium_id);
                eta_scale /= (bsdf_sample->eta * bsdf_sample->eta);
            }

            dir_pdf = bsdf_pdf;
            nee_p_cache = ray.org;
            multi_trans_pdf = make_const_spectrum(1);
            multi_nee_pdf = make_const_spectrum(1);
            // ray.org = isect->position + N * get_intersection_epsilon(scene);

            // Update ray
            ray.org = ray.org;
            ray.dir = dir_out;
            ray.tnear = get_intersection_epsilon(scene);
            // ray = Ray{ ray.org, dir_out,  get_intersection_epsilon(scene), infinity<Real>()};
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
