#include "rdr/integrator.h"

#include <omp.h>

#include "rdr/bsdf.h"
#include "rdr/camera.h"
#include "rdr/canary.h"
#include "rdr/film.h"
#include "rdr/halton.h"
#include "rdr/interaction.h"
#include "rdr/light.h"
#include "rdr/math_aliases.h"
#include "rdr/math_utils.h"
#include "rdr/platform.h"
#include "rdr/properties.h"
#include "rdr/ray.h"
#include "rdr/scene.h"
#include "rdr/sdtree.h"

RDR_NAMESPACE_BEGIN

/* ===================================================================== *
 *
 * Intersection Test Integrator's Implementation
 *
 * ===================================================================== */

void IntersectionTestIntegrator::render(ref<Camera> camera, ref<Scene> scene) {
  // Statistics
  std::atomic<int> cnt = 0;

  const Vec2i &resolution = camera->getFilm()->getResolution();
#pragma omp parallel for schedule(dynamic)
  for (int dx = 0; dx < resolution.x; dx++) {
    ++cnt;
    if (cnt % (resolution.x / 10) == 0)
      Info_("Rendering: {:.02f}%", cnt * 100.0 / resolution.x);
    Sampler sampler;
    for (int dy = 0; dy < resolution.y; dy++) {
      sampler.setPixelIndex2D(Vec2i(dx, dy));
      for (int sample = 0; sample < spp; sample++) {
        const Vec2f pixel_sample = sampler.getPixelSample();
        auto ray =
            camera->generateDifferentialRay(pixel_sample.x, pixel_sample.y);

        // Accumulate radiance
        assert(pixel_sample.x >= dx && pixel_sample.x <= dx + 1);
        assert(pixel_sample.y >= dy && pixel_sample.y <= dy + 1);
        const Vec3f &L = Li(scene, ray, sampler);
        camera->getFilm()->commitSample(pixel_sample, L);
      }
    }
  }
}

Vec3f IntersectionTestIntegrator::Li(
    ref<Scene> scene, DifferentialRay &ray, Sampler &sampler) const {
  SurfaceInteraction interaction;

  for (int depth = 0; depth < max_depth; ++depth) {
    interaction = SurfaceInteraction();
    if (!scene->intersect(ray, interaction)) {
      return Vec3f(0.0);
    }

    interaction.wo   = -ray.direction;
    const BSDF *bsdf = interaction.bsdf;
    if (bsdf == nullptr) return Vec3f(0.0);

    const bool is_ideal_diffuse =
        dynamic_cast<const IdealDiffusion *>(bsdf) != nullptr;
    const bool is_perfect_refraction =
        dynamic_cast<const PerfectRefraction *>(bsdf) != nullptr;

    if (is_perfect_refraction) {
      Float pdf = 0.0f;
      bsdf->sample(interaction, sampler, &pdf);
      Ray refracted_ray = interaction.spawnRay(interaction.wi);
      ray               = refracted_ray;
      continue;
    }

    if (is_ideal_diffuse) {
      return directLighting(scene, interaction);
    }

    // Unsupported surfaces terminate the path in this simple integrator.
    break;
  }

  return Vec3f(0.0);
}

Vec3f IntersectionTestIntegrator::directLighting(
    ref<Scene> scene, SurfaceInteraction &interaction) const {
  Vec3f color(0, 0, 0);
  Vec3f to_light       = point_light_position - interaction.p;
  Float dist2_to_light = Dot(to_light, to_light);
  Float dist_to_light  = std::sqrt(std::max(dist2_to_light, 1e-12f));
  Vec3f light_dir      = to_light / dist_to_light;

  // TODO(HW3): Test for occlusion
  //
  // You should test if there is any intersection between interaction.p and
  // point_light_position using scene->intersect. If so, return an occluded
  // color. (or Vec3f color(0, 0, 0) to be specific)
  //
  // You may find the following variables useful:
  //
  // @see bool Scene::intersect(const Ray &ray, SurfaceInteraction &interaction)
  //    This function tests whether the ray intersects with any geometry in the
  //    scene. And if so, it returns true and fills the interaction with the
  //    intersection information.
  //
  //    You can use iteraction.p to get the intersection position.
  //
  // Spawn a shadow ray toward the point light and clamp its t_max to the
  // distance to the light via spawnRayTo. If any intersection occurs, the
  // light is occluded.
  {
    SurfaceInteraction tmp;
    Ray shadow = interaction.spawnRayTo(point_light_position);
    if (scene->intersect(shadow, tmp)) {
      return color;  // occluded
    }
  }

  // Not occluded, compute the contribution using perfect diffuse diffuse model
  // Perform a quick and dirty check to determine whether the BSDF is ideal
  // diffuse by RTTI
  const BSDF *bsdf      = interaction.bsdf;
  bool is_ideal_diffuse = dynamic_cast<const IdealDiffusion *>(bsdf) != nullptr;

  if (bsdf != nullptr && is_ideal_diffuse) {
    // TODO(HW3): Compute the contribution
    //
    // You can use bsdf->evaluate(interaction) * cos_theta to approximate the
    // albedo. In this homework, we do not need to consider a
    // radiometry-accurate model, so a simple phong-shading-like model is can be
    // used to determine the value of color.

    // The angle between light direction and surface normal
    Float cos_theta =
        std::max(Dot(light_dir, interaction.normal), 0.0f);  // one-sided

    // You should assign the value to color
    interaction.wi    = light_dir;
    const Vec3f &brdf = bsdf->evaluate(interaction);
    // simple inverse-square attenuation
    Float inv_dist2 = 1.0F / std::max(dist2_to_light, 1e-6f);
    color           = point_light_flux * brdf * cos_theta * inv_dist2;
  }

  return color;
}

/* ===================================================================== *
 *
 * Path Integrator's Implementation
 *
 * ===================================================================== */

void PathIntegrator::render(ref<Camera> camera, ref<Scene> scene) {
  // This is left as the next assignment
  UNIMPLEMENTED;
}

Vec3f PathIntegrator::Li(
    ref<Scene> scene, DifferentialRay &ray, Sampler &sampler) const {
  // This is left as the next assignment
  return Vec3f(0.0);
}

Vec3f PathIntegrator::directLighting(
    ref<Scene> scene, SurfaceInteraction &interaction, Sampler &sampler) const {
  // This is left as the next assignment
  return Vec3f(0.0);
}

/* ===================================================================== *
 *
 * New Integrator's Implementation
 *
 * ===================================================================== */

// Instantiate template
// clang-format off
template Vec3f
IncrementalPathIntegrator::Li<Path>(ref<Scene> scene, DifferentialRay &ray, Sampler &sampler) const;
template Vec3f
IncrementalPathIntegrator::Li<PathImmediate>(ref<Scene> scene, DifferentialRay &ray, Sampler &sampler) const;
// clang-format on

// This is exactly a way to separate dec and def
template <typename PathType>
Vec3f IncrementalPathIntegrator::Li(  // NOLINT
    ref<Scene> scene, DifferentialRay &ray, Sampler &sampler) const {
  // This is left as the next assignment
  return Vec3f(0.0);
}

RDR_NAMESPACE_END
