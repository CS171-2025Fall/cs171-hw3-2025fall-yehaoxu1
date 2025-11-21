#include "rdr/accel.h"

#include "rdr/canary.h"
#include "rdr/interaction.h"
#include "rdr/math_aliases.h"
#include "rdr/platform.h"
#include "rdr/shape.h"

RDR_NAMESPACE_BEGIN

/* ===================================================================== *
 *
 * AABB Implementations
 *
 * ===================================================================== */

bool AABB::isOverlap(const AABB &other) const {
  return ((other.low_bnd[0] >= this->low_bnd[0] &&
              other.low_bnd[0] <= this->upper_bnd[0]) ||
             (this->low_bnd[0] >= other.low_bnd[0] &&
                 this->low_bnd[0] <= other.upper_bnd[0])) &&
         ((other.low_bnd[1] >= this->low_bnd[1] &&
              other.low_bnd[1] <= this->upper_bnd[1]) ||
             (this->low_bnd[1] >= other.low_bnd[1] &&
                 this->low_bnd[1] <= other.upper_bnd[1])) &&
         ((other.low_bnd[2] >= this->low_bnd[2] &&
              other.low_bnd[2] <= this->upper_bnd[2]) ||
             (this->low_bnd[2] >= other.low_bnd[2] &&
                 this->low_bnd[2] <= other.upper_bnd[2]));
}

bool AABB::intersect(const Ray &ray, Float *t_in, Float *t_out) const {
  // Slab method using ray.safe_inverse_direction (precomputed reciprocal
  // direction that handles near-zero components).
  // Compute intersection interval [t_enter, t_exit] across all three axes.
  Float t_enter = -Float_INF;
  Float t_exit  = Float_INF;

  const Vec3f &o    = ray.origin;
  const Vec3f &invd = ray.safe_inverse_direction;

  for (int i = 0; i < 3; ++i) {
    // compute t values where ray crosses the two planes for this axis
    Float t1 = (low_bnd[i] - o[i]) * invd[i];
    Float t2 = (upper_bnd[i] - o[i]) * invd[i];

    Float t_near = std::min(t1, t2);
    Float t_far  = std::max(t1, t2);

    t_enter = std::max(t_enter, t_near);
    t_exit  = std::min(t_exit, t_far);

    // Early exit: if intervals no longer overlap
    if (t_enter > t_exit) return false;
  }

  // Clamp against ray time window
  Float t0 = std::max(t_enter, ray.t_min);
  Float t1 = std::min(t_exit, ray.t_max);

  if (t0 <= t1) {
    if (t_in) *t_in = t0;
    if (t_out) *t_out = t1;
    return true;
  }

  return false;
}

/* ===================================================================== *
 *
 * Accelerator Implementations
 *
 * ===================================================================== */

bool TriangleIntersect(Ray &ray, const uint32_t &triangle_index,
    const ref<TriangleMeshResource> &mesh, SurfaceInteraction &interaction) {
  using InternalScalarType = Double;
  using InternalVecType    = Vec<InternalScalarType, 3>;

  AssertAllValid(ray.direction, ray.origin);
  AssertAllNormalized(ray.direction);

  const auto &vertices = mesh->vertices;
  const Vec3u v_idx(&mesh->v_indices[3 * triangle_index]);
  assert(v_idx.x < mesh->vertices.size());
  assert(v_idx.y < mesh->vertices.size());
  assert(v_idx.z < mesh->vertices.size());

  InternalVecType dir = Cast<InternalScalarType>(ray.direction);
  InternalVecType v0  = Cast<InternalScalarType>(vertices[v_idx[0]]);
  InternalVecType v1  = Cast<InternalScalarType>(vertices[v_idx[1]]);
  InternalVecType v2  = Cast<InternalScalarType>(vertices[v_idx[2]]);

  // TODO(HW3): implement ray-triangle intersection test.
  // You should compute the u, v, t as InternalScalarType
  //
  //   InternalScalarType u = ...;
  //   InternalScalarType v = ...;
  //   InternalScalarType t = ...;
  //
  // And exit early with `return false` if there is no intersection.
  //
  // The intersection points is denoted as:
  // (1 - u - v) * v0 + u * v1 + v * v2 == ray.origin + t * ray.direction
  // where the left side is the barycentric interpolation of the triangle
  // vertices, and the right side is the parametric equation of the ray.
  //
  // You should also make sure that:
  // u >= 0, v >= 0, u + v <= 1, and, ray.t_min <= t <= ray.t_max
  //
  // Useful Functions:
  // You can use @see Cross and @see Dot for determinant calculations.

  // Moller-Trumbore algorithm (stable) implemented in double precision.
  const InternalVecType edge1 = v1 - v0;
  const InternalVecType edge2 = v2 - v0;

  const InternalVecType pvec   = Cross(dir, edge2);
  const InternalScalarType det = Dot(edge1, pvec);

  // If the determinant is near zero, ray lies in plane of triangle or is
  // parallel to it.
  if (abs(det) < InternalScalarType(EPS)) return false;

  const InternalScalarType invDet = InternalScalarType(1) / det;

  const InternalVecType tvec  = Cast<InternalScalarType>(ray.origin) - v0;
  const InternalScalarType uu = Dot(tvec, pvec) * invDet;
  if (uu < InternalScalarType(0) || uu > InternalScalarType(1)) return false;

  const InternalVecType qvec  = Cross(tvec, edge1);
  const InternalScalarType vv = Dot(dir, qvec) * invDet;
  if (vv < InternalScalarType(0) || uu + vv > InternalScalarType(1))
    return false;

  const InternalScalarType tt = Dot(edge2, qvec) * invDet;

  // Check ray parametric range
  if (!(tt >= InternalScalarType(ray.t_min) &&
          tt <= InternalScalarType(ray.t_max)))
    return false;

  // assign results back to expected names
  InternalScalarType u = uu;
  InternalScalarType v = vv;
  InternalScalarType t = tt;

  // We will reach here if there is an intersection

  CalculateTriangleDifferentials(interaction,
      {static_cast<Float>(1 - u - v), static_cast<Float>(u),
          static_cast<Float>(v)},
      mesh, triangle_index);
  AssertNear(interaction.p, ray(t));
  assert(ray.withinTimeRange(t));
  ray.setTimeMax(t);
  return true;
}

void Accel::setTriangleMesh(const ref<TriangleMeshResource> &mesh) {
  // Build the bounding box
  AABB bound(Vec3f(Float_INF, Float_INF, Float_INF),
      Vec3f(Float_MINUS_INF, Float_MINUS_INF, Float_MINUS_INF));
  for (auto &vertex : mesh->vertices) {
    bound.low_bnd   = Min(bound.low_bnd, vertex);
    bound.upper_bnd = Max(bound.upper_bnd, vertex);
  }

  this->mesh  = mesh;   // set the pointer
  this->bound = bound;  // set the bounding box
}

void Accel::build() {}

AABB Accel::getBound() const {
  return bound;
}

bool Accel::intersect(Ray &ray, SurfaceInteraction &interaction) const {
  bool success = false;
  for (int i = 0; i < mesh->v_indices.size() / 3; i++)
    success |= TriangleIntersect(ray, i, mesh, interaction);
  return success;
}

RDR_NAMESPACE_END
