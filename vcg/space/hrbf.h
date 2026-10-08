/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2026                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
#ifndef VCG_SPACE_HRBF_H
#define VCG_SPACE_HRBF_H

#include <vcg/space/point3.h>
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace vcg {

/** \addtogroup space */
/*@{*/
/**
 * \brief Hermite radial basis function interpolation: an implicit function with given
 * zeros and given gradients there.
 *
 * Fit() takes points p_i and vectors n_i and finds
 *
 *     f(x) = sum_i  alpha_i phi(|x - p_i|) + beta_i . grad phi(|x - p_i|),   phi(r) = r^3
 *
 * with f(p_i) = 0 and grad f(p_i) = n_i: 4N linear conditions on the N scalars alpha_i
 * and N vectors beta_i, solved as one dense system (Macedo, Gois and Velho, "Hermite
 * radial basis functions implicits", Computer Graphics Forum 30(1), 2011; this follows the
 * recipe of R. Vaillant, rodolphe-vaillant.fr/entry/12, with no polynomial term). The zero
 * set of f is then a smooth surface through the points with the given normals; where the
 * n_i point outwards, f is negative inside. Unlike a points-only RBF it needs no off-surface
 * constraints, and the n_i need only be consistently oriented, not of a given length.
 *
 * The kernel is global: every point influences f everywhere, and f grows away from them, so
 * the zero set is meaningful near the points and may open up far from them. Fitting costs
 * O(N^3) and evaluating O(N): it suits a few hundred to a couple of thousand points.
 * Points much closer to each other than to the rest make the system ill conditioned; a
 * small \a regularization (added to the diagonal) trades exact interpolation for stability.
 */
template <class ScalarType = double>
class HRBF
{
public:
  typedef Point3<ScalarType> CoordType;

  /// Fit to \a points with gradients \a normals. Throws std::runtime_error if the system is singular.
  void Fit(const std::vector<CoordType> &points, const std::vector<CoordType> &normals, ScalarType regularization = 0)
  {
    const int n = int(points.size());
    if (n == 0 || normals.size() != points.size())
      throw std::invalid_argument("HRBF::Fit: as many normals as points are needed.");
    centers = points;
    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> Vector;
    Matrix A(4 * n, 4 * n);
    Vector b(4 * n);
    for (int j = 0; j < n; ++j)
    {
      // Row block j: f(p_j) = 0, grad f(p_j) = n_j. Column block i: alpha_i, beta_i.
      for (int i = 0; i < n; ++i)
      {
        ScalarType phi, g[3], h[3][3];
        Kernel(points[j] - points[i], phi, g, h);
        A(4 * j, 4 * i) = phi;
        for (int c = 0; c < 3; ++c)
        {
          A(4 * j, 4 * i + 1 + c) = g[c];       // f picks up beta_i . grad phi
          A(4 * j + 1 + c, 4 * i) = g[c];       // grad f picks up alpha_i grad phi
          for (int d = 0; d < 3; ++d) A(4 * j + 1 + c, 4 * i + 1 + d) = h[c][d];   // and H beta_i
        }
      }
      for (int c = 0; c < 4; ++c) A(4 * j + c, 4 * j + c) += regularization;
      b(4 * j) = 0;
      for (int c = 0; c < 3; ++c) b(4 * j + 1 + c) = normals[j][c];
    }
    const Eigen::FullPivLU<Matrix> lu(A);
    if (!lu.isInvertible())
      throw std::runtime_error("HRBF::Fit: the system is singular (coincident points?).");
    const Vector x = lu.solve(b);
    alpha.resize(n);
    beta.resize(n);
    for (int i = 0; i < n; ++i)
    {
      alpha[i] = x(4 * i);
      beta[i] = CoordType(x(4 * i + 1), x(4 * i + 2), x(4 * i + 3));
    }
  }

  /// f at \a x.
  ScalarType Value(const CoordType &x) const
  {
    ScalarType f = 0;
    for (size_t i = 0; i < centers.size(); ++i)
    {
      const CoordType d = x - centers[i];
      const ScalarType r = d.Norm();
      f += alpha[i] * r * r * r + 3 * r * (beta[i] * d);
    }
    return f;
  }

  /// The gradient of f at \a x.
  CoordType Gradient(const CoordType &x) const
  {
    CoordType grad(0, 0, 0);
    for (size_t i = 0; i < centers.size(); ++i)
    {
      ScalarType phi, g[3], h[3][3];
      Kernel(x - centers[i], phi, g, h);
      for (int c = 0; c < 3; ++c)
      {
        grad[c] += alpha[i] * g[c];
        for (int d = 0; d < 3; ++d) grad[c] += h[c][d] * beta[i][d];
      }
    }
    return grad;
  }

  int Size() const { return int(centers.size()); }

private:
  std::vector<CoordType> centers, beta;
  std::vector<ScalarType> alpha;

  // phi(r) = r^3 at d = x - p, its gradient 3 r d and its Hessian 3 (r I + d d^T / r);
  // all three vanish at d = 0.
  static void Kernel(const CoordType &d, ScalarType &phi, ScalarType g[3], ScalarType h[3][3])
  {
    const ScalarType r = d.Norm();
    phi = r * r * r;
    for (int c = 0; c < 3; ++c)
    {
      g[c] = 3 * r * d[c];
      for (int e = 0; e < 3; ++e) h[c][e] = (r > 0) ? 3 * ((c == e ? r : 0) + d[c] * d[e] / r) : 0;
    }
  }
};
/*@}*/

} // end namespace vcg

#endif // VCG_SPACE_HRBF_H
