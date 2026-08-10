/*
 * OpenVINS: An Open Platform for Visual-Inertial Research
 * Copyright (C) 2018-2023 Patrick Geneva
 * Copyright (C) 2018-2023 Guoquan Huang
 * Copyright (C) 2018-2023 OpenVINS Contributors
 * Copyright (C) 2018-2019 Kevin Eckenhoff
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef OV_INIT_CERES_JPLQUATLOCAL_H
#define OV_INIT_CERES_JPLQUATLOCAL_H

#include <ceres/ceres.h>

namespace ov_init {

#if CERES_VERSION_MAJOR >= 3 || (CERES_VERSION_MAJOR >= 2 && CERES_VERSION_MINOR >= 2)
// ===================================================================
// Ceres >= 2.2: Manifold API
// ===================================================================

/**
 * @brief JPL quaternion CERES state parameterization (Manifold API).
 *
 * Implements JPL quaternion on-manifold operations for Ceres >= 2.2,
 * where LocalParameterization has been replaced by Manifold.
 *
 * The state is a 4-element quaternion. The tangent space is 3-dimensional
 * (axis-angle perturbation vector). This matches the 3-DOF rotation.
 */
class State_JPLQuatLocal : public ceres::Manifold {
public:
  /**
   * @brief State update: left-multiply by small perturbation quaternion.
   * @f[
   * \bar{q} = norm\Big(\begin{bmatrix} 0.5\mathbf{\theta} \ 1 \end{bmatrix}\Big) \otimes \bar{q}
   * @f]
   * @param x       Current state (4-element JPL quaternion)
   * @param delta   Tangent space perturbation (3-element axis-angle vector)
   * @param x_plus_delta  Output: updated state
   */
  bool Plus(const double *x, const double *delta, double *x_plus_delta) const override;

  /**
   * @brief Jacobian of Plus w.r.t. the tangent space perturbation.
   *
   * Instead of computing dr/dlocal = dr/dglobal * dglobal/dlocal,
   * we directly set dglobal/dlocal = [I_3; 0] so that Ceres can
   * chain its own residual Jacobian correctly.
   */
  bool PlusJacobian(const double *x, double *jacobian) const override;

  /** @brief Ambient space dimension: 4 (quaternion) */
  int AmbientSize() const override { return 4; }

  /** @brief Tangent space dimension: 3 (axis-angle) */
  int TangentSize() const override { return 3; }

  /**
   * @brief Compute tangent space difference: y ⊖ x.
   * Converts the relative rotation y * x^{-1} to an axis-angle vector.
   */
  bool Minus(const double *y, const double *x, double *y_minus_x) const override;

  /**
   * @brief Jacobian of the Minus operation w.r.t. the ambient state x.
   */
  bool MinusJacobian(const double *x, double *jacobian) const override;
};

#else
// ===================================================================
// Ceres < 2.2: LocalParameterization API (backward compatible)
// ===================================================================

/**
 * @brief JPL quaternion CERES local parameterization (legacy API).
 *
 * Identical behavior to the Manifold version above, but using the
 * pre-2.2 LocalParameterization interface for backward compatibility
 * with older Ceres versions (1.x, 2.0, 2.1).
 */
class State_JPLQuatLocal : public ceres::LocalParameterization {
public:
  /**
   * @brief State update (same as Plus in Manifold API).
   */
  bool Plus(const double *x, const double *delta, double *x_plus_delta) const override;

  /**
   * @brief Jacobian of Plus w.r.t. the local parameterization.
   */
  bool ComputeJacobian(const double *x, double *jacobian) const override;

  /** @brief Global/ambient space dimension: 4 */
  int GlobalSize() const override { return 4; }

  /** @brief Local/tangent space dimension: 3 */
  int LocalSize() const override { return 3; }
};

#endif // CERES_VERSION

} // namespace ov_init

#endif // OV_INIT_CERES_JPLQUATLOCAL_H
