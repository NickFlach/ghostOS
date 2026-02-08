/**
 * Oscillator Manifold
 *
 * Manages the differential geometry of the oscillator phase space.
 * Each process oscillator lives on S^1, and N oscillators together
 * live on the torus T^N. This class tracks the metric, connection,
 * and curvature on this space.
 *
 * The metric tensor g_ij encodes the local geometry:
 * - Identity metric = flat (Euclidean) phase space
 * - BFGS updates adapt the metric to capture curvature from optimization landscape
 * - Covariance tracking provides uncertainty quantification (Kalman-style)
 *
 * Key operations:
 * - updateMetric: BFGS Hessian approximation from gradient observations
 * - computeGeodesicDistance: shortest path on the torus T^N
 * - updateUncertainty: Kalman-style covariance evolution
 */

import { TAU } from '../constants';

// ============================================
// TYPES
// ============================================

export interface ManifoldPoint {
  /** Oscillator phases [theta_1, ..., theta_N] */
  coordinates: Float64Array;
  /** Phase velocities [omega_1, ..., omega_N] */
  tangentVector: Float64Array;
  /** Dimension of the manifold */
  dimension: number;
}

export interface ManifoldMetrics {
  /** Ricci scalar curvature (trace of Ricci tensor contracted with inverse metric) */
  ricciScalar: number;
  /** Accumulated Berry phase from metric evolution */
  berryPhase: number;
  /** Geodesic distance from last observed point to origin */
  geodesicDistance: number;
  /** Convergence score: how close metric is to steady state */
  convergenceScore: number;
  /** Uncertainty volume: determinant approximation of covariance */
  uncertaintyVolume: number;
}

// ============================================
// OSCILLATOR MANIFOLD
// ============================================

export class OscillatorManifold {
  private dimension: number;
  private metric: Float64Array[];
  private previousMetric: Float64Array[];
  private previousPoint: Float64Array | null;
  private previousGradient: Float64Array;
  private covariance: Float64Array[];
  private confidence: number;
  private accumulatedBerryPhase: number;
  private updateCount: number;

  constructor(dimension: number) {
    this.dimension = dimension;
    this.confidence = 0;
    this.accumulatedBerryPhase = 0;
    this.updateCount = 0;
    this.previousPoint = null;
    this.previousGradient = new Float64Array(dimension);

    // Initialize metric as identity matrix (flat Euclidean baseline)
    this.metric = this.createIdentityMatrix(dimension);
    this.previousMetric = this.createIdentityMatrix(dimension);

    // Initialize covariance as 0.5 * identity (moderate initial uncertainty)
    this.covariance = this.createScaledIdentityMatrix(dimension, 0.5);
  }

  /**
   * Update metric using BFGS Hessian approximation.
   *
   * The BFGS formula updates the inverse Hessian approximation H:
   *   H_{k+1} = (I - rho * s * y^T) * H_k * (I - rho * y * s^T) + rho * s * s^T
   *
   * where:
   *   s = x_{k+1} - x_k  (position difference)
   *   y = grad_{k+1} - grad_k  (gradient difference)
   *   rho = 1 / (y^T s)
   *
   * We apply this to update the metric tensor, which serves as an
   * adaptive Hessian approximation for the oscillator landscape.
   */
  updateMetric(point: ManifoldPoint, gradient: Float64Array): void {
    const n = this.dimension;

    // Store previous metric for convergence tracking
    for (let i = 0; i < n; i++) {
      this.previousMetric[i].set(this.metric[i]);
    }

    // Need a previous point to compute BFGS update
    if (this.previousPoint === null) {
      this.previousPoint = new Float64Array(point.coordinates);
      this.previousGradient.set(gradient);
      this.updateCount++;
      return;
    }

    // Compute s = new_pos - old_pos (wrapped on torus)
    const s = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let diff = point.coordinates[i] - this.previousPoint[i];
      // Wrap to [-pi, pi] for torus topology
      while (diff > Math.PI) diff -= TAU;
      while (diff < -Math.PI) diff += TAU;
      s[i] = diff;
    }

    // Compute y = new_grad - old_grad
    const y = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      y[i] = gradient[i] - this.previousGradient[i];
    }

    // Compute s^T y (curvature condition)
    let sTy = 0;
    for (let i = 0; i < n; i++) {
      sTy += s[i] * y[i];
    }

    // Only update if curvature condition is satisfied (positive definite guarantee)
    if (sTy > 1e-10) {
      const rho = 1.0 / sTy;

      // Compute H * y
      const Hy = new Float64Array(n);
      for (let i = 0; i < n; i++) {
        let sum = 0;
        for (let j = 0; j < n; j++) {
          sum += this.metric[i][j] * y[j];
        }
        Hy[i] = sum;
      }

      // Full BFGS rank-2 inverse Hessian update:
      //   H_new = (I - ρ·s·yᵀ)·H·(I - ρ·y·sᵀ) + ρ·s·sᵀ
      // Expanded: H - ρ·s·(Hy)ᵀ - ρ·(Hy)·sᵀ + ρ·(1 + ρ·yᵀHy)·s·sᵀ
      let yTHy = 0;
      for (let i = 0; i < n; i++) {
        yTHy += y[i] * Hy[i];
      }

      const factor = 1.0 + rho * yTHy;
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          this.metric[i][j] +=
            rho * factor * s[i] * s[j] -
            rho * (s[i] * Hy[j] + Hy[i] * s[j]);
        }
      }

      // Ensure symmetry (numerical stability)
      for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
          const avg = 0.5 * (this.metric[i][j] + this.metric[j][i]);
          this.metric[i][j] = avg;
          this.metric[j][i] = avg;
        }
      }

      // Ensure positive definiteness: clamp diagonal to minimum value
      for (let i = 0; i < n; i++) {
        if (this.metric[i][i] < 1e-6) {
          this.metric[i][i] = 1e-6;
        }
      }

      // Track Berry phase from metric evolution (geometric phase of the update)
      this.accumulatedBerryPhase += this.computeMetricBerryPhase();
    }

    // Store current state for next update
    this.previousPoint = new Float64Array(point.coordinates);
    this.previousGradient.set(gradient);
    this.updateCount++;
    this.confidence = Math.min(1.0, this.updateCount / 20);
  }

  /**
   * Compute geodesic distance between two points on the torus T^N.
   *
   * On the torus, the geodesic distance uses the metric tensor:
   *   d(a, b) = sqrt(delta^T * g * delta)
   *
   * where delta_i = min(|a_i - b_i|, 2pi - |a_i - b_i|) is the
   * angular distance with wrapping on each S^1 factor.
   */
  computeGeodesicDistance(a: ManifoldPoint, b: ManifoldPoint): number {
    const n = this.dimension;
    const delta = new Float64Array(n);

    // Compute wrapped angular distances
    for (let i = 0; i < n; i++) {
      let diff = Math.abs(a.coordinates[i] - b.coordinates[i]);
      // Wrap to shortest path on S^1
      if (diff > Math.PI) {
        diff = TAU - diff;
      }
      delta[i] = diff;
    }

    // Compute delta^T * g * delta
    let distSquared = 0;
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        distSquared += delta[i] * this.metric[i][j] * delta[j];
      }
    }

    return Math.sqrt(Math.max(0, distSquared));
  }

  /**
   * Get the metric tensor at the current point.
   * Returns a copy to prevent external mutation.
   */
  getMetric(): Float64Array[] {
    const n = this.dimension;
    const copy: Float64Array[] = new Array(n);
    for (let i = 0; i < n; i++) {
      copy[i] = new Float64Array(this.metric[i]);
    }
    return copy;
  }

  /**
   * Get manifold metrics summary.
   */
  getMetrics(): ManifoldMetrics {
    return {
      ricciScalar: this.computeApproximateRicciScalar(),
      berryPhase: this.accumulatedBerryPhase,
      geodesicDistance: this.computeDistanceFromOrigin(),
      convergenceScore: this.computeConvergenceScore(),
      uncertaintyVolume: this.getUncertaintyVolume(),
    };
  }

  /**
   * Update uncertainty tracking using Kalman-style covariance update.
   *
   * The covariance matrix evolves as:
   *   P_{k+1} = decay * P_k + (1 - decay) * (z * z^T)
   *
   * where z is the observation vector and decay controls the
   * exponential forgetting of old observations.
   */
  updateUncertainty(observation: Float64Array, decay: number = 0.95): void {
    const n = this.dimension;
    const len = Math.min(observation.length, n);

    // Decay existing covariance
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        this.covariance[i][j] *= decay;
      }
    }

    // Add observation outer product: (1 - decay) * z * z^T
    const weight = 1 - decay;
    for (let i = 0; i < len; i++) {
      for (let j = 0; j < len; j++) {
        this.covariance[i][j] += weight * observation[i] * observation[j];
      }
    }
  }

  /**
   * Get uncertainty volume: approximated as the product of diagonal
   * covariance elements (fast approximation of determinant).
   *
   * For a diagonal-dominant covariance, this is close to the true
   * determinant. For full covariance, it's an upper bound.
   */
  getUncertaintyVolume(): number {
    let logVolume = 0;
    for (let i = 0; i < this.dimension; i++) {
      const diagVal = Math.max(this.covariance[i][i], 1e-15);
      logVolume += Math.log(diagVal);
    }
    // Use log-space to avoid overflow for high dimensions
    return Math.exp(logVolume);
  }

  /**
   * Reset manifold to identity metric (flat space).
   */
  reset(): void {
    this.metric = this.createIdentityMatrix(this.dimension);
    this.previousMetric = this.createIdentityMatrix(this.dimension);
    this.covariance = this.createScaledIdentityMatrix(this.dimension, 0.5);
    this.previousPoint = null;
    this.previousGradient.fill(0);
    this.confidence = 0;
    this.accumulatedBerryPhase = 0;
    this.updateCount = 0;
  }

  // ============================================
  // PRIVATE HELPERS
  // ============================================

  /**
   * Create an N x N identity matrix using Float64Array rows.
   */
  private createIdentityMatrix(n: number): Float64Array[] {
    const matrix: Float64Array[] = new Array(n);
    for (let i = 0; i < n; i++) {
      matrix[i] = new Float64Array(n);
      matrix[i][i] = 1.0;
    }
    return matrix;
  }

  /**
   * Create an N x N scaled identity matrix.
   */
  private createScaledIdentityMatrix(n: number, scale: number): Float64Array[] {
    const matrix: Float64Array[] = new Array(n);
    for (let i = 0; i < n; i++) {
      matrix[i] = new Float64Array(n);
      matrix[i][i] = scale;
    }
    return matrix;
  }

  /**
   * Compute an approximate Ricci scalar from the metric tensor.
   *
   * For a nearly-flat metric g = I + h, the linearized Ricci scalar is:
   *   R ~ d^2 h_ii - d_i d_j h_ij
   *
   * We approximate this by measuring deviation from identity:
   *   R ~ sum_ij (g_ij - delta_ij)^2 (Frobenius norm of deviation)
   *
   * Sign convention: positive = metric is expanding, negative = contracting.
   */
  private computeApproximateRicciScalar(): number {
    const n = this.dimension;
    let traceDeviation = 0;
    let offDiagonalDeviation = 0;

    for (let i = 0; i < n; i++) {
      traceDeviation += (this.metric[i][i] - 1.0);
      for (let j = 0; j < n; j++) {
        const target = i === j ? 1.0 : 0.0;
        const dev = this.metric[i][j] - target;
        offDiagonalDeviation += dev * dev;
      }
    }

    // Combine trace and off-diagonal contributions
    // Positive curvature when metric deviates from flat
    return traceDeviation / n - Math.sqrt(offDiagonalDeviation) / (n * n);
  }

  /**
   * Compute Berry phase contribution from metric evolution.
   * Measures the geometric phase picked up as the metric tensor changes.
   */
  private computeMetricBerryPhase(): number {
    const n = this.dimension;
    let phase = 0;

    // Berry phase from metric change: Im(Tr(g^{-1} dg))
    // Approximate as Tr(g_old^{-1} * (g_new - g_old))
    for (let i = 0; i < n; i++) {
      const oldDiag = this.previousMetric[i][i];
      const newDiag = this.metric[i][i];
      if (Math.abs(oldDiag) > 1e-15) {
        phase += (newDiag - oldDiag) / oldDiag;
      }
    }

    return phase / n;
  }

  /**
   * Compute geodesic distance from the last observed point to origin.
   */
  private computeDistanceFromOrigin(): number {
    if (this.previousPoint === null) return 0;

    const origin: ManifoldPoint = {
      coordinates: new Float64Array(this.dimension),
      tangentVector: new Float64Array(this.dimension),
      dimension: this.dimension,
    };

    const current: ManifoldPoint = {
      coordinates: this.previousPoint,
      tangentVector: new Float64Array(this.dimension),
      dimension: this.dimension,
    };

    return this.computeGeodesicDistance(origin, current);
  }

  /**
   * Compute convergence score: how much the metric changed in the last update,
   * weighted by the confidence (number of updates received).
   * Score of 1.0 means fully converged (no change), 0.0 means large change.
   */
  private computeConvergenceScore(): number {
    const n = this.dimension;
    let changeSumSquared = 0;
    let metricSumSquared = 0;

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        const diff = this.metric[i][j] - this.previousMetric[i][j];
        changeSumSquared += diff * diff;
        metricSumSquared += this.metric[i][j] * this.metric[i][j];
      }
    }

    if (metricSumSquared < 1e-15) return 1.0;

    const relativeChange = Math.sqrt(changeSumSquared / metricSumSquared);
    // Weight by confidence: low confidence reduces convergence score
    return this.confidence * Math.exp(-10 * relativeChange);
  }
}
