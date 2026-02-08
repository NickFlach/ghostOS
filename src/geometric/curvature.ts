/**
 * Curvature Tensors
 *
 * Computes differential geometric quantities on the oscillator manifold.
 * These measure how "curved" the operational state space is.
 *
 * Higher curvature = more complex dynamics, harder to predict
 * Lower curvature = smoother dynamics, more predictable
 *
 * Tensor hierarchy:
 *   Metric g_ij  -->  Christoffel Gamma^k_ij  -->  Riemann R^l_kij
 *     -->  Ricci R_ij  -->  Ricci scalar R
 *
 * Index conventions:
 *   - All tensors are stored in flat Float64Array with explicit indexing
 *   - For an NxN metric, index (i,j) maps to array[i * N + j]
 *   - For rank-3 tensors: array[k][i * N + j] for Gamma^k_ij
 *   - For rank-4 tensors: array[l][k][i * N + j] for R^l_kij
 *
 * Practical limits:
 *   - Designed for dimensions up to 16
 *   - Higher dimensions possible but O(N^4) scaling for Riemann tensor
 */

// ============================================
// CHRISTOFFEL SYMBOLS
// ============================================

/**
 * Compute Christoffel symbols of the second kind from the metric tensor.
 *
 * Formula:
 *   Gamma^k_ij = (1/2) g^kl (d_i g_jl + d_j g_il - d_l g_ij)
 *
 * These are the connection coefficients for the Levi-Civita connection,
 * which encodes how tangent vectors parallel-transport on the manifold.
 *
 * @param metric - g_ij metric tensor (N x N as Float64Array[])
 * @param metricInverse - g^ij inverse metric tensor
 * @param metricDerivatives - d_k g_ij derivatives (N arrays of N x N)
 *        metricDerivatives[k][i][j] = partial_k g_ij
 * @returns Christoffel symbols Gamma^k_ij (N arrays of N x N)
 *          result[k][i * N + j] = Gamma^k_ij
 */
export function computeChristoffelSymbols(
  metric: Float64Array[],
  metricInverse: Float64Array[],
  metricDerivatives: Float64Array[][]
): Float64Array[][] {
  const n = metric.length;
  const christoffel: Float64Array[][] = new Array(n);

  for (let k = 0; k < n; k++) {
    christoffel[k] = new Array(n);
    for (let i = 0; i < n; i++) {
      christoffel[k][i] = new Float64Array(n);
    }
  }

  // Gamma^k_ij = (1/2) g^kl (d_i g_jl + d_j g_il - d_l g_ij)
  for (let k = 0; k < n; k++) {
    for (let i = 0; i < n; i++) {
      for (let j = i; j < n; j++) {
        let sum = 0;
        for (let l = 0; l < n; l++) {
          // d_i g_jl
          const digjl = metricDerivatives[i][j][l];
          // d_j g_il
          const djgil = metricDerivatives[j][i][l];
          // d_l g_ij
          const dlgij = metricDerivatives[l][i][j];

          sum += metricInverse[k][l] * (digjl + djgil - dlgij);
        }
        const value = 0.5 * sum;
        christoffel[k][i][j] = value;
        // Christoffel symbols are symmetric in lower indices
        christoffel[k][j][i] = value;
      }
    }
  }

  return christoffel;
}

// ============================================
// RIEMANN CURVATURE TENSOR
// ============================================

/**
 * Compute Riemann curvature tensor from Christoffel symbols.
 *
 * Full formula:
 *   R^l_kij = d_i Gamma^l_jk - d_j Gamma^l_ik
 *             + Gamma^l_im Gamma^m_jk - Gamma^l_jm Gamma^m_ik
 *
 * For practical computation, we include both the derivative terms
 * and the product terms. The derivative terms require Christoffel
 * symbols at neighboring points (provided via christoffelDerivatives).
 *
 * @param christoffel - Gamma^k_ij at current point
 * @param christoffelDerivatives - d_m Gamma^k_ij (derivatives of Christoffel)
 *        christoffelDerivatives[m][k][i][j] = d_m Gamma^k_ij
 * @returns Riemann tensor R^l_kij
 *          result[l][k][i][j] = R^l_kij
 */
export function computeRiemannTensor(
  christoffel: Float64Array[][],
  christoffelDerivatives: Float64Array[][][]
): Float64Array[][][] {
  const n = christoffel.length;
  const riemann: Float64Array[][][] = new Array(n);

  for (let l = 0; l < n; l++) {
    riemann[l] = new Array(n);
    for (let k = 0; k < n; k++) {
      riemann[l][k] = new Array(n);
      for (let i = 0; i < n; i++) {
        riemann[l][k][i] = new Float64Array(n);
      }
    }
  }

  for (let l = 0; l < n; l++) {
    for (let k = 0; k < n; k++) {
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          // Derivative terms: d_i Gamma^l_jk - d_j Gamma^l_ik
          let value = 0;
          if (christoffelDerivatives.length > 0) {
            value += christoffelDerivatives[i][l][j][k]
                   - christoffelDerivatives[j][l][i][k];
          }

          // Product terms: Gamma^l_im Gamma^m_jk - Gamma^l_jm Gamma^m_ik
          for (let m = 0; m < n; m++) {
            value += christoffel[l][i][m] * christoffel[m][j][k]
                   - christoffel[l][j][m] * christoffel[m][i][k];
          }

          riemann[l][k][i][j] = value;
        }
      }
    }
  }

  return riemann;
}

// ============================================
// RICCI TENSOR
// ============================================

/**
 * Compute Ricci tensor by contracting the Riemann tensor.
 *
 * R_ij = R^k_ikj (contraction on first and third index of Riemann)
 *
 * The Ricci tensor measures the degree to which the geometry determined
 * by the metric differs from flat Euclidean space. It captures volume
 * distortion.
 *
 * @param riemann - Riemann tensor R^l_kij
 * @param dimension - Manifold dimension
 * @returns Ricci tensor R_ij (N x N as Float64Array[])
 */
export function computeRicciTensor(
  riemann: Float64Array[][][],
  dimension: number
): Float64Array[] {
  const n = dimension;
  const ricci: Float64Array[] = new Array(n);

  for (let i = 0; i < n; i++) {
    ricci[i] = new Float64Array(n);
    for (let j = 0; j < n; j++) {
      // R_ij = R^k_ikj = sum_k R^k_ikj
      let sum = 0;
      for (let k = 0; k < n; k++) {
        sum += riemann[k][i][k][j];
      }
      ricci[i][j] = sum;
    }
  }

  return ricci;
}

// ============================================
// RICCI SCALAR
// ============================================

/**
 * Compute Ricci scalar by contracting Ricci tensor with inverse metric.
 *
 * R = g^ij R_ij
 *
 * The Ricci scalar is the simplest curvature invariant. It captures
 * the total curvature at a point in a single number.
 * R > 0: positive curvature (sphere-like)
 * R = 0: flat (Euclidean)
 * R < 0: negative curvature (hyperbolic)
 *
 * @param ricciTensor - R_ij Ricci tensor
 * @param metricInverse - g^ij inverse metric
 * @returns Ricci scalar R
 */
export function computeRicciScalar(
  ricciTensor: Float64Array[],
  metricInverse: Float64Array[]
): number {
  const n = ricciTensor.length;
  let scalar = 0;

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      scalar += metricInverse[i][j] * ricciTensor[i][j];
    }
  }

  return scalar;
}

// ============================================
// METRIC INVERSE (GAUSS-JORDAN)
// ============================================

/**
 * Compute metric inverse via Gauss-Jordan elimination.
 *
 * Given g_ij, find g^ij such that g_ik g^kj = delta_i^j.
 *
 * Uses standard Gauss-Jordan elimination with partial pivoting
 * for numerical stability. Not performance-critical; used for
 * setup and diagnostic computations.
 *
 * @param metric - g_ij metric tensor (N x N)
 * @returns g^ij inverse metric (N x N), or identity if singular
 */
export function invertMetric(metric: Float64Array[]): Float64Array[] {
  const n = metric.length;

  // Build augmented matrix [A | I]
  const augmented: number[][] = new Array(n);
  for (let i = 0; i < n; i++) {
    augmented[i] = new Array(2 * n);
    for (let j = 0; j < n; j++) {
      augmented[i][j] = metric[i][j];
    }
    for (let j = 0; j < n; j++) {
      augmented[i][n + j] = i === j ? 1.0 : 0.0;
    }
  }

  // Forward elimination with partial pivoting
  for (let col = 0; col < n; col++) {
    // Find pivot row
    let maxVal = Math.abs(augmented[col][col]);
    let maxRow = col;
    for (let row = col + 1; row < n; row++) {
      const val = Math.abs(augmented[row][col]);
      if (val > maxVal) {
        maxVal = val;
        maxRow = row;
      }
    }

    // Swap rows if needed
    if (maxRow !== col) {
      const temp = augmented[col];
      augmented[col] = augmented[maxRow];
      augmented[maxRow] = temp;
    }

    // Check for singular matrix
    const pivot = augmented[col][col];
    if (Math.abs(pivot) < 1e-12) {
      // Matrix is singular or nearly singular; return identity as fallback
      const identity: Float64Array[] = new Array(n);
      for (let i = 0; i < n; i++) {
        identity[i] = new Float64Array(n);
        identity[i][i] = 1.0;
      }
      return identity;
    }

    // Scale pivot row
    const invPivot = 1.0 / pivot;
    for (let j = 0; j < 2 * n; j++) {
      augmented[col][j] *= invPivot;
    }

    // Eliminate column entries in all other rows
    for (let row = 0; row < n; row++) {
      if (row === col) continue;
      const factor = augmented[row][col];
      for (let j = 0; j < 2 * n; j++) {
        augmented[row][j] -= factor * augmented[col][j];
      }
    }
  }

  // Extract inverse from augmented matrix
  const inverse: Float64Array[] = new Array(n);
  for (let i = 0; i < n; i++) {
    inverse[i] = new Float64Array(n);
    for (let j = 0; j < n; j++) {
      inverse[i][j] = augmented[i][n + j];
    }
  }

  return inverse;
}

// ============================================
// METRIC DERIVATIVES
// ============================================

/**
 * Compute metric derivatives via central finite differences.
 *
 * Given metrics evaluated at points displaced by +/- epsilon in each
 * coordinate direction, compute:
 *   d_k g_ij ~ (g_ij(x + eps*e_k) - g_ij(x - eps*e_k)) / (2*eps)
 *
 * @param metricAtPoints - Array of metrics at displaced points.
 *        Layout: [g(x+eps*e_0), g(x-eps*e_0), g(x+eps*e_1), g(x-eps*e_1), ...]
 *        So metricAtPoints[2*k] = metric at +epsilon in direction k
 *           metricAtPoints[2*k+1] = metric at -epsilon in direction k
 * @param stepSize - The epsilon displacement
 * @returns d_k g_ij: metricDerivatives[k][i][j]
 */
export function computeMetricDerivatives(
  metricAtPoints: Float64Array[][],
  stepSize: number
): Float64Array[][] {
  // Number of directions = number of pairs / 2
  const numDirections = Math.floor(metricAtPoints.length / 2);
  if (numDirections === 0) return [];

  const n = metricAtPoints[0].length; // Dimension from first metric
  const invTwoEps = 1.0 / (2.0 * stepSize);

  const derivatives: Float64Array[][] = new Array(numDirections);

  for (let k = 0; k < numDirections; k++) {
    const metricPlus = metricAtPoints[2 * k];      // g(x + eps*e_k)
    const metricMinus = metricAtPoints[2 * k + 1];  // g(x - eps*e_k)

    derivatives[k] = new Array(n);
    for (let i = 0; i < n; i++) {
      derivatives[k][i] = new Float64Array(n);
      for (let j = 0; j < n; j++) {
        derivatives[k][i][j] =
          (metricPlus[i][j] - metricMinus[i][j]) * invTwoEps;
      }
    }
  }

  return derivatives;
}

// ============================================
// CONVENIENCE: SCALAR CURVATURE FROM METRIC
// ============================================

/**
 * All-in-one: compute Ricci scalar from metric tensor.
 *
 * This convenience function chains all the computations:
 *   metric --> inverse --> derivatives --> Christoffel --> Riemann --> Ricci tensor --> Ricci scalar
 *
 * For the Riemann tensor, we use the product terms of the Christoffel symbols
 * (the Gamma*Gamma terms) as the primary contribution. The derivative-of-Christoffel
 * terms are computed via finite differences of the Christoffel symbols at
 * neighboring metric points.
 *
 * @param metric - g_ij metric tensor at the point of interest
 * @param metricNeighbors - Metrics at displaced points:
 *        [g(x+eps*e_0), g(x-eps*e_0), g(x+eps*e_1), g(x-eps*e_1), ...]
 * @param epsilon - Displacement size for finite differences
 * @returns Ricci scalar R
 */
export function computeScalarCurvature(
  metric: Float64Array[],
  metricNeighbors: Float64Array[][],
  epsilon: number
): number {
  const n = metric.length;

  // Step 1: Compute inverse metric
  const gInverse = invertMetric(metric);

  // Step 2: Compute metric derivatives via finite differences
  const metricDerivs = computeMetricDerivatives(metricNeighbors, epsilon);

  // If no neighbor data, we can only compute from the metric itself
  if (metricDerivs.length === 0) {
    return 0;
  }

  // Pad derivatives to full dimension if fewer neighbors were provided
  const fullDerivs: Float64Array[][] = new Array(n);
  for (let k = 0; k < n; k++) {
    if (k < metricDerivs.length) {
      fullDerivs[k] = metricDerivs[k];
    } else {
      // Zero derivatives for directions without neighbor data
      fullDerivs[k] = new Array(n);
      for (let i = 0; i < n; i++) {
        fullDerivs[k][i] = new Float64Array(n);
      }
    }
  }

  // Step 3: Compute Christoffel symbols at center point
  const christoffel = computeChristoffelSymbols(metric, gInverse, fullDerivs);

  // Step 4: Compute Christoffel symbols at neighbor points for derivative terms
  // We compute Christoffel at each +/- displaced point, then take finite differences
  const christoffelDerivatives: Float64Array[][][] = new Array(n);
  const numProvidedDirections = Math.min(Math.floor(metricNeighbors.length / 2), n);

  for (let m = 0; m < n; m++) {
    christoffelDerivatives[m] = new Array(n);
    for (let k = 0; k < n; k++) {
      christoffelDerivatives[m][k] = new Array(n);
      for (let i = 0; i < n; i++) {
        christoffelDerivatives[m][k][i] = new Float64Array(n);
      }
    }
  }

  // For each direction with available neighbor data, compute Christoffel derivatives
  for (let m = 0; m < numProvidedDirections; m++) {
    const metricPlus = metricNeighbors[2 * m];
    const metricMinus = metricNeighbors[2 * m + 1];

    const gInvPlus = invertMetric(metricPlus);
    const gInvMinus = invertMetric(metricMinus);

    // We need metric derivatives at the displaced points too.
    // Approximate: reuse the center-point derivatives (first-order accuracy).
    const christoffelPlus = computeChristoffelSymbols(metricPlus, gInvPlus, fullDerivs);
    const christoffelMinus = computeChristoffelSymbols(metricMinus, gInvMinus, fullDerivs);

    const invTwoEps = 1.0 / (2.0 * epsilon);

    for (let k = 0; k < n; k++) {
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          christoffelDerivatives[m][k][i][j] =
            (christoffelPlus[k][i][j] - christoffelMinus[k][i][j]) * invTwoEps;
        }
      }
    }
  }

  // Step 5: Compute Riemann tensor
  const riemann = computeRiemannTensor(christoffel, christoffelDerivatives);

  // Step 6: Contract to Ricci tensor
  const ricci = computeRicciTensor(riemann, n);

  // Step 7: Contract to Ricci scalar
  return computeRicciScalar(ricci, gInverse);
}
