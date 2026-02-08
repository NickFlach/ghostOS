/**
 * Berry Phase and Connection
 *
 * Computes geometric phases for oscillator systems.
 * The Berry phase gamma = oint A dot dl measures the
 * geometric phase accumulated around a loop in parameter space.
 *
 * For Kuramoto oscillators:
 * - Connection A_mu = sum_n (d_mu theta_n) relates to phase gradients
 * - Curvature F = dA captures non-trivial geometry
 * - Berry phase around a loop = integral of curvature over enclosed area
 *
 * Physical interpretation:
 * When oscillator parameters are varied adiabatically around a closed loop,
 * the system acquires a geometric phase that depends only on the geometry
 * of the parameter space, not on the speed of traversal. This Berry phase
 * detects topological structure in the oscillator dynamics.
 *
 * Key formulas:
 *   A_mu = Im(<psi|d_mu|psi>)          (Berry connection)
 *   F_12 = dA_1/dx_2 - dA_2/dx_1       (Berry curvature)
 *   gamma = oint A dot dl               (Berry phase)
 *   gamma = sum_i arg(<psi_i|psi_{i+1}>) (discrete loop formula)
 */

// ============================================
// TYPES
// ============================================

export interface BerryState {
  /** Berry connection 1-form components */
  connection: Float64Array;
  /** Berry curvature at current point (scalar for 2D parameter space) */
  curvature: number;
  /** Total accumulated Berry phase */
  phase: number;
  /** Holonomy of the connection (Berry phase mod 2pi) */
  holonomy: number;
}

// ============================================
// BERRY CONNECTION
// ============================================

/**
 * Compute Berry connection at a point from neighboring states.
 *
 * The Berry connection is the gauge potential:
 *   A_mu = Im(<psi|d_mu|psi>)
 *
 * For oscillators, we approximate using finite differences:
 *   A_mu[i] ~ Im(psi*[i] * (psi_neighbor[i] - psi[i])) / dt
 *
 * where psi*[i] is the complex conjugate (for real oscillators,
 * we use the overlap-based formula).
 *
 * @param state - Current oscillator state vector
 * @param neighbors - States at nearby parameter values (one per parameter direction)
 * @param dt - Parameter step size for finite differences
 * @returns Berry connection components (one per parameter direction)
 */
export function computeBerryConnection(
  state: Float64Array,
  neighbors: Float64Array[],
  dt: number
): Float64Array {
  const numDirections = neighbors.length;
  const n = state.length;
  const connection = new Float64Array(numDirections);

  // Normalize the reference state
  let stateNorm = 0;
  for (let i = 0; i < n; i++) {
    stateNorm += state[i] * state[i];
  }
  stateNorm = Math.sqrt(stateNorm);
  if (stateNorm < 1e-15) return connection;

  for (let mu = 0; mu < numDirections; mu++) {
    const neighbor = neighbors[mu];
    if (neighbor.length !== n) continue;

    // Normalize neighbor
    let neighborNorm = 0;
    for (let i = 0; i < n; i++) {
      neighborNorm += neighbor[i] * neighbor[i];
    }
    neighborNorm = Math.sqrt(neighborNorm);
    if (neighborNorm < 1e-15) continue;

    // Compute <psi|d_mu psi> via finite difference
    // d_mu psi ~ (psi_neighbor - psi) / dt
    // A_mu = Im(<psi|d_mu psi>) / <psi|psi>
    //
    // For real-valued oscillator states, the connection comes from the
    // antisymmetric part of the overlap derivative:
    //   A_mu = sum_i state[i] * (neighbor[i] - state[i]) / (dt * ||state||^2)
    //
    // The imaginary part for real vectors comes from the cross-coupling
    // between position and velocity (phase space Berry connection).
    let overlap = 0;
    let crossTerm = 0;
    for (let i = 0; i < n; i++) {
      const normalizedState = state[i] / stateNorm;
      const normalizedNeighbor = neighbor[i] / neighborNorm;
      const diff = (normalizedNeighbor - normalizedState) / dt;
      overlap += normalizedState * diff;

      // Cross-coupling between even/odd components (phase space structure)
      if (i + 1 < n) {
        crossTerm += normalizedState * (neighbor[i + 1] / neighborNorm - state[i + 1] / stateNorm) / dt;
        crossTerm -= (state[i + 1] / stateNorm) * diff;
      }
    }

    // Berry connection: real part gives the gauge potential for real vectors,
    // cross-term captures the symplectic (phase space) Berry connection
    connection[mu] = overlap + 0.5 * crossTerm;
  }

  return connection;
}

// ============================================
// BERRY CURVATURE
// ============================================

/**
 * Compute Berry curvature from two connection components.
 *
 * The Berry curvature is the field strength of the Berry connection:
 *   F_12 = dA_1/dx_2 - dA_2/dx_1
 *
 * This is computed via finite differences on the connection 1-forms.
 *
 * @param connectionA - Connection components along direction 1 at two points
 *                      [A_1(x, y), A_1(x, y + stepB)]
 * @param connectionB - Connection components along direction 2 at two points
 *                      [A_2(x, y), A_2(x + stepA, y)]
 * @param stepA - Step size in direction 1
 * @param stepB - Step size in direction 2
 * @returns Berry curvature F_12 (scalar)
 */
export function computeBerryCurvature(
  connectionA: Float64Array,
  connectionB: Float64Array,
  stepA: number,
  stepB: number
): number {
  if (connectionA.length < 2 || connectionB.length < 2) return 0;
  if (Math.abs(stepA) < 1e-15 || Math.abs(stepB) < 1e-15) return 0;

  // F_12 = dA_1/dx_2 - dA_2/dx_1
  // dA_1/dx_2 ~ (A_1(x, y+dy) - A_1(x, y)) / dy
  const dA1_dx2 = (connectionA[1] - connectionA[0]) / stepB;

  // dA_2/dx_1 ~ (A_2(x+dx, y) - A_2(x, y)) / dx
  const dA2_dx1 = (connectionB[1] - connectionB[0]) / stepA;

  return dA1_dx2 - dA2_dx1;
}

// ============================================
// BERRY PHASE (LOOP INTEGRAL)
// ============================================

/**
 * Compute Berry phase around a closed loop of states.
 *
 * Uses the discrete Berry phase formula:
 *   gamma = -Im(sum_i ln(<psi_i|psi_{i+1}>))
 *         = -sum_i arg(<psi_i|psi_{i+1}>)
 *
 * For real oscillator states, the overlap <psi_i|psi_{i+1}> is real,
 * so we use the signed angle formula:
 *   phi_i = atan2(|psi_i x psi_{i+1}|, psi_i . psi_{i+1})
 *
 * where "x" denotes a generalized cross product (antisymmetric part).
 *
 * @param states - Array of state vectors forming a closed loop
 *                 (last state connects back to first)
 * @returns Total Berry phase (in radians)
 */
export function computeBerryPhaseLoop(
  states: Float64Array[]
): number {
  const numStates = states.length;
  if (numStates < 3) return 0;

  let totalPhase = 0;

  for (let k = 0; k < numStates; k++) {
    const current = states[k];
    const next = states[(k + 1) % numStates];
    const n = Math.min(current.length, next.length);

    // Compute normalized overlap (dot product)
    let dot = 0;
    let normCurrent = 0;
    let normNext = 0;

    for (let i = 0; i < n; i++) {
      dot += current[i] * next[i];
      normCurrent += current[i] * current[i];
      normNext += next[i] * next[i];
    }

    normCurrent = Math.sqrt(normCurrent);
    normNext = Math.sqrt(normNext);

    if (normCurrent < 1e-15 || normNext < 1e-15) continue;

    const normalizedDot = dot / (normCurrent * normNext);

    // Compute generalized cross product magnitude
    // |a x b| = sqrt(|a|^2|b|^2 - (a.b)^2)
    const crossMagSquared = 1.0 - normalizedDot * normalizedDot;
    const crossMag = Math.sqrt(Math.max(0, crossMagSquared));

    // Compute the signed angle between consecutive states
    // Use the antisymmetric combination to determine sign
    let signedCross = 0;
    for (let i = 0; i < n - 1; i += 2) {
      // Use pairs of components as "2D planes" for cross product sign
      const ax = current[i] / normCurrent;
      const ay = current[i + 1] / normCurrent;
      const bx = next[i] / normNext;
      const by = next[i + 1] / normNext;
      signedCross += ax * by - ay * bx;
    }

    // Berry phase contribution from this segment
    const phi = Math.atan2(
      Math.sign(signedCross) * crossMag,
      normalizedDot
    );
    totalPhase += phi;
  }

  return totalPhase;
}

// ============================================
// BERRY TRACKER
// ============================================

/**
 * Create a BerryState tracker for ongoing phase tracking.
 *
 * Maintains a running Berry phase computation as the system evolves.
 * Each update provides the new state and the parameter displacement,
 * and the tracker accumulates the Berry connection contributions.
 *
 * @returns Object with update(), getState(), and reset() methods
 */
export function createBerryTracker(): {
  update(state: Float64Array, parameterStep: Float64Array): void;
  getState(): BerryState;
  reset(): void;
} {
  let previousState: Float64Array | null = null;
  let accumulatedPhase = 0;
  let currentConnection = new Float64Array(0);
  let currentCurvature = 0;
  let previousConnection: Float64Array | null = null;
  let stepCount = 0;

  function update(state: Float64Array, parameterStep: Float64Array): void {
    const n = state.length;

    if (previousState === null) {
      previousState = new Float64Array(state);
      currentConnection = new Float64Array(parameterStep.length);
      return;
    }

    // Compute connection from finite difference: A ~ <psi|Delta psi> / |Delta lambda|
    const paramNorm = computeNorm(parameterStep);
    if (paramNorm < 1e-15) {
      previousState = new Float64Array(state);
      return;
    }

    // Compute overlap and phase between consecutive states
    let dot = 0;
    let normPrev = 0;
    let normCurr = 0;
    for (let i = 0; i < n; i++) {
      dot += previousState[i] * state[i];
      normPrev += previousState[i] * previousState[i];
      normCurr += state[i] * state[i];
    }
    normPrev = Math.sqrt(normPrev);
    normCurr = Math.sqrt(normCurr);

    if (normPrev > 1e-15 && normCurr > 1e-15) {
      const normalizedDot = dot / (normPrev * normCurr);

      // Phase from overlap
      const clampedDot = Math.max(-1, Math.min(1, normalizedDot));
      const overlapPhase = Math.acos(clampedDot);

      // Sign from cross-product
      let signedCross = 0;
      for (let i = 0; i < n - 1; i += 2) {
        signedCross += (previousState[i] / normPrev) * (state[i + 1] / normCurr)
                     - (previousState[i + 1] / normPrev) * (state[i] / normCurr);
      }

      const signedPhase = Math.sign(signedCross) * overlapPhase;
      accumulatedPhase += signedPhase;

      // Update connection: A_mu = phase / parameterStep_mu
      const newConnection = new Float64Array(parameterStep.length);
      for (let mu = 0; mu < parameterStep.length; mu++) {
        if (Math.abs(parameterStep[mu]) > 1e-15) {
          newConnection[mu] = signedPhase / parameterStep[mu];
        }
      }

      // Estimate curvature from connection change
      if (previousConnection !== null && previousConnection.length === newConnection.length) {
        let curvatureSum = 0;
        for (let mu = 0; mu < newConnection.length; mu++) {
          if (Math.abs(parameterStep[mu]) > 1e-15) {
            curvatureSum += (newConnection[mu] - previousConnection[mu]) / parameterStep[mu];
          }
        }
        currentCurvature = curvatureSum / Math.max(1, newConnection.length);
      }

      previousConnection = newConnection;
      currentConnection = newConnection;
    }

    previousState = new Float64Array(state);
    stepCount++;
  }

  function getState(): BerryState {
    const TWO_PI = 2 * Math.PI;
    // Holonomy is Berry phase mod 2pi, mapped to [-pi, pi]
    let holonomy = accumulatedPhase % TWO_PI;
    if (holonomy > Math.PI) holonomy -= TWO_PI;
    if (holonomy < -Math.PI) holonomy += TWO_PI;

    return {
      connection: new Float64Array(currentConnection),
      curvature: currentCurvature,
      phase: accumulatedPhase,
      holonomy,
    };
  }

  function reset(): void {
    previousState = null;
    previousConnection = null;
    accumulatedPhase = 0;
    currentConnection = new Float64Array(0);
    currentCurvature = 0;
    stepCount = 0;
  }

  return { update, getState, reset };
}

// ============================================
// INTERNAL HELPERS
// ============================================

/**
 * Compute Euclidean norm of a vector.
 */
function computeNorm(v: Float64Array): number {
  let sum = 0;
  for (let i = 0; i < v.length; i++) {
    sum += v[i] * v[i];
  }
  return Math.sqrt(sum);
}
