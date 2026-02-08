/**
 * Chiral Anomaly Detection
 *
 * Detects anomalous behavior in chiral oscillator systems.
 * The chiral anomaly occurs when left-handed and right-handed
 * modes become imbalanced in a way that breaks the expected
 * symmetry.
 *
 * Analogous to the Atiyah-Singer index theorem:
 *   index = N_right - N_left = topological charge
 *
 * In the oscillator context:
 * - Left modes: oscillators with negative frequency (counter-clockwise)
 * - Right modes: oscillators with positive frequency (clockwise)
 * - Anomaly: net creation/destruction of one chirality due to topology
 *
 * Physical signatures of anomaly:
 * 1. Spectral asymmetry (eta invariant deviates from zero)
 * 2. Energy imbalance between left and right sectors
 * 3. Non-zero topological charge (winding number changes)
 * 4. Persistent mode count imbalance above noise floor
 *
 * The anomaly monitor tracks these signatures over time and flags
 * sustained anomalous behavior (3+ consecutive anomalous readings).
 */

// ============================================
// TYPES
// ============================================

export interface ChiralAnomalyResult {
  /** Anomaly index: N_right - N_left (integer for exact zero modes) */
  anomalyIndex: number;
  /** Spectral asymmetry eta = sum(sign(lambda_n)) */
  spectralAsymmetry: number;
  /** Conserved topological charge Q = (1/2)(N_right - N_left) */
  topologicalCharge: number;
  /** Whether the anomaly exceeds the detection threshold */
  isAnomalous: boolean;
  /** Total energy in left-propagating modes */
  leftModeEnergy: number;
  /** Total energy in right-propagating modes */
  rightModeEnergy: number;
  /** Asymmetry ratio: |E_R - E_L| / (E_R + E_L) in [0, 1] */
  asymmetryRatio: number;
}

// ============================================
// CHIRAL ANOMALY DETECTION
// ============================================

/**
 * Detect chiral anomaly from left and right mode spectra.
 *
 * Algorithm:
 * 1. Compute energy in each sector: E = sum(amplitude^2)
 * 2. Count active modes above noise floor (amplitude > noiseFloor)
 * 3. Anomaly index = count_right - count_left
 * 4. Spectral asymmetry from sign-weighted amplitude sum
 * 5. Topological charge from zero-mode counting
 * 6. Flag as anomalous if |asymmetryRatio| > threshold
 *
 * @param leftModes - Amplitudes of left-propagating modes
 * @param rightModes - Amplitudes of right-propagating modes
 * @param threshold - Asymmetry ratio threshold for anomaly (default: 0.1)
 * @returns Full anomaly analysis result
 */
export function detectChiralAnomaly(
  leftModes: Float64Array,
  rightModes: Float64Array,
  threshold: number = 0.1
): ChiralAnomalyResult {
  // Noise floor: modes below this amplitude are considered inactive
  const noiseFloor = 1e-6;

  // Compute mode energies: E = sum(a_i^2)
  let leftEnergy = 0;
  let rightEnergy = 0;
  let leftCount = 0;
  let rightCount = 0;
  let leftZeroModes = 0;
  let rightZeroModes = 0;

  for (let i = 0; i < leftModes.length; i++) {
    const amplitude = leftModes[i];
    const energy = amplitude * amplitude;
    leftEnergy += energy;
    if (Math.abs(amplitude) > noiseFloor) {
      leftCount++;
    }
    // Zero modes: very small but non-zero (within noise floor * 10)
    if (Math.abs(amplitude) <= noiseFloor * 10 && Math.abs(amplitude) > noiseFloor * 0.1) {
      leftZeroModes++;
    }
  }

  for (let i = 0; i < rightModes.length; i++) {
    const amplitude = rightModes[i];
    const energy = amplitude * amplitude;
    rightEnergy += energy;
    if (Math.abs(amplitude) > noiseFloor) {
      rightCount++;
    }
    if (Math.abs(amplitude) <= noiseFloor * 10 && Math.abs(amplitude) > noiseFloor * 0.1) {
      rightZeroModes++;
    }
  }

  // Anomaly index: difference in active mode counts
  const anomalyIndex = rightCount - leftCount;

  // Spectral asymmetry: sign-weighted sum over all modes
  // eta = sum(sign(a_n) * |a_n|) for both sectors combined
  const spectralAsymmetry = computeSpectralAsymmetryFromModes(leftModes, rightModes);

  // Topological charge from zero-mode counting
  const topologicalCharge = computeTopologicalCharge(leftZeroModes, rightZeroModes);

  // Energy asymmetry ratio
  const totalEnergy = leftEnergy + rightEnergy;
  const asymmetryRatio = totalEnergy > 1e-15
    ? Math.abs(rightEnergy - leftEnergy) / totalEnergy
    : 0;

  // Anomaly detection: sustained asymmetry above threshold
  const isAnomalous = asymmetryRatio > threshold;

  return {
    anomalyIndex,
    spectralAsymmetry,
    topologicalCharge,
    isAnomalous,
    leftModeEnergy: leftEnergy,
    rightModeEnergy: rightEnergy,
    asymmetryRatio,
  };
}

// ============================================
// SPECTRAL ASYMMETRY (ETA INVARIANT)
// ============================================

/**
 * Compute spectral asymmetry (eta invariant).
 *
 * The eta invariant measures the asymmetry of the spectrum:
 *   eta = sum_n sign(lambda_n) * |lambda_n|^{-s}  evaluated at s = 0
 *
 * At s = 0, this reduces to:
 *   eta = sum_n sign(lambda_n)
 *
 * which counts the difference between positive and negative eigenvalues.
 *
 * For regularization, we use a soft sign function to handle
 * near-zero eigenvalues smoothly.
 *
 * @param eigenvalues - Spectrum of the relevant operator
 * @returns Spectral asymmetry eta
 */
export function computeSpectralAsymmetry(
  eigenvalues: Float64Array
): number {
  let eta = 0;
  const regularization = 1e-8;

  for (let i = 0; i < eigenvalues.length; i++) {
    const lambda = eigenvalues[i];
    // Soft sign: lambda / sqrt(lambda^2 + eps^2) for regularization
    const softSign = lambda / Math.sqrt(lambda * lambda + regularization * regularization);
    eta += softSign;
  }

  return eta;
}

// ============================================
// TOPOLOGICAL CHARGE
// ============================================

/**
 * Compute topological charge from zero-mode spectrum.
 *
 * By the Atiyah-Singer index theorem:
 *   Q = (1/2)(N_right - N_left)
 *
 * where N_right and N_left are the numbers of right-handed
 * and left-handed zero modes of the Dirac operator.
 *
 * The factor of 1/2 converts from the index to the topological
 * charge convention.
 *
 * @param leftZeroModes - Number of left-handed zero modes
 * @param rightZeroModes - Number of right-handed zero modes
 * @returns Topological charge Q (half-integer)
 */
export function computeTopologicalCharge(
  leftZeroModes: number,
  rightZeroModes: number
): number {
  return 0.5 * (rightZeroModes - leftZeroModes);
}

// ============================================
// ANOMALY MONITOR
// ============================================

/**
 * Create a monitor for tracking anomaly development over time.
 *
 * The monitor maintains a sliding window of the most recent anomaly
 * readings and reports the system as anomalous when 3 or more
 * consecutive readings exceed the threshold.
 *
 * This temporal filtering prevents transient fluctuations from
 * triggering false anomaly alerts while still detecting sustained
 * anomalous behavior.
 *
 * @param threshold - Asymmetry ratio threshold for anomaly detection (default: 0.1)
 * @returns Monitor object with update(), getHistory(), isAnomalous(), and reset()
 */
export function createAnomalyMonitor(threshold: number = 0.1): {
  update(leftModes: Float64Array, rightModes: Float64Array): ChiralAnomalyResult;
  getHistory(): ChiralAnomalyResult[];
  isAnomalous(): boolean;
  reset(): void;
} {
  const maxHistory = 16;
  let history: ChiralAnomalyResult[] = [];
  const consecutiveThreshold = 3;

  function update(leftModes: Float64Array, rightModes: Float64Array): ChiralAnomalyResult {
    const result = detectChiralAnomaly(leftModes, rightModes, threshold);

    // Add to history, maintain max size
    history.push(result);
    if (history.length > maxHistory) {
      history = history.slice(history.length - maxHistory);
    }

    return result;
  }

  function getHistory(): ChiralAnomalyResult[] {
    // Return a copy to prevent external mutation
    return history.slice();
  }

  function isAnomalous(): boolean {
    if (history.length < consecutiveThreshold) return false;

    // Check if the last N consecutive readings are anomalous
    let consecutiveCount = 0;
    for (let i = history.length - 1; i >= 0; i--) {
      if (history[i].isAnomalous) {
        consecutiveCount++;
        if (consecutiveCount >= consecutiveThreshold) return true;
      } else {
        break; // Must be consecutive from the end
      }
    }

    return false;
  }

  function reset(): void {
    history = [];
  }

  return { update, getHistory, isAnomalous, reset };
}

// ============================================
// INTERNAL HELPERS
// ============================================

/**
 * Compute spectral asymmetry from separate left and right mode amplitudes.
 *
 * Left modes contribute negatively (sign = -1) and right modes contribute
 * positively (sign = +1) to the spectral asymmetry.
 *
 * eta = sum_n(rightAmplitudes) - sum_n(leftAmplitudes)
 *
 * with regularization for numerical stability.
 */
function computeSpectralAsymmetryFromModes(
  leftModes: Float64Array,
  rightModes: Float64Array
): number {
  let eta = 0;

  // Right modes contribute positively
  for (let i = 0; i < rightModes.length; i++) {
    eta += Math.abs(rightModes[i]);
  }

  // Left modes contribute negatively
  for (let i = 0; i < leftModes.length; i++) {
    eta -= Math.abs(leftModes[i]);
  }

  return eta;
}
