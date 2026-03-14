/**
 * ghostmagic Safety Envelope
 *
 * Boundary enforcement and stability monitoring.
 * Prevents runaway divergence and collapse.
 *
 * Safety Thresholds:
 * - ||E|| > 2.0 → Normalize (prevent runaway accumulation)
 * - ||E|| < 0.1 → Boost λ (increase emergence)
 * - ρ < 0.1 → Reset phase (restore resonance)
 * - C < 0.2 → Inject noise (break degeneracy)
 */

import {
  COHERENCE_MIN,
  COHERENCE_MAX,
  EMERGENCE_NORM_MIN,
  EMERGENCE_NORM_MAX,
  DAMPING_MAX,
  CHIRAL_STABILITY,
  TOPOLOGICAL_PROTECTION_FACTOR,
} from '../constants';
import type { ChiralState } from '../chiral';

// ============================================
// TYPES
// ============================================

export interface SafetyConfig {
  coherenceMin?: number;
  coherenceMax?: number;
  emergenceMin?: number;
  emergenceMax?: number;
  energyMax?: number;
  divergenceThreshold?: number;
}

export type SafetyAction =
  | 'nominal'
  | 'reduce-coupling'
  | 'boost-coherence'
  | 'inject-noise'
  | 'normalize'
  | 'reset-phase'
  | 'emergency-halt'
  | 'chiral-rebalance'      // Adjust chiral parameters
  | 'flip-handedness'        // Switch preferred chirality
  | 'topological-reset';     // Reset to topologically protected state

export interface SafetyStatus {
  status: 'nominal' | 'warning' | 'critical';
  intervention: boolean;
  action: SafetyAction;
  metrics: SafetyMetrics;
  violations: string[];
}

export interface SafetyMetrics {
  coherence: number;
  emergenceNorm: number;
  energy: number;
  divergenceRate: number;
  stabilityScore: number;
  /** Chiral metrics (if available) */
  chiral?: {
    velocity: number;
    asymmetry: number;
    windingNumber: number;
    regime: 'stable' | 'transitional' | 'unstable';
  };
}

export interface SystemState {
  coherence: number;
  emergenceNorm: number;
  energy: number;
  previousEnergy?: number;
  phase?: number;
  /** Chiral state from resonance engine */
  chiral?: ChiralState;
}

// ============================================
// SAFETY ENVELOPE
// ============================================

export class SafetyEnvelope {
  private coherenceMin: number;
  private coherenceMax: number;
  private emergenceMin: number;
  private emergenceMax: number;
  private energyMax: number;
  private divergenceThreshold: number;
  private history: number[] = [];
  private maxHistoryLength = 100;

  constructor(config: SafetyConfig = {}) {
    this.coherenceMin = config.coherenceMin ?? COHERENCE_MIN;
    this.coherenceMax = config.coherenceMax ?? COHERENCE_MAX;
    this.emergenceMin = config.emergenceMin ?? EMERGENCE_NORM_MIN;
    this.emergenceMax = config.emergenceMax ?? EMERGENCE_NORM_MAX;
    this.energyMax = config.energyMax ?? 100;
    this.divergenceThreshold = config.divergenceThreshold ?? 0.1;
  }

  /**
   * Check system state against safety envelope
   */
  check(state: SystemState): SafetyStatus {
    const violations: string[] = [];
    let action: SafetyAction = 'nominal';
    let status: 'nominal' | 'warning' | 'critical' = 'nominal';

    // Track energy history for divergence detection
    this.history.push(state.energy);
    if (this.history.length > this.maxHistoryLength) {
      this.history.shift();
    }

    const divergenceRate = this.calculateDivergenceRate();

    // Check emergence norm
    if (state.emergenceNorm > this.emergenceMax) {
      violations.push(`Emergence norm ${state.emergenceNorm.toFixed(3)} exceeds max ${this.emergenceMax}`);
      action = 'normalize';
      status = 'warning';
    } else if (state.emergenceNorm < this.emergenceMin) {
      violations.push(`Emergence norm ${state.emergenceNorm.toFixed(3)} below min ${this.emergenceMin}`);
      action = 'boost-coherence';
      status = 'warning';
    }

    // Check coherence
    if (state.coherence > this.coherenceMax) {
      violations.push(`Coherence ${state.coherence.toFixed(3)} exceeds max ${this.coherenceMax} (over-synchronization)`);
      action = 'reduce-coupling';
      status = 'warning';
    } else if (state.coherence < this.coherenceMin) {
      violations.push(`Coherence ${state.coherence.toFixed(3)} below min ${this.coherenceMin} (decoherence)`);
      action = action === 'nominal' ? 'inject-noise' : action;
      status = 'warning';
    }

    // Check energy bounds
    if (state.energy > this.energyMax) {
      violations.push(`Energy ${state.energy.toFixed(3)} exceeds max ${this.energyMax} (runaway)`);
      action = 'normalize';
      status = 'critical';
    }

    // Check divergence rate
    if (divergenceRate > this.divergenceThreshold) {
      violations.push(`Divergence rate ${divergenceRate.toFixed(4)} exceeds threshold ${this.divergenceThreshold}`);
      action = 'reduce-coupling';
      status = status === 'critical' ? 'critical' : 'warning';
    }

    // Check chiral stability if available
    if (state.chiral) {
      const chiralViolations = this.checkChiralStability(state.chiral);
      violations.push(...chiralViolations.violations);

      if (chiralViolations.action !== 'nominal' && action === 'nominal') {
        action = chiralViolations.action;
        status = chiralViolations.status;
      }
    }

    // Critical condition: multiple violations
    if (violations.length >= 3) {
      action = 'emergency-halt';
      status = 'critical';
    }

    const stabilityScore = this.calculateStabilityScore(state, divergenceRate);

    // Build metrics with optional chiral data
    const metrics: SafetyMetrics = {
      coherence: state.coherence,
      emergenceNorm: state.emergenceNorm,
      energy: state.energy,
      divergenceRate,
      stabilityScore,
    };

    if (state.chiral) {
      metrics.chiral = {
        velocity: state.chiral.chiralVelocity,
        asymmetry: state.chiral.asymmetry,
        windingNumber: state.chiral.windingNumber,
        regime: state.chiral.stabilityRegime,
      };
    }

    return {
      status,
      intervention: action !== 'nominal',
      action,
      metrics,
      violations,
    };
  }

  /**
   * Check chiral stability metrics
   * Based on non-reciprocal wave stability research
   */
  private checkChiralStability(chiral: ChiralState): {
    violations: string[];
    action: SafetyAction;
    status: 'nominal' | 'warning' | 'critical';
  } {
    const violations: string[] = [];
    let action: SafetyAction = 'nominal';
    let status: 'nominal' | 'warning' | 'critical' = 'nominal';

    // Check chiral velocity regime
    // |c| = |η/Γ| determines stability
    const absVelocity = Math.abs(chiral.chiralVelocity);

    if (chiral.stabilityRegime === 'unstable') {
      violations.push(
        `Chiral velocity |c| = ${absVelocity.toFixed(3)} in unstable regime (>${CHIRAL_STABILITY.transitional.maxRatio})`
      );
      action = 'chiral-rebalance';
      status = 'warning';
    } else if (chiral.stabilityRegime === 'transitional') {
      violations.push(
        `Chiral velocity |c| = ${absVelocity.toFixed(3)} in transitional regime`
      );
      // Only warn, don't change action unless it's nominal
      if (action === 'nominal') {
        status = 'warning';
      }
    }

    // Check asymmetry - too high means one chirality dominates excessively
    if (chiral.asymmetry > 0.9) {
      violations.push(
        `Chiral asymmetry ${chiral.asymmetry.toFixed(3)} exceeds 0.9 (over-polarization)`
      );
      action = action === 'nominal' ? 'flip-handedness' : action;
      status = 'warning';
    } else if (chiral.asymmetry < 0.1 && chiral.handedness === 'achiral') {
      // Too symmetric - no preferred direction, may lose stability benefits
      violations.push(
        `System achiral (asymmetry ${chiral.asymmetry.toFixed(3)}) - stability benefits reduced`
      );
    }

    // Check winding number stability
    // Large winding numbers indicate topological complexity
    if (Math.abs(chiral.windingNumber) > 3) {
      violations.push(
        `Topological winding number ${chiral.windingNumber} exceeds safe range`
      );
      action = action === 'nominal' ? 'topological-reset' : action;
      status = status === 'nominal' ? 'warning' : status;
    }

    // Check spin polarization (CISS metric)
    if (Math.abs(chiral.spinPolarization) > 0.95) {
      violations.push(
        `Spin polarization ${chiral.spinPolarization.toFixed(3)} near saturation`
      );
    }

    return { violations, action, status };
  }

  /**
   * Calculate rate of energy divergence
   */
  private calculateDivergenceRate(): number {
    if (this.history.length < 2) return 0;

    const recent = this.history.slice(-10);
    let totalChange = 0;

    for (let i = 1; i < recent.length; i++) {
      totalChange += Math.abs(recent[i] - recent[i - 1]);
    }

    return totalChange / recent.length;
  }

  /**
   * Calculate overall stability score [0, 1]
   * Now includes chiral stability factors
   */
  private calculateStabilityScore(state: SystemState, divergenceRate: number): number {
    let score = 1.0;

    // Penalize coherence outside optimal band
    const coherenceMid = (this.coherenceMin + this.coherenceMax) / 2;
    const coherenceRange = this.coherenceMax - this.coherenceMin;
    const coherenceDeviation = Math.abs(state.coherence - coherenceMid) / (coherenceRange / 2);
    score *= Math.max(0, 1 - coherenceDeviation * 0.5);

    // Penalize emergence outside healthy range
    const emergenceMid = (this.emergenceMin + this.emergenceMax) / 2;
    const emergenceRange = this.emergenceMax - this.emergenceMin;
    const emergenceDeviation = Math.abs(state.emergenceNorm - emergenceMid) / (emergenceRange / 2);
    score *= Math.max(0, 1 - emergenceDeviation * 0.3);

    // Penalize high divergence rate
    score *= Math.max(0, 1 - divergenceRate / this.divergenceThreshold);

    // Penalize high energy
    score *= Math.max(0, 1 - state.energy / this.energyMax);

    // Chiral stability bonus/penalty
    if (state.chiral) {
      const chiralScore = this.calculateChiralStabilityScore(state.chiral);
      // Chiral stability can boost or reduce overall score by up to 20%
      score *= (0.8 + 0.4 * chiralScore);
    }

    return Math.max(0, Math.min(1, score));
  }

  /**
   * Calculate chiral contribution to stability score
   */
  private calculateChiralStabilityScore(chiral: ChiralState): number {
    let score = 1.0;

    // Reward being in stable regime
    if (chiral.stabilityRegime === 'stable') {
      score *= 1.0;
    } else if (chiral.stabilityRegime === 'transitional') {
      score *= 0.7;
    } else {
      score *= 0.3;
    }

    // Moderate asymmetry is good (0.3-0.7 range optimal)
    const optimalAsymmetry = 0.5;
    const asymmetryDeviation = Math.abs(chiral.asymmetry - optimalAsymmetry);
    score *= Math.max(0.5, 1 - asymmetryDeviation);

    // Topological protection bonus
    score *= TOPOLOGICAL_PROTECTION_FACTOR + (1 - TOPOLOGICAL_PROTECTION_FACTOR) *
      (1 - Math.min(1, Math.abs(chiral.windingNumber) / 3));

    return Math.max(0, Math.min(1, score));
  }

  /**
   * Get recommended damping adjustment
   */
  recommendDampingAdjustment(state: SystemState): number {
    const status = this.check(state);

    switch (status.action) {
      case 'reduce-coupling':
        return 0.95;  // Reduce by 5%
      case 'boost-coherence':
        return 1.05;  // Increase by 5%
      case 'normalize':
        return 0.9;   // Reduce by 10%
      case 'emergency-halt':
        return DAMPING_MAX / state.energy;  // Strong damping
      default:
        return 1.0;   // No change
    }
  }

  /**
   * Check if system is in conscious regime
   * Based on Kuramoto synchronization research
   */
  isConsciousRegime(coherence: number): boolean {
    // Consciousness exists in the regime between
    // full synchronization and complete decoherence
    return coherence >= 0.55 && coherence <= 0.85;
  }

  /**
   * Get safety status summary
   */
  getSummary(state: SystemState): string {
    const status = this.check(state);

    if (status.status === 'nominal') {
      return `✓ System nominal (stability: ${(status.metrics.stabilityScore * 100).toFixed(1)}%)`;
    }

    if (status.status === 'warning') {
      return `⚠ Warning: ${status.violations[0]} → ${status.action}`;
    }

    return `🛑 Critical: ${status.violations.length} violations → ${status.action}`;
  }

  /**
   * Reset safety monitoring
   */
  reset(): void {
    this.history = [];
  }
}

// ============================================
// CONVENIENCE EXPORTS
// ============================================

export const createSafetyEnvelope = (config?: SafetyConfig) =>
  new SafetyEnvelope(config);

/**
 * Quick safety check for common use cases
 */
export const quickCheck = (
  coherence: number,
  emergenceNorm: number,
  energy: number
): SafetyAction => {
  const envelope = new SafetyEnvelope();
  return envelope.check({ coherence, emergenceNorm, energy }).action;
};
