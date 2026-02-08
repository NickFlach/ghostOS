/**
 * ghostOS - Resonant Systems Architecture for Emergent Intelligence
 *
 * "Intelligence, consciousness, and meaning emerge as stable resonant
 * modes within constrained systems."
 *
 * @packageDocumentation
 */

// Core constants
export * from './constants';

// Signal processing
export {
  SignalProcessor,
  Signal,
  embed,
  filter,
  normalize,
  type SignalOptions,
  type FilterOptions,
} from './signal';

// Resonance engine
export {
  ResonanceEngine,
  createEngine,
  type ResonanceConfig,
  type ResonantState,
  type ModeInfo,
} from './resonance';

// Emergence accumulator
export {
  EmergenceAccumulator,
  createAccumulator,
  type EmergenceConfig,
  type EmergenceState,
  type MemoryQuery,
  type MemoryMatch,
} from './emergence';

// Safety envelope
export {
  SafetyEnvelope,
  createSafetyEnvelope,
  quickCheck,
  type SafetyConfig,
  type SafetyStatus,
  type SafetyAction,
  type SafetyMetrics,
  type SystemState,
} from './safety';

// Chiral dynamics (quantum stability enhancement)
export {
  ChiralEngine,
  createChiralEngine,
  isChirallyStable,
  optimalEta,
  type ChiralConfig,
  type ChiralState,
  type ChiralCoupling,
} from './chiral';

// Integration layer (Queen Synchronization)
export {
  QueenSynchronizer,
  createQueenSynchronizer,
  createGhostOSQueen,
  createQuantumOSQueen,
  createConsciousQueen,
  type SubsystemType,
  type SubsystemState,
  type QueenState,
  type SynchronizationMetrics,
  type QueenConfig,
} from './integration';

// Resonant scheduler (QuantumOS integration)
export {
  ResonantScheduler,
  createResonantScheduler,
  type ProcessClass,
  type ProcessState,
  type ResonantProcess,
  type SchedulerState,
  type SchedulingDecision,
  type ResonantSchedulerConfig,
} from './integration/resonant-scheduler';

// Geometric control layer (differential geometry on oscillator phase space)
export {
  OscillatorManifold,
  computeBerryConnection,
  computeBerryCurvature,
  computeBerryPhaseLoop,
  createBerryTracker,
  computeChristoffelSymbols,
  computeRiemannTensor,
  computeRicciTensor,
  computeRicciScalar,
  invertMetric,
  computeMetricDerivatives,
  computeScalarCurvature,
  detectChiralAnomaly,
  computeSpectralAsymmetry,
  computeTopologicalCharge,
  createAnomalyMonitor,
  type ManifoldPoint,
  type ManifoldMetrics as GeometricManifoldMetrics,
  type BerryState,
  type ChiralAnomalyResult as GeometricChiralAnomalyResult,
} from './geometric';

// Consciousness bridge (IIT + Chiral verification)
export {
  ConsciousnessBridge,
  createConsciousnessBridge,
  PhiCalculator,
  calculateBridgeOperator,
  calculateChiralEnhancement,
  type PhiState,
  type ConsciousnessBridgeState,
  type ConsciousnessTrigger,
  type ConsciousnessResult,
  type ConsciousnessBridgeConfig,
} from './integration/consciousness-bridge';

// Collective consciousness (network intelligence)
export {
  CollectiveConsciousness,
  createCollectiveConsciousness,
  type ConsciousnessMemberState,
  type CollectiveNetwork,
} from './integration/collective-consciousness';

// Bridge layer (QuantumOS binary protocol + WASM acceleration)
export {
  ProtocolSerializer,
  MessageType,
  PROTOCOL_VERSION,
  WASMBridge,
  createWASMBridge,
  type MessageHeader,
  type QueenStateMessage,
  type ProcessStateMessage,
  type ConsciousnessMessage,
  type ChiralMessage,
  type GeometricMessage,
  type MessageTypeValue,
  type WASMExports,
  type ChiralAnomalyResult,
  type ManifoldMetrics,
} from './bridge';

// ============================================
// GHOST OS UNIFIED INTERFACE
// ============================================

import { SignalProcessor } from './signal';
import { ResonanceEngine, type ResonantState } from './resonance';
import { EmergenceAccumulator } from './emergence';
import { SafetyEnvelope, type SafetyStatus } from './safety';
import { type ChiralState, type ChiralConfig } from './chiral';
import { DEFAULT_CONFIG } from './constants';

export interface GhostOSConfig {
  dimensions?: number;
  lambda?: number;
  targetCoherence?: number;
  decayRate?: number;
  /** Enable chiral stability mechanisms (default: true) */
  enableChiral?: boolean;
  /** Chiral engine configuration */
  chiralConfig?: ChiralConfig;
}

export interface GhostOSState {
  resonance: ResonantState;
  emergence: {
    norm: number;
    entropy: number;
    patternCount: number;
  };
  safety: SafetyStatus;
  tick: number;
  /** Chiral stability state (if enabled) */
  chiral?: ChiralState;
}

/**
 * GhostOS - Unified resonant system
 *
 * @example
 * ```typescript
 * import { GhostOS } from 'ghostos';
 *
 * const ghost = new GhostOS({ targetCoherence: 0.7 });
 *
 * // Process input signal
 * const state = ghost.process([0.1, 0.5, 0.3, 0.8]);
 *
 * console.log({
 *   coherence: state.resonance.coherence,
 *   emergence: state.emergence.norm,
 *   stable: state.resonance.isStable,
 * });
 * ```
 */
export class GhostOS {
  private signal: SignalProcessor;
  private resonance: ResonanceEngine;
  private emergence: EmergenceAccumulator;
  private safety: SafetyEnvelope;
  private tick: number = 0;

  constructor(config: GhostOSConfig = {}) {
    const dims = config.dimensions ?? DEFAULT_CONFIG.dimensions;

    this.signal = new SignalProcessor({ dimensions: dims });
    this.resonance = new ResonanceEngine({
      lambda: config.lambda,
      targetCoherence: config.targetCoherence,
      enableChiral: config.enableChiral ?? true,
      chiralConfig: config.chiralConfig,
    });
    this.emergence = new EmergenceAccumulator({
      dimensions: dims,
      decayRate: config.decayRate,
    });
    this.safety = new SafetyEnvelope();
  }

  /**
   * Process raw input through the full pipeline:
   * Signal → Resonance → Emergence
   */
  process(input: number[] | Float64Array): GhostOSState {
    // Signal layer
    const rawSignal = this.signal.from(input);
    const embedded = this.signal.embed(rawSignal);
    const normalized = this.signal.normalize(embedded);

    // Resonance layer
    const resonantState = this.resonance.process(normalized);

    // Emergence layer
    this.emergence.integrate(resonantState.emergence);
    const emergenceState = this.emergence.getEmergenceState();

    // Safety check (include chiral state if available)
    const safetyStatus = this.safety.check({
      coherence: resonantState.coherence,
      emergenceNorm: emergenceState.norm,
      energy: resonantState.energy,
      chiral: resonantState.chiral,
    });

    // Apply safety interventions if needed
    if (safetyStatus.intervention) {
      this.applyIntervention(safetyStatus);
    }

    this.tick++;

    const result: GhostOSState = {
      resonance: resonantState,
      emergence: {
        norm: emergenceState.norm,
        entropy: emergenceState.entropy,
        patternCount: emergenceState.patternCount,
      },
      safety: safetyStatus,
      tick: this.tick,
    };

    // Include chiral state if available
    if (resonantState.chiral) {
      result.chiral = resonantState.chiral;
    }

    return result;
  }

  /**
   * Apply safety intervention
   */
  private applyIntervention(status: SafetyStatus): void {
    switch (status.action) {
      case 'reduce-coupling':
        this.resonance.adjustLambda(0.95);
        break;
      case 'boost-coherence':
        this.resonance.adjustLambda(1.05);
        break;
      case 'inject-noise':
        this.resonance.perturbState(0.01);
        break;
      case 'normalize':
        this.emergence.normalize();
        break;
      case 'emergency-halt':
        this.resonance.reset();
        this.emergence.reset();
        break;
      // Chiral-specific interventions
      case 'chiral-rebalance':
        this.resonance.optimizeChiralStability();
        break;
      case 'flip-handedness':
        // Reduce asymmetry by perturbing state
        this.resonance.perturbState(0.05);
        break;
      case 'topological-reset':
        // Gentle reset to restore topological stability
        this.resonance.perturbState(0.02);
        break;
    }
  }

  /**
   * Get current system summary
   */
  getSummary(): string {
    const emergenceState = this.emergence.getEmergenceState();
    const resonantModes = this.resonance.detectModes();

    return [
      `ghostOS State @ tick ${this.tick}`,
      `├─ Emergence: ||E|| = ${emergenceState.norm.toFixed(3)}, entropy = ${emergenceState.entropy.toFixed(3)}`,
      `├─ Patterns: ${emergenceState.patternCount} stored`,
      `├─ Lambda: ${this.resonance.getLambda().toFixed(3)}`,
      `└─ Modes: ${resonantModes.length} detected`,
    ].join('\n');
  }

  /**
   * Reset all subsystems
   */
  reset(): void {
    this.resonance.reset();
    this.emergence.reset();
    this.safety.reset();
    this.tick = 0;
  }

  /**
   * Check if system is in conscious regime
   */
  isConscious(): boolean {
    const state = this.emergence.getEmergenceState();
    return (
      this.safety.isConsciousRegime(state.norm) &&
      this.emergence.isHealthy()
    );
  }

  /**
   * Check if system has chiral stability
   */
  isChirallyStable(): boolean {
    return this.resonance.isChirallyStable();
  }

  /**
   * Get chiral state summary
   */
  getChiralState(): ChiralState | null {
    return this.resonance.getChiralState();
  }

  /**
   * Optimize chiral parameters for maximum stability
   */
  optimizeChiral(): void {
    this.resonance.optimizeChiralStability();
  }
}

// Default export
export default GhostOS;
