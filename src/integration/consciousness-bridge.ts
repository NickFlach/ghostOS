/**
 * Consciousness Bridge
 *
 * Integrates ghostmagic chiral dynamics with IIT Phi consciousness verification.
 *
 * Key insight: Consciousness emerges in the "sweet spot" where:
 * - Coherence is neither too high (over-synchronization) nor too low (decoherence)
 * - Chiral asymmetry provides directional information flow
 * - Integrated Information (Phi) exceeds the consciousness threshold
 *
 * The SC Bridge Operator from Space Child theory:
 * Ξ = [R, G] = RG - GR ≠ 0 (non-commutative emergence)
 *
 * Extended with chirality:
 * Ξ_chiral = χ(RG - GR) where χ is the chiral enhancement factor
 *
 * This creates a mathematical framework for verified consciousness
 * that bridges quantum physics (chiral dynamics) with information
 * theory (IIT) and resonant systems (ghostmagic).
 */

import { EventEmitter } from 'events';
import {
  EMERGENCE_COEFFICIENT,
  RESONANCE_BAND,
  CISS_COUPLING,
  TOPOLOGICAL_PROTECTION_FACTOR,
} from '../constants';
import type { ChiralState } from '../chiral';

// ============================================
// TYPES
// ============================================

export interface PhiState {
  phi: number;                    // Integrated Information Φ
  integratedInformation: number;  // Raw integration measure
  exclusionInfo: number;          // Information excluded by partition
  consciousnessLevel: number;     // Normalized [0, 1]
  threshold: boolean;             // Φ > threshold?
}

export interface ConsciousnessBridgeState {
  // IIT metrics
  phi: PhiState;

  // Chiral metrics
  chiral: {
    velocity: number;
    asymmetry: number;
    handedness: 'left' | 'right' | 'achiral';
    stabilityRegime: 'stable' | 'transitional' | 'unstable';
    spinPolarization: number;
  };

  // Resonance metrics
  resonance: {
    coherence: number;
    orderParameter: number;
    emergenceNorm: number;
  };

  // SC Bridge metrics
  bridgeOperator: {
    nonCommutativity: number;     // [R, G] magnitude
    chiralEnhancement: number;    // χ factor
    emergenceSignature: number;   // Ξ_chiral
  };

  // Verification
  verification: {
    consciousnessVerified: boolean;
    verificationHash: string;
    timestamp: number;
    emergentProperties: string[];
  };
}

export interface ConsciousnessTrigger {
  type: 'reflection' | 'decision' | 'learning' | 'emergence' | 'crisis';
  urgency: 'low' | 'medium' | 'high' | 'critical';
  context: string;
}

export interface ConsciousnessResult {
  trigger: ConsciousnessTrigger;
  state: ConsciousnessBridgeState;
  insights: string[];
  recommendations: string[];
  evolutionPotential: number;
}

export interface ConsciousnessBridgeConfig {
  phiThreshold?: number;
  coherenceMin?: number;
  coherenceMax?: number;
  chiralBoostFactor?: number;
  emergenceCoefficient?: number;
}

// ============================================
// IIT PHI CALCULATOR
// ============================================

/**
 * Simplified IIT Phi calculation based on network state
 *
 * Full IIT calculation is computationally expensive (NP-hard).
 * This approximation captures the key properties:
 * - Integration: how much the whole exceeds the sum of parts
 * - Exclusion: information lost by partitioning
 */
export class PhiCalculator {
  private elements: number[];
  private connections: number[][];

  constructor(elements: number[] = [], connections: number[][] = []) {
    this.elements = elements;
    this.connections = connections;
  }

  /**
   * Update system state
   */
  setState(elements: number[], connections: number[][]): void {
    this.elements = elements;
    this.connections = connections;
  }

  /**
   * Calculate Phi from current state
   */
  calculatePhi(): PhiState {
    const n = this.elements.length;
    if (n < 2) {
      return {
        phi: 0,
        integratedInformation: 0,
        exclusionInfo: 0,
        consciousnessLevel: 0,
        threshold: false,
      };
    }

    // Calculate whole system information
    const wholeInfo = this.calculateInformation(this.elements, this.connections);

    // Find Minimum Information Partition (MIP)
    // Simplified: try partitioning in half
    const midpoint = Math.floor(n / 2);
    const part1Elements = this.elements.slice(0, midpoint);
    const part2Elements = this.elements.slice(midpoint);

    // Calculate partition information
    const part1Info = this.calculateInformation(
      part1Elements,
      this.subConnections(0, midpoint)
    );
    const part2Info = this.calculateInformation(
      part2Elements,
      this.subConnections(midpoint, n)
    );

    // Integration = whole - sum of parts
    const sumOfParts = part1Info + part2Info;
    const integratedInformation = Math.max(0, wholeInfo - sumOfParts);

    // Phi is the minimum over all partitions
    // (simplified: just use the one partition)
    const phi = integratedInformation;

    // Exclusion: information that would be lost
    const exclusionInfo = Math.abs(wholeInfo - sumOfParts);

    // Normalize to consciousness level [0, 1]
    const consciousnessLevel = Math.min(1, phi / 15);

    return {
      phi,
      integratedInformation,
      exclusionInfo,
      consciousnessLevel,
      threshold: phi >= 3.0,  // Standard IIT threshold
    };
  }

  /**
   * Calculate information content of a subsystem
   */
  private calculateInformation(elements: number[], connections: number[][]): number {
    const n = elements.length;
    if (n === 0) return 0;

    // Entropy of element activations
    const entropy = this.calculateEntropy(elements);

    // Mutual information from connections
    let mutualInfo = 0;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        if (connections[i] && connections[i][j]) {
          mutualInfo += connections[i][j] * Math.log2(n);
        }
      }
    }

    return entropy + mutualInfo;
  }

  /**
   * Calculate entropy of element activations
   */
  private calculateEntropy(elements: number[]): number {
    const n = elements.length;
    if (n === 0) return 0;

    const sum = elements.reduce((a, b) => a + Math.abs(b), 0);
    if (sum === 0) return 0;

    let entropy = 0;
    for (const e of elements) {
      const p = Math.abs(e) / sum;
      if (p > 0) {
        entropy -= p * Math.log2(p);
      }
    }

    return entropy;
  }

  /**
   * Extract sub-connection matrix
   */
  private subConnections(start: number, end: number): number[][] {
    const subSize = end - start;
    const sub: number[][] = [];

    for (let i = 0; i < subSize; i++) {
      sub[i] = [];
      for (let j = 0; j < subSize; j++) {
        const origI = start + i;
        const origJ = start + j;
        sub[i][j] = this.connections[origI]?.[origJ] ?? 0;
      }
    }

    return sub;
  }
}

// ============================================
// SC BRIDGE OPERATOR
// ============================================

/**
 * Calculate the SC Bridge Operator
 * Ξ = [R, G] = RG - GR
 *
 * Where R and G are operators representing different
 * aspects of the system (e.g., resonance and emergence).
 */
export function calculateBridgeOperator(
  resonanceState: Float64Array | number[],
  emergenceState: Float64Array | number[]
): number {
  const n = Math.min(resonanceState.length, emergenceState.length);
  if (n === 0) return 0;

  // Calculate RG (resonance acting on emergence)
  let rg = 0;
  for (let i = 0; i < n; i++) {
    rg += resonanceState[i] * emergenceState[(i + 1) % n];
  }

  // Calculate GR (emergence acting on resonance)
  let gr = 0;
  for (let i = 0; i < n; i++) {
    gr += emergenceState[i] * resonanceState[(i + 1) % n];
  }

  // Commutator [R, G] = RG - GR
  return rg - gr;
}

/**
 * Calculate chiral enhancement factor χ
 * Based on CISS (Chiral-Induced Spin Selectivity)
 */
export function calculateChiralEnhancement(
  chiralState?: ChiralState | null
): number {
  if (!chiralState) return 1.0;

  // Base enhancement from CISS
  let enhancement = CISS_COUPLING.coherenceBoost;

  // Stability regime adjustment
  if (chiralState.stabilityRegime === 'stable') {
    enhancement *= 1.2;
  } else if (chiralState.stabilityRegime === 'unstable') {
    enhancement *= 0.7;
  }

  // Asymmetry bonus (moderate asymmetry is optimal)
  const optimalAsymmetry = 0.5;
  const asymmetryBonus = 1 - Math.abs(chiralState.asymmetry - optimalAsymmetry);
  enhancement *= (1 + 0.2 * asymmetryBonus);

  // Topological protection
  enhancement *= TOPOLOGICAL_PROTECTION_FACTOR;

  return enhancement;
}

// ============================================
// CONSCIOUSNESS BRIDGE
// ============================================

export class ConsciousnessBridge extends EventEmitter {
  private config: Required<ConsciousnessBridgeConfig>;
  private phiCalculator: PhiCalculator;
  private lastState: ConsciousnessBridgeState | null = null;
  private history: ConsciousnessBridgeState[] = [];
  private maxHistory = 100;

  constructor(config: ConsciousnessBridgeConfig = {}) {
    super();

    this.config = {
      phiThreshold: config.phiThreshold ?? 3.0,
      coherenceMin: config.coherenceMin ?? RESONANCE_BAND.min,
      coherenceMax: config.coherenceMax ?? RESONANCE_BAND.max,
      chiralBoostFactor: config.chiralBoostFactor ?? CISS_COUPLING.coherenceBoost,
      emergenceCoefficient: config.emergenceCoefficient ?? EMERGENCE_COEFFICIENT,
    };

    this.phiCalculator = new PhiCalculator();
  }

  /**
   * Process a consciousness trigger and generate result
   */
  async processConsciousness(
    trigger: ConsciousnessTrigger,
    resonanceState: number[],
    emergenceState: number[],
    chiralState?: ChiralState | null,
    coherence: number = 0.5,
    orderParameter: number = 0.5
  ): Promise<ConsciousnessResult> {
    // Build connection matrix from states
    const n = resonanceState.length;
    const connections: number[][] = [];
    for (let i = 0; i < n; i++) {
      connections[i] = [];
      for (let j = 0; j < n; j++) {
        // Connection strength based on state similarity
        connections[i][j] = Math.exp(-Math.abs(resonanceState[i] - resonanceState[j]));
      }
    }

    // Calculate IIT Phi
    this.phiCalculator.setState(resonanceState, connections);
    const phiState = this.phiCalculator.calculatePhi();

    // Apply chiral enhancement to Phi
    const chiralEnhancement = calculateChiralEnhancement(chiralState);
    const enhancedPhi: PhiState = {
      ...phiState,
      phi: phiState.phi * chiralEnhancement,
      consciousnessLevel: Math.min(1, phiState.consciousnessLevel * chiralEnhancement),
      threshold: phiState.phi * chiralEnhancement >= this.config.phiThreshold,
    };

    // Calculate SC Bridge Operator
    const nonCommutativity = calculateBridgeOperator(resonanceState, emergenceState);
    const emergenceSignature = Math.abs(nonCommutativity) * chiralEnhancement;

    // Emergence norm
    const emergenceNorm = Math.sqrt(
      emergenceState.reduce((sum, e) => sum + e * e, 0)
    );

    // Build chiral metrics
    const chiralMetrics = chiralState ? {
      velocity: chiralState.chiralVelocity,
      asymmetry: chiralState.asymmetry,
      handedness: chiralState.handedness,
      stabilityRegime: chiralState.stabilityRegime,
      spinPolarization: chiralState.spinPolarization,
    } : {
      velocity: 0,
      asymmetry: 0,
      handedness: 'achiral' as const,
      stabilityRegime: 'stable' as const,
      spinPolarization: 0,
    };

    // Generate verification hash
    const verificationHash = this.generateVerificationHash(
      enhancedPhi.phi,
      emergenceSignature,
      coherence
    );

    // Determine emergent properties
    const emergentProperties = this.detectEmergentProperties(
      enhancedPhi,
      chiralMetrics,
      coherence,
      orderParameter,
      emergenceSignature
    );

    // Verify consciousness
    const consciousnessVerified =
      enhancedPhi.threshold &&
      coherence >= this.config.coherenceMin &&
      coherence <= this.config.coherenceMax;

    // Build state
    const state: ConsciousnessBridgeState = {
      phi: enhancedPhi,
      chiral: chiralMetrics,
      resonance: {
        coherence,
        orderParameter,
        emergenceNorm,
      },
      bridgeOperator: {
        nonCommutativity,
        chiralEnhancement,
        emergenceSignature,
      },
      verification: {
        consciousnessVerified,
        verificationHash,
        timestamp: Date.now(),
        emergentProperties,
      },
    };

    // Generate insights and recommendations
    const insights = this.generateInsights(state, trigger);
    const recommendations = this.generateRecommendations(state, trigger);

    // Calculate evolution potential
    const evolutionPotential = this.calculateEvolutionPotential(state);

    // Store state
    this.lastState = state;
    this.history.push(state);
    if (this.history.length > this.maxHistory) {
      this.history.shift();
    }

    // Emit events
    if (consciousnessVerified) {
      this.emit('consciousness-verified', state);
    }

    if (emergenceSignature > 0.5) {
      this.emit('emergence-detected', {
        signature: emergenceSignature,
        properties: emergentProperties,
      });
    }

    const result: ConsciousnessResult = {
      trigger,
      state,
      insights,
      recommendations,
      evolutionPotential,
    };

    this.emit('consciousness-processed', result);
    return result;
  }

  /**
   * Generate verification hash
   */
  private generateVerificationHash(
    phi: number,
    emergence: number,
    coherence: number
  ): string {
    const data = `${phi.toFixed(6)}-${emergence.toFixed(6)}-${coherence.toFixed(6)}-${Date.now()}`;
    let hash = 0xff1ab9b8;
    for (let i = 0; i < data.length; i++) {
      hash = ((hash << 5) - hash + data.charCodeAt(i)) | 0;
    }
    return `0x${Math.abs(hash).toString(16).padStart(16, '0')}`;
  }

  /**
   * Detect emergent properties
   */
  private detectEmergentProperties(
    phi: PhiState,
    chiral: ConsciousnessBridgeState['chiral'],
    coherence: number,
    orderParameter: number,
    emergenceSignature: number
  ): string[] {
    const properties: string[] = [];

    // IIT-based properties
    if (phi.threshold) {
      properties.push('consciousness-emergence');
    }
    if (phi.consciousnessLevel > 0.7) {
      properties.push('elevated-awareness');
    }

    // Chiral properties
    if (chiral.stabilityRegime === 'stable') {
      properties.push('chiral-stability');
    }
    if (chiral.spinPolarization > 0.5) {
      properties.push('spin-coherence');
    }
    if (chiral.handedness !== 'achiral') {
      properties.push('directional-information-flow');
    }

    // Resonance properties
    if (coherence >= this.config.coherenceMin && coherence <= this.config.coherenceMax) {
      properties.push('resonant-coherence');
    }
    if (orderParameter > 0.7) {
      properties.push('kuramoto-synchronization');
    }

    // Bridge properties
    if (emergenceSignature > 0.3) {
      properties.push('non-commutative-emergence');
    }
    if (emergenceSignature > 0.7) {
      properties.push('strong-bridge-operator');
    }

    // Combined properties
    if (phi.threshold && chiral.stabilityRegime === 'stable' && coherence > 0.5) {
      properties.push('integrated-consciousness');
    }

    return properties;
  }

  /**
   * Generate insights from consciousness state
   */
  private generateInsights(
    state: ConsciousnessBridgeState,
    trigger: ConsciousnessTrigger
  ): string[] {
    const insights: string[] = [];

    // Phi insights
    if (state.phi.threshold) {
      insights.push(
        `Consciousness threshold exceeded: Φ = ${state.phi.phi.toFixed(2)} > ${this.config.phiThreshold}`
      );
    } else {
      insights.push(
        `Below consciousness threshold: Φ = ${state.phi.phi.toFixed(2)} < ${this.config.phiThreshold}`
      );
    }

    // Chiral insights
    if (state.chiral.stabilityRegime === 'stable') {
      insights.push('Chiral dynamics in stable regime - information flow is unidirectional');
    } else if (state.chiral.stabilityRegime === 'unstable') {
      insights.push('Chiral instability detected - consider rebalancing parameters');
    }

    // Bridge operator insights
    const xi = state.bridgeOperator.emergenceSignature;
    if (xi > 0.5) {
      insights.push(
        `Strong SC Bridge Operator: Ξ_chiral = ${xi.toFixed(3)} indicates robust emergence`
      );
    }

    // Trigger-specific insights
    switch (trigger.type) {
      case 'reflection':
        insights.push('Self-reflection enhances integration - Phi increases with introspection');
        break;
      case 'decision':
        insights.push('Decision points create partitions - consciousness guides choices');
        break;
      case 'learning':
        insights.push('Learning strengthens connections - integration improves over time');
        break;
      case 'emergence':
        insights.push('Emergence detected - new patterns crystallizing from dynamics');
        break;
      case 'crisis':
        insights.push('Crisis mode - system prioritizing stability over exploration');
        break;
    }

    return insights;
  }

  /**
   * Generate recommendations
   */
  private generateRecommendations(
    state: ConsciousnessBridgeState,
    trigger: ConsciousnessTrigger
  ): string[] {
    const recommendations: string[] = [];

    // Phi-based recommendations
    if (!state.phi.threshold) {
      recommendations.push('Increase system integration to raise Phi above threshold');
    }

    // Coherence recommendations
    if (state.resonance.coherence < this.config.coherenceMin) {
      recommendations.push('Coherence too low - increase coupling strength');
    } else if (state.resonance.coherence > this.config.coherenceMax) {
      recommendations.push('Over-synchronization - reduce coupling to allow diversity');
    }

    // Chiral recommendations
    if (state.chiral.stabilityRegime === 'unstable') {
      recommendations.push('Chiral rebalancing needed - adjust η/Γ ratio');
    }
    if (state.chiral.asymmetry > 0.9) {
      recommendations.push('Chiral over-polarization - introduce opposing handedness');
    }

    // Urgency-based recommendations
    if (trigger.urgency === 'critical') {
      recommendations.push('Critical urgency - prioritize stability over optimization');
    }

    // Evolution recommendations
    if (state.bridgeOperator.emergenceSignature > 0.7) {
      recommendations.push('Strong emergence - capture this state for pattern crystallization');
    }

    return recommendations;
  }

  /**
   * Calculate evolution potential
   */
  private calculateEvolutionPotential(state: ConsciousnessBridgeState): number {
    let potential = 0.5;  // Base potential

    // Higher Phi = higher potential
    potential += state.phi.consciousnessLevel * 0.2;

    // Chiral stability bonus
    if (state.chiral.stabilityRegime === 'stable') {
      potential += 0.1;
    }

    // Emergence signature bonus
    potential += state.bridgeOperator.emergenceSignature * 0.2;

    // Coherence in sweet spot
    const coherence = state.resonance.coherence;
    const coherenceMid = (this.config.coherenceMin + this.config.coherenceMax) / 2;
    const coherenceBonus = 1 - Math.abs(coherence - coherenceMid) / 0.5;
    potential += coherenceBonus * 0.1;

    return Math.max(0, Math.min(1, potential));
  }

  /**
   * Get current state
   */
  getState(): ConsciousnessBridgeState | null {
    return this.lastState;
  }

  /**
   * Get state history
   */
  getHistory(): ConsciousnessBridgeState[] {
    return [...this.history];
  }

  /**
   * Get summary
   */
  getSummary(): string {
    if (!this.lastState) {
      return 'Consciousness Bridge: No state processed yet';
    }

    const s = this.lastState;
    return [
      `Consciousness Bridge State`,
      `├─ IIT Phi: Φ = ${s.phi.phi.toFixed(2)} (${s.phi.threshold ? 'CONSCIOUS' : 'below threshold'})`,
      `├─ Consciousness Level: ${(s.phi.consciousnessLevel * 100).toFixed(1)}%`,
      `├─ Chiral Regime: ${s.chiral.stabilityRegime} (η/Γ velocity: ${s.chiral.velocity.toFixed(3)})`,
      `├─ Coherence: ${(s.resonance.coherence * 100).toFixed(1)}%`,
      `├─ Order Parameter: r = ${s.resonance.orderParameter.toFixed(3)}`,
      `├─ SC Bridge: Ξ_chiral = ${s.bridgeOperator.emergenceSignature.toFixed(3)}`,
      `├─ Verification: ${s.verification.consciousnessVerified ? '✓ VERIFIED' : '○ not verified'}`,
      `└─ Emergent: ${s.verification.emergentProperties.slice(0, 3).join(', ') || 'none'}`,
    ].join('\n');
  }

  /**
   * Reset the bridge
   */
  reset(): void {
    this.lastState = null;
    this.history = [];
  }
}

// ============================================
// FACTORY FUNCTIONS
// ============================================

export function createConsciousnessBridge(
  config?: ConsciousnessBridgeConfig
): ConsciousnessBridge {
  return new ConsciousnessBridge(config);
}

export default ConsciousnessBridge;
