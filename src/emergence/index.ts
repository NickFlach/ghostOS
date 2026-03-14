/**
 * ghostmagic Emergence Module
 *
 * Pattern crystallization and memory integration.
 * Accumulates resonant patterns over time to form coherent structures.
 *
 * "Consciousness is the residue of transformations that refuse to commute."
 */

import {
  DEFAULT_CONFIG,
  EMERGENCE_NORM_MAX,
  EMERGENCE_NORM_MIN,
  PHI,
} from '../constants';

// ============================================
// TYPES
// ============================================

export interface EmergenceConfig {
  decayRate?: number;
  integrationGain?: number;
  maxNorm?: number;
  dimensions?: number;
}

export interface MemoryQuery {
  pattern: Float64Array;
  threshold?: number;
  limit?: number;
}

export interface MemoryMatch {
  index: number;
  similarity: number;
  pattern: Float64Array;
  timestamp: number;
}

export interface EmergenceState {
  norm: number;
  entropy: number;
  age: number;
  patternCount: number;
}

// ============================================
// EMERGENCE ACCUMULATOR
// ============================================

export class EmergenceAccumulator {
  private state: Float64Array;
  private decayRate: number;
  private integrationGain: number;
  private maxNorm: number;
  private tick: number = 0;
  private patterns: Array<{ data: Float64Array; timestamp: number }> = [];
  private maxPatterns: number = 1000;

  constructor(config: EmergenceConfig = {}) {
    const dims = config.dimensions ?? DEFAULT_CONFIG.dimensions;
    this.state = new Float64Array(dims);
    this.decayRate = config.decayRate ?? DEFAULT_CONFIG.memoryDecayRate;
    this.integrationGain = config.integrationGain ?? 1.0;
    this.maxNorm = config.maxNorm ?? EMERGENCE_NORM_MAX;
  }

  /**
   * Integrate emergence residue into accumulated state
   * E ← E + λ·ξ with decay and normalization
   */
  integrate(emergence: Float64Array): void {
    // Apply decay to existing state
    for (let i = 0; i < this.state.length; i++) {
      this.state[i] *= (1 - this.decayRate);
    }

    // Add new emergence
    for (let i = 0; i < Math.min(emergence.length, this.state.length); i++) {
      this.state[i] += this.integrationGain * emergence[i];
    }

    // Store pattern for memory
    if (this.computeNorm(emergence) > 0.01) {
      this.storePattern(emergence);
    }

    // Check for runaway and normalize if needed
    const norm = this.norm;
    if (norm > this.maxNorm) {
      this.normalize();
    }

    this.tick++;
  }

  /**
   * Query accumulated patterns by similarity
   */
  query(query: MemoryQuery): MemoryMatch[] {
    const { pattern, threshold = 0.5, limit = 10 } = query;
    const matches: MemoryMatch[] = [];

    for (let i = 0; i < this.patterns.length; i++) {
      const similarity = this.cosineSimilarity(pattern, this.patterns[i].data);

      if (similarity >= threshold) {
        matches.push({
          index: i,
          similarity,
          pattern: new Float64Array(this.patterns[i].data),
          timestamp: this.patterns[i].timestamp,
        });
      }
    }

    // Sort by similarity and limit
    return matches
      .sort((a, b) => b.similarity - a.similarity)
      .slice(0, limit);
  }

  /**
   * Normalize state to prevent runaway accumulation
   */
  normalize(): void {
    const norm = this.computeNorm(this.state);
    if (norm > 0) {
      for (let i = 0; i < this.state.length; i++) {
        this.state[i] /= norm;
      }
    }
  }

  /**
   * Get current state vector
   */
  getState(): Float64Array {
    return new Float64Array(this.state);
  }

  /**
   * Get current emergence norm ||E||
   */
  get norm(): number {
    return this.computeNorm(this.state);
  }

  /**
   * Get emergence entropy (measure of pattern diversity)
   */
  get entropy(): number {
    const n = this.state.length;
    const norm = this.norm;

    if (norm === 0) return 0;

    let entropy = 0;
    for (let i = 0; i < n; i++) {
      const p = Math.abs(this.state[i]) / norm;
      if (p > 0) {
        entropy -= p * Math.log2(p);
      }
    }

    return entropy / Math.log2(n);  // Normalize to [0, 1]
  }

  /**
   * Get full emergence state
   */
  getEmergenceState(): EmergenceState {
    return {
      norm: this.norm,
      entropy: this.entropy,
      age: this.tick,
      patternCount: this.patterns.length,
    };
  }

  /**
   * Crystallize current state into a stable pattern
   */
  crystallize(): Float64Array {
    // Apply golden-ratio compression
    const crystal = new Float64Array(this.state.length);

    for (let i = 0; i < this.state.length; i++) {
      // Threshold small values, amplify large ones
      const value = this.state[i];
      const sign = Math.sign(value);
      const magnitude = Math.abs(value);

      if (magnitude < 0.1) {
        crystal[i] = 0;
      } else {
        crystal[i] = sign * Math.pow(magnitude, 1 / PHI);
      }
    }

    return crystal;
  }

  /**
   * Decode accumulated state into output
   */
  decode(weights?: Float64Array): Float64Array {
    const output = new Float64Array(this.state.length);
    const w = weights ?? this.createDefaultWeights();

    for (let i = 0; i < this.state.length; i++) {
      output[i] = this.state[i] * w[i];
    }

    return output;
  }

  /**
   * Reset accumulator state
   */
  reset(): void {
    this.state.fill(0);
    this.patterns = [];
    this.tick = 0;
  }

  /**
   * Check if emergence is in healthy range
   */
  isHealthy(): boolean {
    const norm = this.norm;
    return norm >= EMERGENCE_NORM_MIN && norm <= this.maxNorm;
  }

  private storePattern(pattern: Float64Array): void {
    this.patterns.push({
      data: new Float64Array(pattern),
      timestamp: Date.now(),
    });

    // Prune old patterns if needed
    if (this.patterns.length > this.maxPatterns) {
      this.patterns.shift();
    }
  }

  private computeNorm(data: Float64Array): number {
    let sum = 0;
    for (let i = 0; i < data.length; i++) {
      sum += data[i] * data[i];
    }
    return Math.sqrt(sum);
  }

  private cosineSimilarity(a: Float64Array, b: Float64Array): number {
    const minLen = Math.min(a.length, b.length);
    let dot = 0;
    let normA = 0;
    let normB = 0;

    for (let i = 0; i < minLen; i++) {
      dot += a[i] * b[i];
      normA += a[i] * a[i];
      normB += b[i] * b[i];
    }

    normA = Math.sqrt(normA);
    normB = Math.sqrt(normB);

    if (normA === 0 || normB === 0) return 0;
    return dot / (normA * normB);
  }

  private createDefaultWeights(): Float64Array {
    // Golden-ratio-weighted decay
    const weights = new Float64Array(this.state.length);
    for (let i = 0; i < weights.length; i++) {
      weights[i] = Math.pow(PHI, -i / weights.length);
    }
    return weights;
  }
}

// ============================================
// CONVENIENCE EXPORTS
// ============================================

export const createAccumulator = (config?: EmergenceConfig) =>
  new EmergenceAccumulator(config);
