/**
 * ghostmagic Signal Processing Module
 *
 * Handles input processing, embedding, filtering, and normalization.
 * Transforms raw signals into the vector space where resonance occurs.
 */

import { DEFAULT_CONFIG, PHI, TAU } from '../constants';

// ============================================
// TYPES
// ============================================

export interface SignalOptions {
  dimensions?: number;
  sampleRate?: number;
}

export interface FilterOptions {
  lowCutoff?: number;
  highCutoff?: number;
  type: 'lowpass' | 'highpass' | 'bandpass';
}

export interface Signal {
  data: Float64Array;
  dimensions: number;
  timestamp: number;
  metadata?: Record<string, unknown>;
}

// ============================================
// SIGNAL CLASS
// ============================================

export class SignalProcessor {
  readonly dimensions: number;
  readonly sampleRate: number;

  constructor(options: SignalOptions = {}) {
    this.dimensions = options.dimensions ?? DEFAULT_CONFIG.dimensions;
    this.sampleRate = options.sampleRate ?? DEFAULT_CONFIG.sampleRate;
  }

  /**
   * Create a Signal from raw input data
   */
  from(input: number[] | Float64Array): Signal {
    const data = input instanceof Float64Array
      ? input
      : new Float64Array(input);

    return {
      data,
      dimensions: data.length,
      timestamp: Date.now(),
    };
  }

  /**
   * Embed raw data into a higher-dimensional vector space
   * Uses random Fourier features for approximate kernel embedding
   */
  embed(signal: Signal, targetDimensions?: number): Signal {
    const dims = targetDimensions ?? this.dimensions;
    const embedded = new Float64Array(dims);

    // Random Fourier feature embedding
    for (let i = 0; i < dims; i++) {
      let value = 0;
      const freq = (i + 1) * PHI;  // Golden-ratio-spaced frequencies

      for (let j = 0; j < signal.data.length; j++) {
        value += Math.cos(freq * signal.data[j] + (i * Math.PI) / dims);
      }

      embedded[i] = value / signal.data.length;
    }

    return {
      data: embedded,
      dimensions: dims,
      timestamp: signal.timestamp,
      metadata: { ...signal.metadata, embedded: true },
    };
  }

  /**
   * Apply frequency-domain filtering
   */
  filter(signal: Signal, options: FilterOptions): Signal {
    const { lowCutoff = 0, highCutoff = Infinity, type } = options;
    const filtered = new Float64Array(signal.data.length);

    // Simple IIR filter approximation
    const alpha = this.calculateAlpha(lowCutoff, highCutoff, type);

    let prev = 0;
    for (let i = 0; i < signal.data.length; i++) {
      const current = signal.data[i];

      switch (type) {
        case 'lowpass':
          filtered[i] = prev + alpha * (current - prev);
          break;
        case 'highpass':
          filtered[i] = alpha * (prev + current - (i > 0 ? signal.data[i - 1] : 0));
          break;
        case 'bandpass':
          // Cascade lowpass then highpass
          const lp = prev + alpha * (current - prev);
          filtered[i] = alpha * (lp - (i > 0 ? filtered[i - 1] : 0));
          break;
      }

      prev = filtered[i];
    }

    return {
      data: filtered,
      dimensions: signal.dimensions,
      timestamp: signal.timestamp,
      metadata: { ...signal.metadata, filtered: true, filterType: type },
    };
  }

  /**
   * Normalize signal to unit sphere
   */
  normalize(signal: Signal): Signal {
    const norm = this.computeNorm(signal.data);

    if (norm === 0) {
      return signal;  // Avoid division by zero
    }

    const normalized = new Float64Array(signal.data.length);
    for (let i = 0; i < signal.data.length; i++) {
      normalized[i] = signal.data[i] / norm;
    }

    return {
      data: normalized,
      dimensions: signal.dimensions,
      timestamp: signal.timestamp,
      metadata: { ...signal.metadata, normalized: true, originalNorm: norm },
    };
  }

  /**
   * Compute L2 norm of a vector
   */
  computeNorm(data: Float64Array): number {
    let sum = 0;
    for (let i = 0; i < data.length; i++) {
      sum += data[i] * data[i];
    }
    return Math.sqrt(sum);
  }

  /**
   * Compute dot product between two signals
   */
  dot(a: Signal, b: Signal): number {
    if (a.dimensions !== b.dimensions) {
      throw new Error('Signal dimensions must match');
    }

    let sum = 0;
    for (let i = 0; i < a.data.length; i++) {
      sum += a.data[i] * b.data[i];
    }
    return sum;
  }

  /**
   * Compute cosine similarity between two signals
   */
  cosineSimilarity(a: Signal, b: Signal): number {
    const dotProduct = this.dot(a, b);
    const normA = this.computeNorm(a.data);
    const normB = this.computeNorm(b.data);

    if (normA === 0 || normB === 0) return 0;
    return dotProduct / (normA * normB);
  }

  private calculateAlpha(low: number, high: number, type: string): number {
    const dt = 1 / this.sampleRate;
    const rc = type === 'lowpass'
      ? 1 / (TAU * high)
      : type === 'highpass'
        ? 1 / (TAU * low)
        : 1 / (TAU * Math.sqrt(low * high));

    return dt / (rc + dt);
  }
}

// ============================================
// CONVENIENCE EXPORTS
// ============================================

const defaultProcessor = new SignalProcessor();

export const embed = (signal: Signal, options?: { dimensions?: number }) =>
  defaultProcessor.embed(signal, options?.dimensions);

export const filter = (signal: Signal, options: FilterOptions) =>
  defaultProcessor.filter(signal, options);

export const normalize = (signal: Signal) =>
  defaultProcessor.normalize(signal);

export const Signal = {
  from: (data: number[] | Float64Array) => defaultProcessor.from(data),
};
