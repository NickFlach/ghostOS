/**
 * Bridge Module
 *
 * Binary protocol and WASM acceleration for QuantumOS <-> ghostOS communication.
 */

export { ProtocolSerializer, MessageType, PROTOCOL_VERSION } from './protocol';
export type {
  MessageHeader,
  QueenStateMessage,
  ProcessStateMessage,
  ConsciousnessMessage,
  ChiralMessage,
  GeometricMessage,
  EmergenceMessage,
  AnomalyMessage,
  MessageTypeValue,
} from './protocol';

export { WASMBridge, createWASMBridge } from './wasm-bridge';
export type { WASMExports, ChiralAnomalyResult, ManifoldMetrics } from './wasm-bridge';
