# Quantum Circuit Engine (QCE)

[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)]()
[![Build System](https://img.shields.io/badge/build-CMake-green.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)]()

A modular quantum circuit simulation framework implemented in modern C++17, providing explicit statevector-based modeling of quantum computation.

---

## Overview

Quantum Circuit Engine (QCE) is a from-first-principles implementation of quantum circuit simulation using explicit complex statevector evolution. The framework models qubit registers, unitary transformations, and probabilistic measurement through linear algebra operations over complex amplitudes.

The architecture cleanly separates core simulation primitives from higher-level algorithm modules, enabling structured experimentation with quantum algorithms, error correction schemes, and variational optimization methods. The system is designed with modular compilation units, backend extensibility, and performance-aware memory management.

---

## Core Capabilities

### 1. Statevector Simulation

- Explicit 2^n-dimensional complex amplitude representation  
- Deterministic unitary evolution via matrix–vector transformations  
- Probabilistic measurement with state collapse  
- Bitwise qubit indexing and amplitude manipulation  

### 2. Quantum Gates

- Single-qubit gates: X, Y, Z, H  
- Multi-qubit gates: CNOT  
- Parameterized rotation gates  
- Support for extensible custom gate definitions  

### 3. Implemented Algorithms

- Grover’s Search  
- Quantum Fourier Transform (QFT)  
- Quantum Teleportation  

### 4. Error Correction

- Three-qubit bit-flip code  
- Three-qubit phase-flip code  
- Modular framework for extending to larger codes  

### 5. Variational & Hybrid Methods

- Variational Quantum Eigensolver (VQE)  
- Quantum Approximate Optimization Algorithm (QAOA)  
- Parameterized quantum circuits for hybrid classical–quantum experimentation  

---

## Project Structure

```text
qce/
│
├── core/              # Core primitives (QubitRegister, QuantumGate, utilities)
├── simulators/        # Backend implementations (statevector, experimental)
├── algorithms/        # Algorithm implementations (Grover, QFT, teleportation)
├── error_correction/  # Error correction modules
├── optimization/      # Variational methods (VQE, QAOA)
├── examples/          # Demonstration circuits
├── tests/             # Unit and integration tests
├── cmake/             # CMake configuration modules
└── CMakeLists.txt
---

## Design Principles

- Explicit linear algebra–based state evolution  
- Clear separation of simulation backend and algorithm layer  
- Backend extensibility by design  
- Performance-conscious memory handling in C++17  
- Modular architecture with CMake build system  

---

## Mathematical Foundation

QCE represents an n-qubit system as a complex vector in a 2^n-dimensional Hilbert space. Gate application corresponds to multiplication by unitary matrices, while measurement is modeled via probabilistic amplitude collapse. Variational modules implement parameterized unitary layers optimized through classical feedback loops.

---

## Build Instructions

```bash
mkdir build
cd build
cmake ..
make
```
