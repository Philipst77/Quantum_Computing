# QCSim-Lite

A modular **quantum computing simulator** written in C++17.  
This project is an independent reimplementation inspired by existing quantum simulators, designed for learning and extensibility.  

---

## 🚀 Features
- **Core simulator**
  - Qubit register and statevector representation
  - Support for common quantum gates (X, Y, Z, H, CNOT, rotations)
  - Probabilistic measurement
- **Algorithms**
  - Grover’s search
  - Quantum Fourier Transform (QFT)
  - Quantum teleportation
  - (More coming soon: Shor, BB84, Deutsch–Jozsa, VQE, QAOA)
- **Error Correction**
  - 3-qubit bit-flip and phase-flip codes
  - Shor’s 9-qubit code (planned)
- **Optimization**
  - Variational Quantum Eigensolver (VQE)
  - Quantum Approximate Optimization Algorithm (QAOA)

---

## 📂 Project Structure
```text
qcsim/
│
├── core/              # Core classes (QubitRegister, QuantumGate, Utils)
├── simulators/        # Different backends (statevector, MPS, Clifford)
├── algorithms/        # Quantum algorithms (Grover, QFT, Shor, etc.)
├── error_correction/  # Error correction codes
├── optimization/      # VQE, QAOA
├── examples/          # Example circuits (Bell state, teleportation, Grover)
├── tests/             # Unit & integration tests
├── cmake/             # CMake modules (optional)
└── CMakeLists.txt
