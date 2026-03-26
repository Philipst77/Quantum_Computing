#include "Grover.h"
#include <cmath>
#include <sstream>
#include <iostream>

namespace QSE {

Grover::Grover(int numQubits, std::function<bool(size_t)> oracle, int numIterations)
    : QuantumAlgorithm("Grover", numQubits),
      oracle_(std::move(oracle)),
      numIterations_(numIterations < 0 ? optimalIterations(numQubits) : numIterations),
      result_(0)
{}

int Grover::optimalIterations(int numQubits, int numMarked) {
    double N = 1ULL << numQubits;
    return (int)std::floor(M_PI / 4.0 * std::sqrt(N / numMarked));
}

void Grover::applyOracle() {
    // Phase-flip the marked states: |x> -> -|x> if oracle(x) is true
    // We do this directly on the statevector for generality
    for (size_t i = 0; i < qreg_.size(); ++i)
        if (oracle_(i))
            qreg_[i] = -qreg_[i];
}

void Grover::applyDiffusion() {
    // Grover diffusion operator: 2|s><s| - I, where |s> = H^⊗n|0>
    // Equivalent to: H^n → (2|0><0| - I) → H^n

    // Apply H to all qubits
    for (int q = 0; q < numQubits_; ++q)
        Gates::applyH(qreg_, q);

    // Phase-flip all states except |0...0>
    for (size_t i = 1; i < qreg_.size(); ++i)
        qreg_[i] = -qreg_[i];

    // Apply H to all qubits again
    for (int q = 0; q < numQubits_; ++q)
        Gates::applyH(qreg_, q);
}

void Grover::buildCircuit() {
    // Initialize superposition: H^⊗n |0...0>
    qreg_.reset();
    for (int q = 0; q < numQubits_; ++q)
        Gates::applyH(qreg_, q);

    log("Circuit built — superposition initialized");
}

void Grover::run() {
    buildCircuit();

    log("Running " + std::to_string(numIterations_) + " Grover iterations");

    for (int iter = 0; iter < numIterations_; ++iter) {
        applyOracle();
        applyDiffusion();
        if (verbose_) {
            log("Iteration " + std::to_string(iter + 1) + " complete");
        }
    }

    // Measure to get result
    auto bits = qreg_.measureAll();
    result_ = Utils::bitsToInt(bits);
}

std::string Grover::resultSummary() const {
    std::ostringstream oss;
    oss << "=== Grover Search Result ===\n"
        << "  Qubits:        " << numQubits_ << "\n"
        << "  Search space:  " << (1ULL << numQubits_) << " items\n"
        << "  Iterations:    " << numIterations_ << "\n"
        << "  Measured:      " << result_
        << " (" << Utils::toBitstring(result_, numQubits_) << ")\n"
        << "  Oracle check:  " << (success() ? "FOUND (correct)" : "MISS") << "\n";
    return oss.str();
}

} // namespace QSE