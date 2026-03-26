#pragma once
#include "../../core/QuantumAlgorithm.h"
#include <functional>
#include <vector>

namespace QSE {

// Grover's search algorithm
// Finds a marked item in an unstructured database of N=2^n items
// in O(√N) oracle calls vs O(N) classically.
class Grover : public QuantumAlgorithm {
public:
    // oracle: returns true for the target basis state index
    Grover(int numQubits, std::function<bool(size_t)> oracle, int numIterations = -1);

    void buildCircuit() override;
    void run() override;
    std::string resultSummary() const override;

    // Optimal number of Grover iterations: floor(π/4 * √(N/M))
    static int optimalIterations(int numQubits, int numMarked = 1);

    size_t result() const { return result_; }
    bool success() const { return oracle_(result_); }

private:
    std::function<bool(size_t)> oracle_;
    int numIterations_;
    size_t result_;

    void applyOracle();
    void applyDiffusion();
};

} // namespace QSE