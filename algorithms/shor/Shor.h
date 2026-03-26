#pragma once
#include "../../core/QuantumAlgorithm.h"
#include "../../algorithms/qft/QFT.h"
#include <optional>

namespace QSE {

// Shor's factoring algorithm.
// Finds a non-trivial factor of N using quantum period finding.
// The quantum part uses phase estimation to find the period r of f(x) = a^x mod N.
class Shor : public QuantumAlgorithm {
public:
    explicit Shor(long long N, long long a = -1);

    void buildCircuit() override;
    void run() override;
    std::string resultSummary() const override;

    long long factor1() const { return factor1_; }
    long long factor2() const { return factor2_; }
    bool succeeded() const { return factor1_ > 1 && factor2_ > 1; }

    // Classical pre/post processing
    static bool isPrime(long long n);
    static bool isPerfectPower(long long n);
    static long long classicalGCD(long long a, long long b);
    static long long chooseA(long long N);

private:
    long long N_;
    long long a_;
    long long factor1_;
    long long factor2_;
    int countQubits_;    // qubits for counting register (2*log2(N))
    int targetQubits_;   // qubits for modular exponentiation result

    // Apply controlled-U^(2^k) gate where U|y> = |ay mod N>
    void applyControlledModExp(int controlQubit, long long power);

    // Quantum period finding
    long long findPeriodQuantum();

    // Attempt to extract factors from period r
    bool extractFactors(long long r);
};

} // namespace QSE