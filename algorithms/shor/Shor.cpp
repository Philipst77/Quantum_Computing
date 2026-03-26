#include "Shor.h"
#include <cmath>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <random>

namespace QSE {

Shor::Shor(long long N, long long a)
    : QuantumAlgorithm("Shor", (int)(2 * std::ceil(std::log2(N + 1)) + (int)std::ceil(std::log2(N + 1)))),
      N_(N), a_(a), factor1_(1), factor2_(1)
{
    if (N < 4) throw std::invalid_argument("N must be >= 4");
    int logN = (int)std::ceil(std::log2(N + 1));
    countQubits_  = 2 * logN;
    targetQubits_ = logN;

    if (a_ < 0) a_ = chooseA(N_);
}

bool Shor::isPrime(long long n) {
    if (n < 2) return false;
    if (n == 2) return true;
    if (n % 2 == 0) return false;
    for (long long i = 3; i * i <= n; i += 2)
        if (n % i == 0) return false;
    return true;
}

bool Shor::isPerfectPower(long long n) {
    for (int b = 2; b <= (int)std::log2(n) + 1; ++b) {
        long long a = (long long)std::round(std::pow(n, 1.0 / b));
        for (long long aa = std::max(2LL, a - 1); aa <= a + 1; ++aa)
            if ((long long)std::pow(aa, b) == n) return true;
    }
    return false;
}

long long Shor::classicalGCD(long long a, long long b) {
    return Utils::gcd(a, b);
}

long long Shor::chooseA(long long N) {
    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<long long> dist(2, N - 1);
    long long a;
    do { a = dist(rng); } while (classicalGCD(a, N) != 1);
    return a;
}

void Shor::applyControlledModExp(int controlQubit, long long power) {
    // In a full implementation this would apply the quantum modular
    // exponentiation circuit. For simulation purposes we implement it
    // as a classically-controlled phase based on the period structure.
    // A production implementation would decompose a^(2^k) mod N into
    // quantum arithmetic gates (Draper adders, etc.).
    (void)controlQubit; (void)power; // placeholder
}

long long Shor::findPeriodQuantum() {
    // Simulate quantum period finding using phase estimation.
    // The counting register starts in superposition via QFT.
    // We simulate the period finding classically as the full
    // quantum modular exponentiation circuit is resource-intensive.

    long long r = 1;
    long long val = a_ % N_;
    while (val != 1 && r < N_) {
        val = val * a_ % N_;
        ++r;
    }
    return (val == 1) ? r : 0;
}

bool Shor::extractFactors(long long r) {
    if (r == 0 || r % 2 != 0) return false;
    long long x = Utils::modpow(a_, r / 2, N_);
    if (x == N_ - 1) return false;

    long long f1 = classicalGCD(x + 1, N_);
    long long f2 = classicalGCD(x - 1, N_);

    if (f1 > 1 && f1 < N_) { factor1_ = f1; factor2_ = N_ / f1; return true; }
    if (f2 > 1 && f2 < N_) { factor1_ = f2; factor2_ = N_ / f2; return true; }
    return false;
}

void Shor::buildCircuit() {
    log("Building Shor circuit for N=" + std::to_string(N_) + ", a=" + std::to_string(a_));
    log("Count register: " + std::to_string(countQubits_) + " qubits");
    log("Target register: " + std::to_string(targetQubits_) + " qubits");
}

void Shor::run() {
    // Classical pre-checks
    if (N_ % 2 == 0) { factor1_ = 2; factor2_ = N_ / 2; return; }
    if (isPrime(N_))  { log("N is prime — no factors"); return; }
    if (isPerfectPower(N_)) { log("N is a perfect power"); return; }

    long long g = classicalGCD(a_, N_);
    if (g > 1) { factor1_ = g; factor2_ = N_ / g; return; }

    buildCircuit();

    long long r = findPeriodQuantum();
    if (!extractFactors(r)) {
        log("Period finding failed — try different a");
    }
}

std::string Shor::resultSummary() const {
    std::ostringstream oss;
    oss << "=== Shor's Algorithm ===\n"
        << "  N:        " << N_ << "\n"
        << "  a:        " << a_ << "\n"
        << "  Success:  " << (succeeded() ? "YES" : "NO") << "\n";
    if (succeeded())
        oss << "  Factors:  " << factor1_ << " x " << factor2_ << "\n"
            << "  Verify:   " << factor1_ * factor2_
            << (factor1_ * factor2_ == N_ ? " ✓" : " ✗") << "\n";
    return oss.str();
}

} // namespace QSE