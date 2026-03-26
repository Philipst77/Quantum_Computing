#pragma once
#include "../core/QuantumAlgorithm.h"
#include <complex>
#include <utility>
#include <string>
#include <random>

namespace QSE {

class ThreeQubitFlip : public QuantumAlgorithm {
public:
    explicit ThreeQubitFlip(unsigned seed = 42);

    void buildCircuit() override {}
    void run() override;
    std::string resultSummary() const override { return description(); }
    std::string description() const;

    void encode(Complex alpha, Complex beta);
    void injectError(int qubit);
    void detectAndCorrect();
    std::pair<Complex, Complex> decode();

    bool success() const { return m_success; }
    void printResult() const;

private:
    bool    m_success;
    Complex m_alpha, m_beta;
    std::mt19937 m_rng;

    int syndromeZ12();
    int syndromeZ02();
};

} // namespace QSE