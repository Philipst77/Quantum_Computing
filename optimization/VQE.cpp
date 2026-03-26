#include "VQE.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <complex>

namespace QSE {

static const std::vector<std::vector<Complex>> X_GATE = {
    {{0,0},{1,0}},
    {{1,0},{0,0}}
};
static const std::vector<std::vector<Complex>> Y_GATE = {
    {{0,0},{0,-1}},
    {{0,1},{0,0}}
};
static const std::vector<std::vector<Complex>> Z_GATE = {
    {{1,0},{0,0}},
    {{0,0},{-1,0}}
};

std::vector<std::vector<Complex>> VQE::rxMatrix(double theta) {
    return {
        {{std::cos(theta/2), 0}, {0, -std::sin(theta/2)}},
        {{0, -std::sin(theta/2)}, {std::cos(theta/2), 0}}
    };
}

std::vector<std::vector<Complex>> VQE::ryMatrix(double theta) {
    return {
        {{std::cos(theta/2), 0}, {-std::sin(theta/2), 0}},
        {{std::sin(theta/2), 0}, { std::cos(theta/2), 0}}
    };
}

std::vector<std::vector<Complex>> VQE::rzMatrix(double theta) {
    return {
        {{std::cos(theta/2), -std::sin(theta/2)}, {0,0}},
        {{0,0}, {std::cos(theta/2), std::sin(theta/2)}}
    };
}

VQE::VQE(int numQubits, int layers, unsigned seed)
    : QuantumAlgorithm("VQE", numQubits),
      m_layers(layers),
      m_params(numQubits * layers * 3, 0.1),
      m_groundEnergy(1e18),
      m_optIters(0),
      m_rng(seed)
{
    std::uniform_real_distribution<double> d(0, 2 * M_PI);
    for (auto& p : m_params) p = d(m_rng);
}

std::string VQE::description() const {
    return "VQE with " + std::to_string(m_layers) + "-layer ansatz on " +
           std::to_string(numQubits()) + " qubits.";
}

void VQE::setHamiltonian(const std::vector<PauliTerm>& H) {
    m_hamiltonian = H;
}

void VQE::setH2Hamiltonian() {
    m_hamiltonian = {
        {-0.8105479, "IIII"},
        { 0.1721839, "ZZII"},
        {-0.2257534, "IIZI"},
        { 0.1209362, "IIZZ"},
        {-0.2257534, "IZII"},
        { 0.1686890, "ZZII"},
        { 0.1209362, "ZIZI"},
        { 0.0453220, "YYYY"},
        { 0.0453220, "XXYY"},
        { 0.0453220, "YYXX"},
        { 0.0453220, "XXXX"}
    };
}

void VQE::setParameters(const std::vector<double>& params) {
    m_params = params;
}

void VQE::applyAnsatz() {
    int n = numQubits();
    qreg_.reset();
    int idx = 0;
    for (int layer = 0; layer < m_layers; ++layer) {
        for (int q = 0; q < n; ++q) {
            qreg_.applyGate(q, rxMatrix(m_params[idx++]));
            qreg_.applyGate(q, ryMatrix(m_params[idx++]));
            qreg_.applyGate(q, rzMatrix(m_params[idx++]));
        }
        for (int q = 0; q < n - 1; ++q)
            qreg_.applyControlledGate(q, q+1, X_GATE);
        if (n > 2)
            qreg_.applyControlledGate(n-1, 0, X_GATE);
    }
}

double VQE::measurePauliString(const std::string& pauliStr) {
    QubitRegister tmp = qreg_;
    int n = numQubits();
    for (int q = 0; q < n && q < (int)pauliStr.size(); ++q) {
        char p = pauliStr[pauliStr.size() - 1 - q];
        if      (p == 'X') tmp.applyGate(q, X_GATE);
        else if (p == 'Y') tmp.applyGate(q, Y_GATE);
        else if (p == 'Z') tmp.applyGate(q, Z_GATE);
    }
    Complex exp_val = 0;
    for (size_t i = 0; i < qreg_.size(); ++i)
        exp_val += std::conj(qreg_[i]) * tmp[i];
    return exp_val.real();
}

double VQE::energyExpectation() {
    applyAnsatz();
    double E = 0;
    for (auto& [coeff, pauli] : m_hamiltonian)
        E += coeff * measurePauliString(pauli);
    return E;
}

double VQE::paramShiftGradient(int idx) {
    const double shift = M_PI / 2.0;
    m_params[idx] += shift;
    double Ep = energyExpectation();
    m_params[idx] -= 2 * shift;
    double Em = energyExpectation();
    m_params[idx] += shift;
    return (Ep - Em) / 2.0;
}

void VQE::gradientDescent(int maxIter, double lr) {
    for (int iter = 0; iter < maxIter; ++iter) {
        double E = energyExpectation();
        if (E < m_groundEnergy) m_groundEnergy = E;
        std::vector<double> grads(m_params.size());
        for (size_t i = 0; i < m_params.size(); ++i)
            grads[i] = paramShiftGradient(i);
        for (size_t i = 0; i < m_params.size(); ++i)
            m_params[i] -= lr * grads[i];
        m_optIters++;
    }
}

void VQE::run() {
    if (m_hamiltonian.empty()) setH2Hamiltonian();
    gradientDescent(50, 0.1);
    m_groundEnergy = energyExpectation();
}

void VQE::printResult() const {
    std::cout << "=== VQE ===\n";
    std::cout << "Ground energy estimate : " << m_groundEnergy << " Ha\n";
    std::cout << "Optimisation iterations: " << m_optIters << "\n";
}

} // namespace QSE