#include "QAOA.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <complex>

namespace QSE {

static const std::vector<std::vector<Complex>> H_GATE = {
    {{ 1.0/std::sqrt(2.0), 0}, { 1.0/std::sqrt(2.0), 0}},
    {{ 1.0/std::sqrt(2.0), 0}, {-1.0/std::sqrt(2.0), 0}}
};
static const std::vector<std::vector<Complex>> X_GATE = {
    {{0,0},{1,0}},
    {{1,0},{0,0}}
};

std::vector<std::vector<Complex>> QAOA::rzMatrix(double theta) {
    return {
        {{std::cos(theta/2), -std::sin(theta/2)}, {0,0}},
        {{0,0}, {std::cos(theta/2),  std::sin(theta/2)}}
    };
}

std::vector<std::vector<Complex>> QAOA::rxMatrix(double theta) {
    return {
        {{std::cos(theta/2), 0}, {0, -std::sin(theta/2)}},
        {{0, -std::sin(theta/2)}, {std::cos(theta/2), 0}}
    };
}

QAOA::QAOA(int numQubits, int p, unsigned seed)
    : QuantumAlgorithm("QAOA", numQubits),
      m_p(p), m_gammas(p, 0.5), m_betas(p, 0.5),
      m_bestSolution(0), m_bestCost(-1e18), m_rng(seed) {}

std::string QAOA::description() const {
    return "QAOA with p=" + std::to_string(m_p) + " layers on " +
           std::to_string(numQubits()) + " qubits.";
}

void QAOA::setMaxCutProblem(const std::vector<std::pair<int,int>>& edges) {
    m_edges = edges;
    m_cost = [this](size_t bitstr) -> double {
        double cut = 0;
        for (auto& [u, v] : m_edges) {
            int bu = (bitstr >> u) & 1;
            int bv = (bitstr >> v) & 1;
            if (bu != bv) cut += 1.0;
        }
        return cut;
    };
}

void QAOA::setCostFunction(std::function<double(size_t)> cost) {
    m_cost = std::move(cost);
}

void QAOA::setParameters(const std::vector<double>& gammas,
                          const std::vector<double>& betas) {
    if ((int)gammas.size() != m_p || (int)betas.size() != m_p)
        throw std::invalid_argument("Parameter vectors must have length p");
    m_gammas = gammas;
    m_betas  = betas;
}

void QAOA::prepareInitialState() {
    qreg_.reset();
    for (int q = 0; q < numQubits(); ++q)
        qreg_.applyGate(q, H_GATE);
}

void QAOA::applyProblemUnitary(double gamma) {
    for (auto& [u, v] : m_edges) {
        qreg_.applyControlledGate(u, v, X_GATE);
        qreg_.applyGate(v, rzMatrix(2 * gamma));
        qreg_.applyControlledGate(u, v, X_GATE);
    }
}

void QAOA::applyMixerUnitary(double beta) {
    for (int q = 0; q < numQubits(); ++q)
        qreg_.applyGate(q, rxMatrix(2 * beta));
}

double QAOA::evaluateExpectation() {
    if (!m_cost) throw std::runtime_error("No cost function set");
    auto probs = qreg_.probabilities();
    double E = 0;
    for (size_t i = 0; i < probs.size(); ++i)
        E += probs[i] * m_cost(i);
    return E;
}

double QAOA::expectationValue() {
    prepareInitialState();
    for (int layer = 0; layer < m_p; ++layer) {
        applyProblemUnitary(m_gammas[layer]);
        applyMixerUnitary(m_betas[layer]);
    }
    return evaluateExpectation();
}

void QAOA::optimizeParameters(int maxIter) {
    std::normal_distribution<double> perturb(0, 0.1);
    std::uniform_int_distribution<int> layerD(0, m_p - 1);
    double bestE = -1e18;

    for (int iter = 0; iter < maxIter; ++iter) {
        auto gTry = m_gammas, bTry = m_betas;
        int layer = layerD(m_rng);
        gTry[layer] += perturb(m_rng);
        bTry[layer] += perturb(m_rng);

        auto oldG = m_gammas, oldB = m_betas;
        m_gammas = gTry; m_betas = bTry;
        double E = expectationValue();
        if (E > bestE) { bestE = E; }
        else { m_gammas = oldG; m_betas = oldB; }
    }
}

void QAOA::run() {
    if (!m_cost) throw std::runtime_error("No cost function set.");
    optimizeParameters(200);

    prepareInitialState();
    for (int layer = 0; layer < m_p; ++layer) {
        applyProblemUnitary(m_gammas[layer]);
        applyMixerUnitary(m_betas[layer]);
    }

    auto probs = qreg_.probabilities();
    m_bestCost = -1e18;
    for (size_t i = 0; i < probs.size(); ++i) {
        double c = m_cost(i);
        if (c > m_bestCost) { m_bestCost = c; m_bestSolution = i; }
    }
}

void QAOA::printResult() const {
    std::cout << "=== QAOA ===\n";
    std::cout << "Best solution: " << m_bestSolution
              << " (cost=" << m_bestCost << ")\n";
}

} // namespace QSE