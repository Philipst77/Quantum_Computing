#pragma once
#include "../core/QuantumAlgorithm.h"
#include <vector>
#include <functional>
#include <random>
#include <string>

namespace QSE {

class QAOA : public QuantumAlgorithm {
public:
    QAOA(int numQubits, int p, unsigned seed = 42);

    void buildCircuit() override {}
    void run() override;
    std::string resultSummary() const override { return description(); }
    std::string description() const;

    void setMaxCutProblem(const std::vector<std::pair<int,int>>& edges);
    void setCostFunction(std::function<double(size_t)> cost);
    void setParameters(const std::vector<double>& gammas,
                       const std::vector<double>& betas);
    void optimizeParameters(int maxIter = 100);

    size_t              getBestSolution() const { return m_bestSolution; }
    double              getBestCost()     const { return m_bestCost; }
    std::vector<double> getGammas()       const { return m_gammas; }
    std::vector<double> getBetas()        const { return m_betas; }

    double expectationValue();
    void   printResult() const;

private:
    int                              m_p;
    std::vector<double>              m_gammas, m_betas;
    std::function<double(size_t)>    m_cost;
    std::vector<std::pair<int,int>>  m_edges;
    size_t                           m_bestSolution;
    double                           m_bestCost;
    std::mt19937                     m_rng;

    void prepareInitialState();
    void applyProblemUnitary(double gamma);
    void applyMixerUnitary(double beta);
    double evaluateExpectation();

    // Gate helpers using qreg_.applyGate / applyControlledGate
    static std::vector<std::vector<Complex>> rzMatrix(double theta);
    static std::vector<std::vector<Complex>> rxMatrix(double theta);
};

} // namespace QSE