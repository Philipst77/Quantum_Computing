#pragma once
#include "../../core/QubitRegister.h"
#include "../../core/QuantumGate.h"
#include "../../core/Utils.h"
#include <map>
#include <string>
#include <vector>

namespace QSE {

// Full statevector simulator — exact simulation of n-qubit systems.
// Memory: O(2^n) complex doubles. Practical limit ~28 qubits.
// Fluent interface: sim.h(0).cx(0,1).measure(0)
class StatevectorSimulator {
public:
    explicit StatevectorSimulator(int numQubits);

    // ── Single-qubit gates ───────────────────────────────────
    StatevectorSimulator& h(int qubit);
    StatevectorSimulator& x(int qubit);
    StatevectorSimulator& y(int qubit);
    StatevectorSimulator& z(int qubit);
    StatevectorSimulator& s(int qubit);
    StatevectorSimulator& t(int qubit);
    StatevectorSimulator& sdg(int qubit);
    StatevectorSimulator& tdg(int qubit);
    StatevectorSimulator& rx(int qubit, double theta);
    StatevectorSimulator& ry(int qubit, double theta);
    StatevectorSimulator& rz(int qubit, double theta);
    StatevectorSimulator& p(int qubit, double phi);
    StatevectorSimulator& u(int qubit, double theta, double phi, double lam);

    // ── Two-qubit gates ──────────────────────────────────────
    StatevectorSimulator& cx(int control, int target);
    StatevectorSimulator& cy(int control, int target);
    StatevectorSimulator& cz(int control, int target);
    StatevectorSimulator& ch(int control, int target);
    StatevectorSimulator& swap(int q1, int q2);
    StatevectorSimulator& crz(int control, int target, double theta);
    StatevectorSimulator& cp(int control, int target, double phi);

    // ── Three-qubit gates ────────────────────────────────────
    StatevectorSimulator& ccx(int c1, int c2, int target);   // Toffoli
    StatevectorSimulator& cswap(int control, int t1, int t2); // Fredkin

    // ── Custom gate ──────────────────────────────────────────
    StatevectorSimulator& applyGate(int qubit, const Gates::Matrix2& matrix);

    // ── Measurement ──────────────────────────────────────────
    int measure(int qubit);
    std::vector<int> measureAll();
    std::map<std::string, int> sample(int shots = 1024);

    // ── State access ─────────────────────────────────────────
    const QubitRegister&    qreg()          const { return qreg_; }
    const StateVector&      stateVector()   const { return qreg_.state(); }
    double                  probability(size_t basisState) const;
    std::vector<double>     probabilities() const;

    // ── Utilities ────────────────────────────────────────────
    void        reset();
    void        printState()         const;
    void        printProbabilities() const;
    std::string mostLikelyOutcome()  const;
    int         numQubits()          const { return numQubits_; }

private:
    int            numQubits_;
    QubitRegister  qreg_;
};

} // namespace QSE