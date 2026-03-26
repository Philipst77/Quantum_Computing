#include "StatevectorSimulator.h"

namespace QSE {

StatevectorSimulator::StatevectorSimulator(int numQubits)
    : numQubits_(numQubits), qreg_(numQubits)
{}

// ── Single-qubit gates ───────────────────────────────────────

StatevectorSimulator& StatevectorSimulator::h(int qubit)
    { Gates::applyH(qreg_, qubit); return *this; }

StatevectorSimulator& StatevectorSimulator::x(int qubit)
    { Gates::applyX(qreg_, qubit); return *this; }

StatevectorSimulator& StatevectorSimulator::y(int qubit)
    { Gates::applyY(qreg_, qubit); return *this; }

StatevectorSimulator& StatevectorSimulator::z(int qubit)
    { Gates::applyZ(qreg_, qubit); return *this; }

StatevectorSimulator& StatevectorSimulator::s(int qubit)
    { Gates::applyS(qreg_, qubit); return *this; }

StatevectorSimulator& StatevectorSimulator::t(int qubit)
    { Gates::applyT(qreg_, qubit); return *this; }

StatevectorSimulator& StatevectorSimulator::sdg(int qubit)
    { Gates::applySdg(qreg_, qubit); return *this; }

StatevectorSimulator& StatevectorSimulator::tdg(int qubit)
    { Gates::applyTdg(qreg_, qubit); return *this; }

StatevectorSimulator& StatevectorSimulator::rx(int qubit, double theta)
    { Gates::applyRx(qreg_, qubit, theta); return *this; }

StatevectorSimulator& StatevectorSimulator::ry(int qubit, double theta)
    { Gates::applyRy(qreg_, qubit, theta); return *this; }

StatevectorSimulator& StatevectorSimulator::rz(int qubit, double theta)
    { Gates::applyRz(qreg_, qubit, theta); return *this; }

StatevectorSimulator& StatevectorSimulator::p(int qubit, double phi)
    { Gates::applyP(qreg_, qubit, phi); return *this; }

StatevectorSimulator& StatevectorSimulator::u(int qubit,
                                               double theta,
                                               double phi,
                                               double lam)
    { Gates::applyU(qreg_, qubit, theta, phi, lam); return *this; }

// ── Two-qubit gates ──────────────────────────────────────────

StatevectorSimulator& StatevectorSimulator::cx(int control, int target)
    { Gates::applyCX(qreg_, control, target); return *this; }

StatevectorSimulator& StatevectorSimulator::cy(int control, int target)
    { Gates::applyCY(qreg_, control, target); return *this; }

StatevectorSimulator& StatevectorSimulator::cz(int control, int target)
    { Gates::applyCZ(qreg_, control, target); return *this; }

StatevectorSimulator& StatevectorSimulator::ch(int control, int target)
    { Gates::applyCH(qreg_, control, target); return *this; }

StatevectorSimulator& StatevectorSimulator::swap(int q1, int q2)
    { Gates::applySWAP(qreg_, q1, q2); return *this; }

StatevectorSimulator& StatevectorSimulator::crz(int control, int target, double theta)
    { Gates::applyCRz(qreg_, control, target, theta); return *this; }

StatevectorSimulator& StatevectorSimulator::cp(int control, int target, double phi)
    { Gates::applyCP(qreg_, control, target, phi); return *this; }

// ── Three-qubit gates ────────────────────────────────────────

StatevectorSimulator& StatevectorSimulator::ccx(int c1, int c2, int target)
    { Gates::applyToffoli(qreg_, c1, c2, target); return *this; }

StatevectorSimulator& StatevectorSimulator::cswap(int control, int t1, int t2)
    { Gates::applyFredkin(qreg_, control, t1, t2); return *this; }

StatevectorSimulator& StatevectorSimulator::applyGate(int qubit,
                                                       const Gates::Matrix2& matrix)
    { qreg_.applyGate(qubit, matrix); return *this; }

// ── Measurement ──────────────────────────────────────────────

int StatevectorSimulator::measure(int qubit)
    { return qreg_.measure(qubit); }

std::vector<int> StatevectorSimulator::measureAll()
    { return qreg_.measureAll(); }

std::map<std::string, int> StatevectorSimulator::sample(int shots)
    { return Utils::sampleCounts(qreg_, shots); }

// ── State access ─────────────────────────────────────────────

double StatevectorSimulator::probability(size_t basisState) const
    { return qreg_.probability(basisState); }

std::vector<double> StatevectorSimulator::probabilities() const
    { return qreg_.probabilities(); }

// ── Utilities ────────────────────────────────────────────────

void StatevectorSimulator::reset()
    { qreg_.reset(); }

void StatevectorSimulator::printState() const
    { Utils::printStateVector(qreg_); }

void StatevectorSimulator::printProbabilities() const
    { Utils::printProbabilities(qreg_); }

std::string StatevectorSimulator::mostLikelyOutcome() const
    { return Utils::mostLikelyOutcome(qreg_); }

} // namespace QSE