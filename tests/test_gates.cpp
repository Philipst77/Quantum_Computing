// test_gates.cpp — Tests for all gate types, parametric gates, and StatevectorSimulator
#include "../core/QubitRegister.h"
#include "../core/QuantumGate.h"
#include "../simulators/statevector/StatevectorSimulator.h"
#include <iostream>
#include <cmath>
#include <vector>

using namespace QSE;

int passed = 0, failed = 0;

#define TEST(name, expr) do { \
    bool _ok = (expr); \
    if (_ok) { std::cout << "  [PASS] " << name << "\n"; ++passed; } \
    else      { std::cout << "  [FAIL] " << name << "\n"; ++failed; } \
} while(0)

#define NEAR(a, b)  (std::abs((double)(a) - (double)(b)) < 1e-9)
#define CNEAR(a, b) (std::abs((a) - (b)) < 1e-9)

// ── Rotation gates ───────────────────────────────────────────

void testRotationGates() {
    std::cout << "\n-- Rotation gates --\n";

    // Rx(2π) = -I  (global phase)
    QubitRegister q(1);
    Gates::applyRx(q, 0, 2 * M_PI);
    TEST("Rx(2π)|0> amplitude magnitude = 1", NEAR(std::abs(q[0]), 1.0));

    // Ry(π)|0> = |1>
    QubitRegister q2(1);
    Gates::applyRy(q2, 0, M_PI);
    TEST("Ry(π)|0>: |1> amplitude ≈ 1", NEAR(std::abs(q2[1]), 1.0));
    TEST("Ry(π)|0>: |0> amplitude ≈ 0", NEAR(std::abs(q2[0]), 0.0));

    // Ry(π/2)|0> = (|0>+|1>)/√2
    QubitRegister q3(1);
    Gates::applyRy(q3, 0, M_PI / 2);
    TEST("Ry(π/2)|0>: p(0) = 0.5", NEAR(q3.probability(0), 0.5));
    TEST("Ry(π/2)|0>: p(1) = 0.5", NEAR(q3.probability(1), 0.5));
    TEST("Ry(π/2)|0>: normalized",  q3.isNormalized());

    // Rz(π)|+> = |->  (relative phase flip)
    QubitRegister q4(1);
    Gates::applyH(q4, 0);
    Gates::applyRz(q4, 0, M_PI);
    TEST("Rz(π)|+>[0] ≈ 1/√2",  NEAR(std::abs(q4[0]), M_SQRT1_2));
    TEST("Rz(π)|+>[1] ≈ -1/√2", NEAR(q4[1].real(), -M_SQRT1_2));

    // Rx and Ry self-inverse at angle π (up to global phase)
    QubitRegister q5(1);
    Gates::applyRx(q5, 0, M_PI / 4);
    Gates::applyRx(q5, 0, M_PI / 4);
    Gates::applyRx(q5, 0, -M_PI / 2);
    TEST("Rx(π/4)·Rx(π/4)·Rx(-π/2) = I", NEAR(std::abs(q5[0]), 1.0));
}

// ── Phase gates ──────────────────────────────────────────────

void testPhaseGates() {
    std::cout << "\n-- Phase gates --\n";

    // S|0> = |0>  (phase only affects |1>)
    QubitRegister q(1);
    Gates::applyS(q, 0);
    TEST("S|0> = |0>", NEAR(std::abs(q[0]), 1.0));

    // S|1> = i|1>
    QubitRegister q2(1);
    Gates::applyX(q2, 0);
    Gates::applyS(q2, 0);
    TEST("S|1>: real = 0", NEAR(q2[1].real(), 0.0));
    TEST("S|1>: imag = 1", NEAR(q2[1].imag(), 1.0));

    // S² = Z
    QubitRegister q3(1), q4(1);
    Gates::applyX(q3, 0); Gates::applyS(q3, 0); Gates::applyS(q3, 0);
    Gates::applyX(q4, 0); Gates::applyZ(q4, 0);
    TEST("S²|1> = Z|1>", CNEAR(q3[1], q4[1]));

    // T|1> = e^{iπ/4}|1>
    QubitRegister q5(1);
    Gates::applyX(q5, 0);
    Gates::applyT(q5, 0);
    TEST("T|1>: real = 1/√2", NEAR(q5[1].real(), M_SQRT1_2));
    TEST("T|1>: imag = 1/√2", NEAR(q5[1].imag(), M_SQRT1_2));

    // T² = S
    QubitRegister q6(1);
    Gates::applyX(q6, 0);
    Gates::applyT(q6, 0); Gates::applyT(q6, 0);
    QubitRegister q7(1);
    Gates::applyX(q7, 0); Gates::applyS(q7, 0);
    TEST("T²|1> = S|1>: imag equal", NEAR(q6[1].imag(), q7[1].imag()));

    // P(φ) gate
    double phi = M_PI / 3;
    QubitRegister q8(1);
    Gates::applyX(q8, 0);
    Gates::applyP(q8, 0, phi);
    TEST("P(π/3)|1>: real = cos(π/3)", NEAR(q8[1].real(), std::cos(phi)));
    TEST("P(π/3)|1>: imag = sin(π/3)", NEAR(q8[1].imag(), std::sin(phi)));
}

// ── Controlled gates ─────────────────────────────────────────

void testControlledGates() {
    std::cout << "\n-- Controlled gates --\n";

    // CRz: only applies Rz when control = |1>
    QubitRegister q(2);
    Gates::applyX(q, 0);
    Gates::applyH(q, 1);
    double theta = M_PI;
    Gates::applyCRz(q, 0, 1, theta);
    TEST("CRz(π) with control=|1>: normalized", q.isNormalized());

    // CP gate: controlled phase
    QubitRegister q2(2);
    Gates::applyX(q2, 0); Gates::applyX(q2, 1); // |11>
    Gates::applyCP(q2, 0, 1, M_PI);
    TEST("CP(π)|11> phase = -1", NEAR(q2[3].real(), -1.0));

    // CY gate
    QubitRegister q3(2);
    Gates::applyX(q3, 0); // |01> -> index 1
    Gates::applyCY(q3, 0, 1);
    // CY|01>: control=q0=1, Y applied to q1 -> i|11> = index 3
    TEST("CY|10>: target flipped", NEAR(std::abs(q3[3]), 1.0));
    TEST("CY|10>: phase is i", NEAR(q3[3].imag(), 1.0));
}

// ── Multi-qubit entanglement ──────────────────────────────────

void testEntanglement() {
    std::cout << "\n-- Entanglement --\n";

    // GHZ state: H(0) + CX(0,1) + CX(0,2)
    QubitRegister q(3);
    Gates::applyH(q, 0);
    Gates::applyCX(q, 0, 1);
    Gates::applyCX(q, 0, 2);
    TEST("GHZ: p(|000>) = 0.5", NEAR(q.probability(0), 0.5));
    TEST("GHZ: p(|111>) = 0.5", NEAR(q.probability(7), 0.5));
    TEST("GHZ: p(|001>) = 0.0", NEAR(q.probability(1), 0.0));
    TEST("GHZ: p(|110>) = 0.0", NEAR(q.probability(6), 0.0));
    TEST("GHZ: normalized", q.isNormalized());

    // W state: (|100> + |010> + |001>) / √3
    QubitRegister w(3);
    StateVector sv(8, {0,0});
    double amp = 1.0 / std::sqrt(3.0);
    sv[1] = sv[2] = sv[4] = {amp, 0.0};
    w.setState(sv);
    TEST("W state: p(|100>) = 1/3", NEAR(w.probability(4), 1.0/3.0));
    TEST("W state: p(|010>) = 1/3", NEAR(w.probability(2), 1.0/3.0));
    TEST("W state: p(|001>) = 1/3", NEAR(w.probability(1), 1.0/3.0));
    TEST("W state: normalized", w.isNormalized());
}

// ── StatevectorSimulator fluent API ──────────────────────────

void testStatevectorSimulator() {
    std::cout << "\n-- StatevectorSimulator (fluent API) --\n";

    // Bell state
    StatevectorSimulator sim(2);
    sim.h(0).cx(0, 1);
    auto probs = sim.probabilities();
    TEST("Bell: p(|00>) = 0.5", NEAR(probs[0], 0.5));
    TEST("Bell: p(|11>) = 0.5", NEAR(probs[3], 0.5));
    TEST("Bell: p(|01>) = 0.0", NEAR(probs[1], 0.0));
    TEST("Bell: p(|10>) = 0.0", NEAR(probs[2], 0.0));

    // GHZ state via simulator
    StatevectorSimulator sim2(3);
    sim2.h(0).cx(0,1).cx(0,2);
    auto p2 = sim2.probabilities();
    TEST("GHZ via sim: p(|000>) = 0.5", NEAR(p2[0], 0.5));
    TEST("GHZ via sim: p(|111>) = 0.5", NEAR(p2[7], 0.5));

    // Rotation circuit
    StatevectorSimulator sim3(1);
    sim3.ry(0, M_PI / 2);
    auto p3 = sim3.probabilities();
    TEST("Ry(π/2) via sim: p(0)=0.5", NEAR(p3[0], 0.5));
    TEST("Ry(π/2) via sim: p(1)=0.5", NEAR(p3[1], 0.5));

    // Toffoli (CCX) via simulator
    StatevectorSimulator sim4(3);
    sim4.x(0).x(1).ccx(0, 1, 2); // qubits 0,1 set -> index 3, after CCX -> index 7
    TEST("CCX via sim: |110>->|111>", NEAR(sim4.probabilities()[7], 1.0));

    // SWAP via simulator
    // x(0) sets qubit 0 (LSB) -> index 1
    // swap(0,1) moves amplitude to qubit 1 (LSB) -> index 2
    StatevectorSimulator sim5(2);
    sim5.x(0).swap(0, 1);
    TEST("SWAP via sim: |10>->|01>", NEAR(sim5.probabilities()[2], 1.0));  // FIXED: index 2

    // mostLikelyOutcome
    StatevectorSimulator sim6(2);
    sim6.x(1); // |10>
    TEST("mostLikelyOutcome |10> = 10", sim6.mostLikelyOutcome() == "10");

    // Reset
    sim6.reset();
    TEST("After reset: p(|00>) = 1.0", NEAR(sim6.probabilities()[0], 1.0));

    // Sampling returns correct keys
    StatevectorSimulator sim7(2);
    sim7.h(0).cx(0, 1); // Bell state
    auto counts = sim7.sample(200);
    TEST("Bell sample: only 00 and 11 appear",
         counts.count("01") == 0 && counts.count("10") == 0);
    TEST("Bell sample: both 00 and 11 seen", counts.count("00") > 0 && counts.count("11") > 0);
}

// ── Tensor product ───────────────────────────────────────────

void testTensorProduct() {
    std::cout << "\n-- Tensor product --\n";

    QubitRegister q0(1), q1(1);
    Gates::applyX(q1, 0); // q1 = |1>
    auto combined = QubitRegister::tensorProduct(q0, q1);
    TEST("|0>⊗|1> = |01>: amp[1]=1", NEAR(std::abs(combined[1]), 1.0));
    TEST("|0>⊗|1> = |01>: amp[0]=0", NEAR(std::abs(combined[0]), 0.0));
    TEST("|0>⊗|1>: 2 qubits", combined.numQubits() == 2);
}

int main() {
    std::cout << "=== QSE Gate Tests ===\n";
    testRotationGates();
    testPhaseGates();
    testControlledGates();
    testEntanglement();
    testStatevectorSimulator();
    testTensorProduct();

    std::cout << "\n=== Results: " << passed << " passed, " << failed << " failed ===\n";
    return (failed > 0) ? 1 : 0;
}