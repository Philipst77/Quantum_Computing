// test_core.cpp — Unit tests for QubitRegister, QuantumGate, and Utils
#include "../core/QubitRegister.h"
#include "../core/QuantumGate.h"
#include "../core/Utils.h"
#include <iostream>
#include <cmath>

using namespace QSE;

int passed = 0, failed = 0;

#define TEST(name, expr) do { \
    bool _ok = (expr); \
    if (_ok) { std::cout << "  [PASS] " << name << "\n"; ++passed; } \
    else      { std::cout << "  [FAIL] " << name << "\n"; ++failed; } \
} while(0)

static inline bool near(double a, double b)        { return std::abs(a - b) < 1e-9; }
static inline bool ceq(Complex a, Complex b)       { return std::abs(a - b) < 1e-9; }

// ── QubitRegister ────────────────────────────────────────────

void testQubitRegister() {
    std::cout << "\n-- QubitRegister --\n";

    QubitRegister q1(1);
    TEST("1-qubit: starts in |0>",       ceq(q1[0], {1,0}) && ceq(q1[1], {0,0}));
    TEST("1-qubit: size = 2",            q1.size() == 2);
    TEST("1-qubit: normalized",          q1.isNormalized());

    QubitRegister q3(3);
    TEST("3-qubit: size = 8",            q3.size() == 8);
    TEST("3-qubit: numQubits = 3",       q3.numQubits() == 3);
    TEST("3-qubit: amplitude |000> = 1", ceq(q3[0], {1,0}));

    q3.setComputationalBasis(5);
    TEST("setComputationalBasis(5): amp[5]=1", ceq(q3[5], {1,0}));
    TEST("setComputationalBasis(5): amp[0]=0", ceq(q3[0], {0,0}));
    TEST("setComputationalBasis normalized",   q3.isNormalized());

    q3.reset();
    TEST("reset: back to |000>", ceq(q3[0], {1,0}));

    StateVector sv(8, {0,0});
    sv[7] = {1,0};
    q3.setState(sv);
    TEST("setState: amp[7]=1", ceq(q3[7], {1,0}));

    QubitRegister q2(1);
    TEST("norm of |0> = 1.0",           near(q2.norm(), 1.0));
    TEST("probability(0) of |0> = 1.0", near(q2.probability(0), 1.0));
    TEST("probability(1) of |0> = 0.0", near(q2.probability(1), 0.0));
}

// ── Gate matrices ────────────────────────────────────────────

void testGateMatrices() {
    std::cout << "\n-- Gate matrices --\n";

    TEST("H is unitary",        Gates::isUnitary(Gates::H));
    TEST("X is unitary",        Gates::isUnitary(Gates::X));
    TEST("Y is unitary",        Gates::isUnitary(Gates::Y));
    TEST("Z is unitary",        Gates::isUnitary(Gates::Z));
    TEST("S is unitary",        Gates::isUnitary(Gates::S));
    TEST("T is unitary",        Gates::isUnitary(Gates::T));
    TEST("Sdg is unitary",      Gates::isUnitary(Gates::Sdg));
    TEST("Tdg is unitary",      Gates::isUnitary(Gates::Tdg));
    TEST("Rx(π/3) is unitary",  Gates::isUnitary(Gates::Rx(M_PI / 3)));
    TEST("Ry(π/4) is unitary",  Gates::isUnitary(Gates::Ry(M_PI / 4)));
    TEST("Rz(2π/5) is unitary", Gates::isUnitary(Gates::Rz(2 * M_PI / 5)));

    TEST("X is Hermitian",     Gates::isHermitian(Gates::X));
    TEST("Y is Hermitian",     Gates::isHermitian(Gates::Y));
    TEST("Z is Hermitian",     Gates::isHermitian(Gates::Z));
    TEST("H is Hermitian",     Gates::isHermitian(Gates::H));
    TEST("S is NOT Hermitian", !Gates::isHermitian(Gates::S));

    auto XtX = Gates::matMul(Gates::X, Gates::X);
    TEST("X*X = I", ceq(XtX[0][0], {1,0}) && ceq(XtX[1][1], {1,0})
                 && ceq(XtX[0][1], {0,0}) && ceq(XtX[1][0], {0,0}));

    auto HtH = Gates::matMul(Gates::H, Gates::H);
    TEST("H*H = I", ceq(HtH[0][0], {1,0}) && ceq(HtH[1][1], {1,0}));
}

// ── Single-qubit gates ───────────────────────────────────────

void testSingleQubitGates() {
    std::cout << "\n-- Single-qubit gate application --\n";

    QubitRegister q(1);
    Gates::applyX(q, 0);
    TEST("X|0> = |1>", ceq(q[0], {0,0}) && near(std::abs(q[1]), 1.0));

    Gates::applyX(q, 0);
    TEST("X(X|0>) = |0>", near(std::abs(q[0]), 1.0) && ceq(q[1], {0,0}));

    QubitRegister q2(1);
    Gates::applyH(q2, 0);
    TEST("H|0>[0] = 1/√2",  near(q2[0].real(), M_SQRT1_2));
    TEST("H|0>[1] = 1/√2",  near(q2[1].real(), M_SQRT1_2));
    TEST("H|0> normalized", q2.isNormalized());

    Gates::applyH(q2, 0);
    TEST("H(H|0>) = |0>", near(std::abs(q2[0]), 1.0));

    QubitRegister q3(1);
    Gates::applyH(q3, 0);
    Gates::applyZ(q3, 0);
    TEST("Z|+>[0] =  1/√2", near(q3[0].real(),  M_SQRT1_2));
    TEST("Z|+>[1] = -1/√2", near(q3[1].real(), -M_SQRT1_2));

    QubitRegister q4(1);
    Gates::applyY(q4, 0);
    TEST("Y|0>[0] = 0", ceq(q4[0], {0,0}));
    TEST("Y|0>[1] = i", near(q4[1].imag(), 1.0));

    QubitRegister q5(1);
    Gates::applyX(q5, 0);
    Gates::applyS(q5, 0);
    TEST("S|1>: real=0", near(q5[1].real(), 0.0));
    TEST("S|1>: imag=1", near(q5[1].imag(), 1.0));

    Gates::applySdg(q5, 0);
    TEST("Sdg(S|1>) = |1>", near(std::abs(q5[1]), 1.0) && near(q5[1].imag(), 0.0));

    QubitRegister q6(1);
    Gates::applyX(q6, 0);
    Gates::applyT(q6, 0);
    TEST("T|1>: real = 1/√2", near(q6[1].real(), M_SQRT1_2));
    TEST("T|1>: imag = 1/√2", near(q6[1].imag(), M_SQRT1_2));

    QubitRegister q7(1);
    Gates::applyRx(q7, 0, M_PI);
    TEST("Rx(π)|0>: |amp[1]| ≈ 1", near(std::abs(q7[1]), 1.0));

    QubitRegister q8(1);
    Gates::applyRy(q8, 0, M_PI / 2);
    TEST("Ry(π/2)|0> normalized",   q8.isNormalized());
    TEST("Ry(π/2)|0>: both nonzero",
         std::abs(q8[0]) > 0.1 && std::abs(q8[1]) > 0.1);
}

// ── Two-qubit gates ──────────────────────────────────────────

void testTwoQubitGates() {
    std::cout << "\n-- Two-qubit gates --\n";

    // qubit 0 = LSB, so applyX(q,0) sets index 1 (0b01),
    // applyX(q,1) sets index 2 (0b10)

    QubitRegister q(2);
    Gates::applyX(q, 0);               // |01> -> index 1
    Gates::applyCX(q, 0, 1);           // control=q0=1, flip q1 -> |11> = index 3
    TEST("CX|10> = |11>", near(std::abs(q[3]), 1.0));

    QubitRegister q2(2);
    Gates::applyCX(q2, 0, 1);          // |00>, control=0, no change
    TEST("CX|00> = |00>", near(std::abs(q2[0]), 1.0));

    QubitRegister q3(2);
    Gates::applyX(q3, 1);              // |10> -> index 2 (qubit1=1, qubit0=0)
    Gates::applyCX(q3, 0, 1);          // control=q0=0, no change
    TEST("CX|01> unchanged", near(std::abs(q3[2]), 1.0));  // FIXED: index 2

    QubitRegister q4(2);
    Gates::applyX(q4, 0); Gates::applyX(q4, 1);  // |11> = index 3
    Gates::applyCZ(q4, 0, 1);
    TEST("CZ|11> phase = -1", near(q4[3].real(), -1.0));

    QubitRegister q5(2);
    Gates::applyX(q5, 0);              // |01> -> index 1
    Gates::applySWAP(q5, 0, 1);        // swap -> |10> = index 2
    TEST("SWAP|10> = |01>", near(std::abs(q5[2]), 1.0));  // FIXED: index 2

    QubitRegister q6(2);
    Gates::applyH(q6, 0);
    Gates::applyCX(q6, 0, 1);
    TEST("Bell p(|00>) = 0.5", near(q6.probability(0), 0.5));
    TEST("Bell p(|11>) = 0.5", near(q6.probability(3), 0.5));
    TEST("Bell p(|01>) = 0.0", near(q6.probability(1), 0.0));
}

// ── Three-qubit gates ────────────────────────────────────────

void testThreeQubitGates() {
    std::cout << "\n-- Three-qubit gates --\n";

    // qubit 0 = LSB. applyX(q,0) -> bit0=1, applyX(q,1) -> bit1=1, applyX(q,2) -> bit2=1
    // |110> means qubit2=1,qubit1=1,qubit0=0 -> index = 0b110 = 6... but
    // with LSB convention: applyX(q,0) and applyX(q,1) -> index = 0b011 = 3... 
    // Let's be precise: index = sum of (1 << qubitN) for each set qubit

    QubitRegister q(3);
    Gates::applyX(q, 0); Gates::applyX(q, 1);  // index = 1+2 = 3
    Gates::applyToffoli(q, 0, 1, 2);            // both controls set, flip q2 -> index = 1+2+4 = 7
    TEST("Toffoli|110> = |111>", near(std::abs(q[7]), 1.0));

    QubitRegister q2(3);
    Gates::applyX(q2, 0);                       // index = 1 (only qubit0 set)
    Gates::applyToffoli(q2, 0, 1, 2);           // control q1=0, no change -> stays at index 1
    TEST("Toffoli|100> unchanged", near(std::abs(q2[1]), 1.0));  // FIXED: index 1

    QubitRegister q3(3);
    Gates::applyX(q3, 1);                       // index = 2 (only qubit1 set)
    Gates::applyToffoli(q3, 0, 1, 2);           // control q0=0, no change -> stays at index 2
    TEST("Toffoli|010> unchanged", near(std::abs(q3[2]), 1.0));

    QubitRegister q4(3);
    Gates::applyX(q4, 1);                       // index = 2
    Gates::applyFredkin(q4, 0, 1, 2);           // control q0=0, no change
    TEST("Fredkin|010> ctrl=0: unchanged", near(std::abs(q4[2]), 1.0));

    QubitRegister q5(3);
    Gates::applyX(q5, 0); Gates::applyX(q5, 1); // index = 3
    Gates::applyFredkin(q5, 0, 1, 2);            // control q0=1, swap q1 and q2 -> index = 1+4 = 5
    TEST("Fredkin|110> ctrl=1: |101>", near(std::abs(q5[5]), 1.0));
}

// ── Measurement ──────────────────────────────────────────────

void testMeasurement() {
    std::cout << "\n-- Measurement --\n";

    QubitRegister q(1);
    TEST("Measure |0> = 0", q.measure(0) == 0);

    QubitRegister q2(1);
    Gates::applyX(q2, 0);
    TEST("Measure |1> = 1", q2.measure(0) == 1);

    QubitRegister q3(1);
    Gates::applyH(q3, 0);
    int m = q3.measure(0);
    TEST("Measure |+>: 0 or 1",            m == 0 || m == 1);
    TEST("State normalized after measure",  q3.isNormalized());
    TEST("State definite after measure",
         near(q3.probability(m), 1.0) && near(q3.probability(1 - m), 0.0));

    QubitRegister q4(3);
    auto bits = q4.measureAll();
    TEST("measureAll returns 3 bits",  (int)bits.size() == 3);
    TEST("measureAll |000> = 000",     bits[0]==0 && bits[1]==0 && bits[2]==0);
}

// ── Utils ────────────────────────────────────────────────────

void testUtils() {
    std::cout << "\n-- Utils --\n";

    TEST("toBitstring(0,4)  = 0000",  Utils::toBitstring(0,  4) == "0000");
    TEST("toBitstring(5,4)  = 0101",  Utils::toBitstring(5,  4) == "0101");
    TEST("toBitstring(15,4) = 1111",  Utils::toBitstring(15, 4) == "1111");
    TEST("fromBitstring(0101) = 5",   Utils::fromBitstring("0101") == 5);
    TEST("fromBitstring(1111) = 15",  Utils::fromBitstring("1111") == 15);

    auto bits = Utils::intToBits(13, 4);
    TEST("intToBits(13,4)[0] = 1", bits[0] == 1);
    TEST("intToBits(13,4)[3] = 1", bits[3] == 1);
    TEST("bitsToInt roundtrip",    (int)Utils::bitsToInt(bits) == 13);

    TEST("modpow(2,10,1000) = 24", Utils::modpow(2, 10, 1000) == 24);
    TEST("modpow(3,0,7) = 1",      Utils::modpow(3, 0, 7) == 1);
    TEST("gcd(48,18) = 6",         Utils::gcd(48, 18) == 6);
    TEST("gcd(100,75) = 25",       Utils::gcd(100, 75) == 25);
    TEST("popcount(0) = 0",        Utils::popcount(0) == 0);
    TEST("popcount(7) = 3",        Utils::popcount(7) == 3);
    TEST("popcount(255) = 8",      Utils::popcount(255) == 8);

    QubitRegister a(2), b(2);
    TEST("Fidelity identical = 1.0", near(Utils::fidelity(a, b), 1.0));

    QubitRegister c(1), d(1);
    Gates::applyX(d, 0);
    TEST("Fidelity orthogonal = 0.0", near(Utils::fidelity(c, d), 0.0));

    QubitRegister e(2);
    Gates::applyX(e, 1);
    TEST("mostLikelyOutcome |10> = 10", Utils::mostLikelyOutcome(e) == "10");
}

int main() {
    std::cout << "=== QSE Core Tests ===\n";
    testQubitRegister();
    testGateMatrices();
    testSingleQubitGates();
    testTwoQubitGates();
    testThreeQubitGates();
    testMeasurement();
    testUtils();

    std::cout << "\n=== Results: " << passed << " passed, " << failed << " failed ===\n";
    return (failed > 0) ? 1 : 0;
}