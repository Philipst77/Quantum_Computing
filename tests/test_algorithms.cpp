// test_algorithms.cpp — Tests for Grover, QFT, Teleportation, BB84,
//                       QAOA, ThreeQubitFlip, ShorCode
// NOTE: VQE tests intentionally omitted until real VQE.h is provided.

#include "../algorithms/grover/Grover.h"
#include "../algorithms/qft/QFT.h"
#include "../algorithms/teleportation/Teleportation.h"
#include "../algorithms/cryptography/BB84.h"
#include "../optimization/QAOA.h"
#include "../error_correction/ThreeQubitFlip.h"
#include "../error_correction/ShorCode.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <string>

using namespace QSE;

int passed = 0, failed = 0;

#define TEST(name, expr) do { \
    bool _ok = (expr); \
    if (_ok) { std::cout << "  [PASS] " << name << "\n"; ++passed; } \
    else      { std::cout << "  [FAIL] " << name << "\n"; ++failed; } \
} while(0)

static inline bool near(double a, double b, double eps = 1e-9) {
    return std::abs(a - b) < eps;
}

static inline bool ceq(Complex a, Complex b, double eps = 1e-9) {
    return std::abs(a - b) < eps;
}

// ── Grover's search ──────────────────────────────────────────

void testGrover() {
    std::cout << "\n-- Grover's search --\n";

    for (int target : {0, 1, 7, 10, 15}) {
        Grover g(4, [target](size_t s){ return (int)s == target; });
        g.run();
        TEST(("4-qubit search: find target=" + std::to_string(target)).c_str(),
             g.success());
    }

    TEST("optimalIterations(4) = 3",  Grover::optimalIterations(4) == 3);
    TEST("optimalIterations(6) = 6",  Grover::optimalIterations(6) == 6);
    TEST("optimalIterations(8) = 12", Grover::optimalIterations(8) == 12);

    int target = 9;
    Grover g2(4, [target](size_t s){ return (int)s == target; });
    g2.run();
    auto probs = g2.qreg().probabilities();
    double targetProb = probs[target];
    bool isMax = true;
    for (int i = 0; i < 16; ++i)
        if (i != target && probs[i] > targetProb) isMax = false;
    TEST("Target has highest probability after search", isMax);
}

// ── QFT ──────────────────────────────────────────────────────

void testQFT() {
    std::cout << "\n-- QFT --\n";

    QubitRegister q(3);
    QFT::applyFull(q, false);
    auto probs = q.probabilities();
    double expected = 1.0 / 8.0;
    bool isUniform = true;
    for (double p : probs)
        if (std::abs(p - expected) > 1e-9) isUniform = false;
    TEST("QFT|000> = uniform superposition", isUniform);

    QubitRegister q2(4);
    Gates::applyX(q2, 0);
    QFT::applyFull(q2, false);
    QFT::applyFull(q2, true);
    TEST("QFT then IQFT: |1000> recovered", std::abs(q2[1]) > 0.99);

    QubitRegister q3(4);
    Gates::applyX(q3, 1);
    Gates::applyX(q3, 2);
    QFT::applyFull(q3, false);
    QFT::applyFull(q3, true);
    TEST("QFT then IQFT: |0110> recovered", std::abs(q3[6]) > 0.99);

    QubitRegister q4(3);
    Gates::applyH(q4, 0);
    Gates::applyH(q4, 1);
    QFT::applyFull(q4, false);
    TEST("QFT preserves normalization", q4.isNormalized());
}

// ── Teleportation ────────────────────────────────────────────

void testTeleportation() {
    std::cout << "\n-- Quantum Teleportation --\n";

    struct Case { const char* name; Complex alpha, beta; };
    std::vector<Case> cases = {
        {"|0>", {1,0},           {0,0}},
        {"|1>", {0,0},           {1,0}},
        {"|+>", {M_SQRT1_2,0},   {M_SQRT1_2,0}},
        {"|->", {M_SQRT1_2,0},   {-M_SQRT1_2,0}},
        {"|i>", {M_SQRT1_2,0},   {0,M_SQRT1_2}},
        {"|T>", {M_SQRT1_2,0},   {0.5,0.5}},
    };

    for (auto& c : cases) {
        Teleportation t(c.alpha, c.beta);
        t.run();
        TEST((std::string("Teleport ") + c.name + ": fidelity=1").c_str(),
             t.fidelity() > 0.999);
        TEST((std::string("Teleport ") + c.name + ": success").c_str(),
             t.success());
    }

    auto t2 = Teleportation::fromAngle(M_PI / 7);
    t2.run();
    TEST("Teleport fromAngle(pi/7): success", t2.success());
}

// ── BB84 — fixed to actual header API ────────────────────────
// Actual BB84 API you pasted:
//   BB84(int keyLength, unsigned seed = ...)
//   setEavesdropping(bool, double)
//   run()
//   isKeySecure(), getSiftedKey(), getErrorRate()

void testBB84() {
    std::cout << "\n-- BB84 QKD --\n";

    BB84 clean(300);
    clean.run();
    TEST("BB84 no Eve: sifted key non-empty", !clean.getSiftedKey().empty());
    TEST("BB84 no Eve: error rate < 5%", clean.getErrorRate() < 0.05);
    TEST("BB84 no Eve: key is secure", clean.isKeySecure());

    BB84 eve(600);
    eve.setEavesdropping(true, 1.0);
    eve.run();
    TEST("BB84 with Eve: error rate > 10%", eve.getErrorRate() > 0.10);
    TEST("BB84 with Eve: key NOT secure", !eve.isKeySecure());

    BB84 evePartial(800);
    evePartial.setEavesdropping(true, 0.5);
    evePartial.run();
    TEST("BB84 partial Eve: some errors", evePartial.getErrorRate() > 0.0);

    BB84 small(50);
    small.run();
    BB84 large(200);
    large.run();
    TEST("Larger BB84 gives longer sifted key",
         large.getSiftedKey().size() > small.getSiftedKey().size());
}

// ── QAOA — fixed to actual header API ────────────────────────
// Actual QAOA API you pasted:
//   QAOA(int numQubits, int p, unsigned seed = 42)
//   setMaxCutProblem(vector<pair<int,int>>)
//   setCostFunction(function<double(size_t)>)
//   setParameters(gammas, betas)
//   optimizeParameters(int)
//   run()
//   getBestSolution(), getBestCost(), getGammas(), getBetas()
//   expectationValue()

static double maxCutCostFromEdges(
    size_t bitstring,
    const std::vector<std::pair<int,int>>& edges
) {
    double cut = 0.0;
    for (auto [u, v] : edges) {
        int bu = (bitstring >> u) & 1;
        int bv = (bitstring >> v) & 1;
        if (bu != bv) cut += 1.0;
    }
    return cut;
}

void testQAOA() {
    std::cout << "\n-- QAOA --\n";

    {
        std::vector<std::pair<int,int>> edges = {
            {0, 1}, {1, 2}, {2, 3}
        };

        QAOA qaoa(4, 1);
        qaoa.setMaxCutProblem(edges);
        qaoa.setCostFunction([edges](size_t s) {
            return maxCutCostFromEdges(s, edges);
        });
        qaoa.setParameters({0.5}, {0.5});
        qaoa.run();

        TEST("QAOA path: best cost >= 0", qaoa.getBestCost() >= 0.0);
        TEST("QAOA path: best solution in range", qaoa.getBestSolution() < (1ULL << 4));
        TEST("QAOA path: gamma count == p", qaoa.getGammas().size() == 1);
        TEST("QAOA path: beta count == p", qaoa.getBetas().size() == 1);
        TEST("QAOA path: expectation finite", std::isfinite(qaoa.expectationValue()));
    }

    {
        std::vector<std::pair<int,int>> triangle = {
            {0, 1}, {1, 2}, {0, 2}
        };

        QAOA qaoa2(3, 1);
        qaoa2.setMaxCutProblem(triangle);
        qaoa2.setCostFunction([triangle](size_t s) {
            return maxCutCostFromEdges(s, triangle);
        });
        qaoa2.optimizeParameters(25);
        qaoa2.run();

        TEST("QAOA triangle: best cost >= 1", qaoa2.getBestCost() >= 1.0);
        TEST("QAOA triangle: best solution in range", qaoa2.getBestSolution() < (1ULL << 3));
        TEST("QAOA triangle: expectation finite", std::isfinite(qaoa2.expectationValue()));
    }
}

// ── ThreeQubitFlip — fixed to actual header API ──────────────
// Actual ThreeQubitFlip API you pasted:
//   ThreeQubitFlip(unsigned seed = 42)
//   encode(alpha, beta)
//   injectError(int qubit)
//   detectAndCorrect()
//   decode() -> pair<Complex, Complex>
//   success()

void testThreeQubitFlip() {
    std::cout << "\n-- Three-qubit bit-flip code --\n";

    {
        ThreeQubitFlip code;
        code.encode(Complex{1.0, 0.0}, Complex{0.0, 0.0});
        code.detectAndCorrect();
        auto [a, b] = code.decode();

        TEST("ThreeQubitFlip no error: alpha preserved", ceq(a, Complex{1.0, 0.0}));
        TEST("ThreeQubitFlip no error: beta preserved",  ceq(b, Complex{0.0, 0.0}));
        TEST("ThreeQubitFlip no error: success", code.success());
    }

    for (int q = 0; q < 3; ++q) {
        ThreeQubitFlip code;
        Complex alpha{M_SQRT1_2, 0.0};
        Complex beta {M_SQRT1_2, 0.0};

        code.encode(alpha, beta);
        code.injectError(q);
        code.detectAndCorrect();
        auto [a, b] = code.decode();

        TEST(("ThreeQubitFlip: error on qubit " + std::to_string(q) + " corrected (alpha)").c_str(),
             ceq(a, alpha, 1e-6));
        TEST(("ThreeQubitFlip: error on qubit " + std::to_string(q) + " corrected (beta)").c_str(),
             ceq(b, beta, 1e-6));
        TEST(("ThreeQubitFlip: success on qubit " + std::to_string(q)).c_str(),
             code.success());
    }
}

// ── ShorCode — fixed to actual header API ────────────────────
// Actual ShorCode API you pasted:
//   ShorCode(unsigned seed = 42)
//   encode()
//   injectError(int qubit, char errorType)
//   decode()
//   corrected()

void testShorCode() {
    std::cout << "\n-- Shor 9-qubit error correction --\n";

    {
        ShorCode sc;
        sc.encode();
        sc.decode();
        TEST("ShorCode no error: corrected", sc.corrected());
    }

    for (int q : {0, 1, 4, 8}) {
        ShorCode sc;
        sc.encode();
        sc.injectError(q, 'X');
        sc.decode();
        TEST(("ShorCode: X on qubit " + std::to_string(q) + " corrected").c_str(),
             sc.corrected());
    }

    for (int q : {0, 3, 6}) {
        ShorCode sc;
        sc.encode();
        sc.injectError(q, 'Z');
        sc.decode();
        TEST(("ShorCode: Z on qubit " + std::to_string(q) + " corrected").c_str(),
             sc.corrected());
    }

    for (int q : {2, 5, 7}) {
        ShorCode sc;
        sc.encode();
        sc.injectError(q, 'Y');
        sc.decode();
        TEST(("ShorCode: Y on qubit " + std::to_string(q) + " corrected").c_str(),
             sc.corrected());
    }
}

int main() {
    std::cout << "=== QSE Algorithm Tests ===\n";
    testGrover();
    testQFT();
    testTeleportation();
    testBB84();
    testQAOA();
    testThreeQubitFlip();
    testShorCode();

    std::cout << "\n=== Results: " << passed << " passed, " << failed << " failed ===\n";
    return (failed > 0) ? 1 : 0;
}