#include <iostream>
#include <random>
#include <cmath>
#include "../algorithms/teleportation/Teleportation.h"

using namespace QSE;

int main() {
    std::cout << "=== Quantum Teleportation Demo ===\n\n";

    // ── Example 1: teleport |+> = (|0>+|1>)/sqrt(2) ──
    {
        double inv_sqrt2 = 1.0 / std::sqrt(2.0);
        Teleportation tele(Complex(inv_sqrt2, 0), Complex(inv_sqrt2, 0));
        std::cout << "Teleporting |+> = (|0>+|1>)/sqrt(2)\n";
        tele.run();
        std::cout << tele.resultSummary() << "\n";
        std::cout << "Fidelity: " << tele.fidelity() << "\n";
    }

    // ── Example 2: teleport |1> ──
    std::cout << "\n--- Teleporting |1> ---\n";
    {
        Teleportation tele(Complex(0, 0), Complex(1, 0));
        tele.run();
        std::cout << tele.resultSummary() << "\n";
        std::cout << "Fidelity: " << tele.fidelity() << "\n";
    }

    // ── Example 3: teleport arbitrary state ──
    std::cout << "\n--- Teleporting 0.6|0> + 0.8|1> ---\n";
    {
        Teleportation tele(Complex(0.6, 0), Complex(0.8, 0));
        tele.run();
        std::cout << tele.resultSummary() << "\n";
        std::cout << "Fidelity: " << tele.fidelity() << "\n";
    }

    // ── Run 100 random teleportations and report success rate ──
    std::cout << "\n--- 100 Random Teleportations ---\n";
    std::mt19937 rng(1234);
    std::uniform_real_distribution<double> d(-1.0, 1.0);
    int successes = 0;
    for (int i = 0; i < 100; ++i) {
        double a = d(rng), b = d(rng);
        double norm = std::sqrt(a*a + b*b);
        if (norm < 1e-6) { a = 1; b = 0; norm = 1; }
        Teleportation tele(Complex(a/norm, 0), Complex(b/norm, 0));
        tele.run();
        if (tele.success()) ++successes;
    }
    std::cout << "Success rate: " << successes << "/100\n";

    return 0;
}