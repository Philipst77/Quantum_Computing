#include <iostream>
#include "../core/QubitRegister.h"
#include "../core/QuantumGate.h"

using namespace QSE;

int main() {
    std::cout << "=== Bell State Demo ===\n\n";

    // Gate matrices
    static const std::vector<std::vector<Complex>> H_GATE = {
        {{ 1.0/std::sqrt(2.0), 0}, { 1.0/std::sqrt(2.0), 0}},
        {{ 1.0/std::sqrt(2.0), 0}, {-1.0/std::sqrt(2.0), 0}}
    };
    static const std::vector<std::vector<Complex>> X_GATE = {
        {{0,0},{1,0}},
        {{1,0},{0,0}}
    };

    // ── Bell state |Φ+> = (|00> + |11>) / sqrt(2) ──
    {
        std::cout << "--- |Phi+> = (|00> + |11>)/sqrt(2) ---\n";
        QubitRegister reg(2);
        reg.reset();
        reg.applyGate(0, H_GATE);
        reg.applyControlledGate(0, 1, X_GATE);
        reg.print();
    }

    // ── Bell state |Φ-> = (|00> - |11>) / sqrt(2) ──
    {
        std::cout << "\n--- |Phi-> = (|00> - |11>)/sqrt(2) ---\n";
        QubitRegister reg(2);
        reg.reset();
        reg.applyGate(0, H_GATE);
        reg.applyControlledGate(0, 1, X_GATE);
        reg.applyGate(0, X_GATE); // Z gate effect via X on qubit 0
        reg.print();
    }

    // ── Bell state |Ψ+> = (|01> + |10>) / sqrt(2) ──
    {
        std::cout << "\n--- |Psi+> = (|01> + |10>)/sqrt(2) ---\n";
        QubitRegister reg(2);
        reg.reset();
        reg.applyGate(1, X_GATE); // flip qubit 1 first
        reg.applyGate(0, H_GATE);
        reg.applyControlledGate(0, 1, X_GATE);
        reg.print();
    }

    // ── Bell state |Ψ-> = (|01> - |10>) / sqrt(2) ──
    {
        std::cout << "\n--- |Psi-> = (|01> - |10>)/sqrt(2) ---\n";
        QubitRegister reg(2);
        reg.reset();
        reg.applyGate(1, X_GATE);
        reg.applyGate(0, H_GATE);
        reg.applyControlledGate(0, 1, X_GATE);
        reg.applyGate(0, X_GATE);
        reg.print();
    }

    // ── Sample measurements from |Φ+> ──
    std::cout << "\n--- Sampling 10 measurements from |Phi+> ---\n";
    {
        std::mt19937 rng(42);
        for (int i = 0; i < 10; ++i) {
            QubitRegister reg(2);
            reg.reset();
            reg.applyGate(0, H_GATE);
            reg.applyControlledGate(0, 1, X_GATE);
            auto bits = reg.measureAll();
            std::cout << "  Shot " << i+1 << ": |" << bits[0] << bits[1] << ">\n";
        }
    }

    return 0;
}