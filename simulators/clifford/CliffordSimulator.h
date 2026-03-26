#pragma once
#include "../../core/QubitRegister.h"
#include <vector>
#include <string>

namespace QSE {

// Clifford (stabilizer) simulator using the Aaronson-Gottesman tableau.
// Efficiently simulates circuits composed only of H, S, CX (and derived gates).
// Memory: O(n^2) — handles thousands of qubits where statevector simulation can't.
//
// Tableau layout: 2n rows x (2n+1) columns  [x-block | z-block | phase]
//   Rows  0..n-1  = destabilizer generators
//   Rows  n..2n-1 = stabilizer generators
class CliffordSimulator {
public:
    explicit CliffordSimulator(int numQubits);

    // ── Clifford gates (complete generating set) ─────────────
    void h(int qubit);                 // Hadamard
    void s(int qubit);                 // Phase (S)
    void sdg(int qubit);               // S-dagger
    void cx(int control, int target);  // CNOT
    void x(int qubit);                 // Pauli X
    void y(int qubit);                 // Pauli Y
    void z(int qubit);                 // Pauli Z = SS
    void cz(int control, int target);  // CZ
    void cy(int control, int target);  // CY
    void swap(int q1, int q2);         // SWAP via 3x CNOT

    // ── Measurement ──────────────────────────────────────────
    // Returns 0 or 1; collapses the tableau.
    int measure(int qubit);

    // ── Inspection ───────────────────────────────────────────
    std::string stabilizers() const;
    std::string destabilizers() const;
    std::string fullTableau() const;

    int numQubits() const { return n_; }
    void reset();

private:
    int n_;
    // tab_[i][j]:  j in [0, n)   -> x part of qubit j in generator i
    //              j in [n, 2n)  -> z part of qubit (j-n) in generator i
    //              j == 2n       -> phase bit (0 = +, 1 = -)
    std::vector<std::vector<int>> tab_;

    void rowsum(int h, int i);
    int  g(int x1, int z1, int x2, int z2) const;
};

} // namespace QSE