#pragma once
#include "QubitRegister.h"
#include <string>
#include <cmath>

namespace QSE {

// Pre-built 2x2 gate matrices and multi-qubit gate applications
namespace Gates {

// ── Type alias ──────────────────────────────────────────────
using Matrix2 = std::vector<std::vector<Complex>>;

// ── Single-qubit gates ───────────────────────────────────────
const Matrix2 I  = {{{1,0},{0,0}}, {{0,0},{1,0}}};
const Matrix2 X  = {{{0,0},{1,0}}, {{1,0},{0,0}}};  // Pauli-X (NOT)
const Matrix2 Y  = {{{0,0},{0,-1}},{{0,1},{0,0}}};  // Pauli-Y
const Matrix2 Z  = {{{1,0},{0,0}}, {{0,0},{-1,0}}}; // Pauli-Z
const Matrix2 H  = {{{M_SQRT1_2,0},{M_SQRT1_2,0}},
                    {{M_SQRT1_2,0},{-M_SQRT1_2,0}}};  // Hadamard
const Matrix2 S  = {{{1,0},{0,0}}, {{0,0},{0,1}}};  // Phase (S)
const Matrix2 T  = {{{1,0},{0,0}}, {{0,0},{M_SQRT1_2,M_SQRT1_2}}}; // π/8
const Matrix2 Sdg= {{{1,0},{0,0}}, {{0,0},{0,-1}}}; // S†
const Matrix2 Tdg= {{{1,0},{0,0}}, {{0,0},{M_SQRT1_2,-M_SQRT1_2}}}; // T†
const Matrix2 SX = {{{0.5,0.5},{0.5,-0.5}},
                    {{0.5,-0.5},{0.5,0.5}}};           // √X

// ── Rotation gates (parametric) ──────────────────────────────
Matrix2 Rx(double theta);   // Rotation around X-axis
Matrix2 Ry(double theta);   // Rotation around Y-axis
Matrix2 Rz(double theta);   // Rotation around Z-axis
Matrix2 P(double phi);      // Phase gate
Matrix2 U(double theta, double phi, double lam); // General U gate

// ── Application helpers ──────────────────────────────────────
void applyH(QubitRegister& qreg, int qubit);
void applyX(QubitRegister& qreg, int qubit);
void applyY(QubitRegister& qreg, int qubit);
void applyZ(QubitRegister& qreg, int qubit);
void applyS(QubitRegister& qreg, int qubit);
void applyT(QubitRegister& qreg, int qubit);
void applySdg(QubitRegister& qreg, int qubit);
void applyTdg(QubitRegister& qreg, int qubit);
void applyRx(QubitRegister& qreg, int qubit, double theta);
void applyRy(QubitRegister& qreg, int qubit, double theta);
void applyRz(QubitRegister& qreg, int qubit, double theta);
void applyP(QubitRegister& qreg, int qubit, double phi);
void applyU(QubitRegister& qreg, int qubit, double theta, double phi, double lam);

// ── Two-qubit gates ──────────────────────────────────────────
void applyCX(QubitRegister& qreg, int control, int target);   // CNOT
void applyCY(QubitRegister& qreg, int control, int target);
void applyCZ(QubitRegister& qreg, int control, int target);
void applyCH(QubitRegister& qreg, int control, int target);
void applySWAP(QubitRegister& qreg, int qubit1, int qubit2);
void applyISWAP(QubitRegister& qreg, int qubit1, int qubit2);
void applyCRz(QubitRegister& qreg, int control, int target, double theta);
void applyCP(QubitRegister& qreg, int control, int target, double phi);

// ── Three-qubit gates ────────────────────────────────────────
void applyToffoli(QubitRegister& qreg, int c1, int c2, int target); // CCX
void applyFredkin(QubitRegister& qreg, int control, int t1, int t2); // CSWAP

// ── Matrix utilities ─────────────────────────────────────────
bool isUnitary(const Matrix2& m, double tol = 1e-10);
bool isHermitian(const Matrix2& m, double tol = 1e-10);
Matrix2 matMul(const Matrix2& a, const Matrix2& b);
Matrix2 conjugateTranspose(const Matrix2& m);
std::string matToString(const Matrix2& m);

} // namespace Gates
} // namespace QSE