#include "QuantumGate.h"
#include <sstream>
#include <cmath>
#include <stdexcept>

namespace QSE {
namespace Gates {

// ── Parametric gate constructors ─────────────────────────────

Matrix2 Rx(double theta) {
    Complex c = std::cos(theta / 2.0);
    Complex s = Complex(0, -std::sin(theta / 2.0));
    return {
        {c, s},
        {s, c}
    };
}

Matrix2 Ry(double theta) {
    double c = std::cos(theta / 2.0);
    double s = std::sin(theta / 2.0);
    return {
        {{c, 0}, {-s, 0}},
        {{s, 0}, {c, 0}}
    };
}

Matrix2 Rz(double theta) {
    Complex e0 = std::exp(Complex(0, -theta/2));
    Complex e1 = std::exp(Complex(0,  theta/2));
    return {
        {e0, {0,0}},
        {{0,0}, e1}
    };
}

Matrix2 P(double phi) {
    return {
        {{1,0}, {0,0}},
        {{0,0}, {std::cos(phi), std::sin(phi)}}
    };
}

Matrix2 U(double theta, double phi, double lam) {
    double c = std::cos(theta / 2.0);
    double s = std::sin(theta / 2.0);
    return {
        {{c, 0}, {-std::cos(lam) * s, -std::sin(lam) * s}},
        {{std::cos(phi) * s, std::sin(phi) * s},
         {std::cos(phi + lam) * c, std::sin(phi + lam) * c}}
    };
}

// ── Single-qubit helpers ─────────────────────────────────────

void applyH(QubitRegister& qreg, int qubit)   { qreg.applyGate(qubit, H); }
void applyX(QubitRegister& qreg, int qubit)   { qreg.applyGate(qubit, X); }
void applyY(QubitRegister& qreg, int qubit)   { qreg.applyGate(qubit, Y); }
void applyZ(QubitRegister& qreg, int qubit)   { qreg.applyGate(qubit, Z); }
void applyS(QubitRegister& qreg, int qubit)   { qreg.applyGate(qubit, S); }
void applyT(QubitRegister& qreg, int qubit)   { qreg.applyGate(qubit, T); }
void applySdg(QubitRegister& qreg, int qubit) { qreg.applyGate(qubit, Sdg); }
void applyTdg(QubitRegister& qreg, int qubit) { qreg.applyGate(qubit, Tdg); }

void applyRx(QubitRegister& qreg, int qubit, double theta) { qreg.applyGate(qubit, Rx(theta)); }
void applyRy(QubitRegister& qreg, int qubit, double theta) { qreg.applyGate(qubit, Ry(theta)); }
void applyRz(QubitRegister& qreg, int qubit, double theta) { qreg.applyGate(qubit, Rz(theta)); }
void applyP(QubitRegister& qreg, int qubit, double phi)    { qreg.applyGate(qubit, P(phi)); }

void applyU(QubitRegister& qreg, int qubit,
            double theta, double phi, double lam) {
    qreg.applyGate(qubit, U(theta, phi, lam));
}

// ── Two-qubit gates ──────────────────────────────────────────

void applyCX(QubitRegister& qreg, int control, int target) {
    qreg.applyControlledGate(control, target, X);
}

void applyCY(QubitRegister& qreg, int control, int target) {
    qreg.applyControlledGate(control, target, Y);
}

void applyCZ(QubitRegister& qreg, int control, int target) {
    qreg.applyControlledGate(control, target, Z);
}

void applyCH(QubitRegister& qreg, int control, int target) {
    qreg.applyControlledGate(control, target, H);
}

void applySWAP(QubitRegister& qreg, int q1, int q2) {
    applyCX(qreg, q1, q2);
    applyCX(qreg, q2, q1);
    applyCX(qreg, q1, q2);
}

void applyCRz(QubitRegister& qreg, int control, int target, double theta) {
    qreg.applyControlledGate(control, target, Rz(theta));
}

void applyCP(QubitRegister& qreg, int control, int target, double phi) {
    qreg.applyControlledGate(control, target, P(phi));
}

// ── Three-qubit gates ────────────────────────────────────────

void applyToffoli(QubitRegister& qreg, int c1, int c2, int target) {
    qreg.applyToffoliGate(c1, c2, target);
}

void applyFredkin(QubitRegister& qreg, int control, int t1, int t2) {
    applyCX(qreg, t2, t1);
    applyToffoli(qreg, control, t1, t2);
    applyCX(qreg, t2, t1);
}

} // namespace Gates
} // namespace QSE