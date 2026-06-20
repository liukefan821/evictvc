#!/usr/bin/env python3
"""
EvictVC -- Step 1: KZG-based updatable vector commitment (correctness demo).
Backend: BLS12-381 via py_ecc (pure Python).
NOTE: setup() samples a secret tau locally for prototyping. In a real
deployment tau comes from a trusted-setup ceremony, never known to one party.
"""
import secrets
import time
from py_ecc.optimized_bls12_381 import (
    G1, G2, Z1, add, multiply, pairing, curve_order, normalize,
)

R = curve_order  # scalar field order

# ---------- polynomial arithmetic over F_R (coeffs low -> high) ----------
def poly_mul(a, b):
    res = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai % R:
            for j, bj in enumerate(b):
                res[i + j] = (res[i + j] + ai * bj) % R
    return res

def poly_eval(p, x):
    acc = 0
    for c in reversed(p):
        acc = (acc * x + c) % R
    return acc

def poly_div_by_linear(p, z):
    """Divide p(X) by (X - z). Returns (quotient, remainder)."""
    d = len(p) - 1
    if d < 1:
        return [], p[0] % R
    q = [0] * d
    q[d - 1] = p[d] % R
    for k in range(d - 2, -1, -1):
        q[k] = (p[k + 1] + z * q[k + 1]) % R
    rem = (p[0] + z * q[0]) % R
    return q, rem

def lagrange_basis_coeffs(xs):
    """Coefficient-vectors for each Lagrange basis poly l_i(X) over nodes xs."""
    M = [1]
    for x in xs:
        M = poly_mul(M, [(-x) % R, 1])           # M *= (X - x)
    basis = []
    for i in range(len(xs)):
        Mi, _ = poly_div_by_linear(M, xs[i])     # M / (X - xs[i])
        denom = poly_eval(Mi, xs[i])             # = M'(xs[i])
        inv = pow(denom, -1, R)
        basis.append([(c * inv) % R for c in Mi])
    return basis

# ---------- KZG group operations ----------
def g1_multiexp(scalars, points):
    acc = Z1
    for s, P in zip(scalars, points):
        s %= R
        if s:
            acc = add(acc, multiply(P, s))
    return acc

def setup(n):
    """Trusted setup for length-n vectors over domain {0,1,...,n-1}."""
    tau = secrets.randbelow(R - 1) + 1
    srs1, p = [], 1
    for _ in range(n):                           # [1, tau, tau^2, ...]_1
        srs1.append(multiply(G1, p))
        p = (p * tau) % R
    xs = list(range(n))
    basis = lagrange_basis_coeffs(xs)
    L = [g1_multiexp(bi, srs1) for bi in basis]  # L_i = [l_i(tau)]_1
    return {"n": n, "srs1": srs1, "g2_one": G2, "g2_tau": multiply(G2, tau),
            "L": L, "xs": xs, "basis": basis}

def commit_vector(vec, pp):
    """C = sum_i vec[i] * L_i  ( = [phi(tau)]_1, phi interpolates vec )."""
    return g1_multiexp(vec, pp["L"])

def update(C, i, delta, pp):
    """New commitment after vec[i] += delta. ONE group operation."""
    return add(C, multiply(pp["L"][i], delta % R))

def open_position(vec, i, pp):
    """Prove vec[i]. Returns (value, proof)."""
    n = pp["n"]
    phi = [0] * n
    for vi, bi in zip(vec, pp["basis"]):
        for k in range(n):
            phi[k] = (phi[k] + vi * bi[k]) % R
    z = pp["xs"][i]
    y = poly_eval(phi, z)
    num = phi[:]
    num[0] = (num[0] - y) % R
    q, rem = poly_div_by_linear(num, z)
    assert rem == 0, "phi(i) != y (interpolation bug)"
    return y, g1_multiexp(q, pp["srs1"])

def verify(C, i, y, proof, pp):
    """Check e(C - y*G1, [1]_2) == e(proof, [tau - i]_2)."""
    z = pp["xs"][i]
    lhs_g1 = add(C, multiply(G1, (-y) % R))                       # C - y*G1
    rhs_g2 = add(pp["g2_tau"], multiply(pp["g2_one"], (-z) % R))  # [tau - i]_2
    return pairing(pp["g2_one"], lhs_g1) == pairing(rhs_g2, proof)

# ---------- demo ----------
def same_point(P, Q):
    return normalize(P) == normalize(Q)

if __name__ == "__main__":
    n = 8
    print(f"[setup] domain size n = {n} (sampling SRS + Lagrange basis)...")
    t = time.perf_counter()
    pp = setup(n)
    print(f"[setup] done in {time.perf_counter() - t:.3f}s\n")

    vec = [11, 23, 47, 5, 88, 200, 17, 99]   # distinct -> non-degenerate phi
    C = commit_vector(vec, pp)
    print(f"vector            = {vec}")
    print("commitment C      = (one G1 point)\n")

    i = 3
    y, proof = open_position(vec, i, pp)
    ok = verify(C, i, y, proof, pp)
    print(f"[open]   position {i} -> value {y}")
    print(f"[verify] honest opening accepted : {ok}")

    ok_bad = verify(C, i, (y + 1) % R, proof, pp)
    print(f"[verify] tampered opening rejected: {not ok_bad}\n")

    j, delta = 5, 1000
    t = time.perf_counter()
    C2 = update(C, j, delta, pp)
    upd_ms = (time.perf_counter() - t) * 1000
    vec2 = vec[:]; vec2[j] = (vec2[j] + delta) % R
    y2, proof2 = open_position(vec2, j, pp)
    ok_upd = verify(C2, j, y2, proof2, pp)
    fresh = commit_vector(vec2, pp)
    print(f"[update] vec[{j}] += {delta} via ONE group op ({upd_ms:.3f} ms)")
    print(f"[verify] updated commitment opens correctly      : {ok_upd}")
    print(f"[check]  updated C == fresh commit of new vector  : {same_point(C2, fresh)}")
    print("\nAll core operations working. Next: add evict() + scale to a real benchmark.")
