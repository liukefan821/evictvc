#!/usr/bin/env python3
"""
EvictVC -- Step 2: eviction + the honest Figure-1 benchmark.

Adds to Step 1:
  * evict(i): remove a position from the LIVE commitment in ONE group op,
    bind the evicted record into a tamper-evident hash-chain HISTORY, and
    carry a POLICY opening proving the evicted position was eligible (its
    committed attention score was below the eviction threshold).
  * verify_eviction(): verifier checks openings + policy + live-update + history.
  * open_subvector(): one-group-element proof for an arbitrary set of positions.
  * benchmark(): batch-audit PROOF SIZE vs number of audited positions k.
        KZG = 48 bytes (one group element), constant in k.
        Merkle = k authentication paths, linear in k.
    Saves evictvc_fig1.png.

HONEST NOTE: KZG does NOT win single-update wall-clock (a group op >> a few
sha256). The real KZG advantage is proof size / aggregation, which is what the
benchmark measures. Naive proof *generation* here is O(n) (interpolation);
O(1)-amortized via precomputed proofs (Feist-Khovratovich), added later.
"""

import secrets, time, hashlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from py_ecc.optimized_bls12_381 import (
    G1, G2, Z1, Z2, add, multiply, neg, pairing, curve_order, normalize,
)

R = curve_order

# ===================== polynomial arithmetic over F_R =====================
def poly_mul(a, b):
    res = [0]*(len(a)+len(b)-1)
    for i, ai in enumerate(a):
        if ai % R:
            for j, bj in enumerate(b):
                res[i+j] = (res[i+j] + ai*bj) % R
    return res

def poly_sub(a, b):
    n = max(len(a), len(b)); res = [0]*n
    for i in range(len(a)): res[i] = a[i] % R
    for i in range(len(b)): res[i] = (res[i] - b[i]) % R
    while len(res) > 1 and res[-1] % R == 0: res.pop()
    return res

def poly_eval(p, x):
    acc = 0
    for c in reversed(p): acc = (acc*x + c) % R
    return acc

def poly_div_by_linear(p, z):
    d = len(p)-1
    if d < 1: return [], p[0] % R
    q = [0]*d
    q[d-1] = p[d] % R
    for k in range(d-2, -1, -1): q[k] = (p[k+1] + z*q[k+1]) % R
    return q, (p[0] + z*q[0]) % R

def poly_div(num, den):
    """Exact polynomial division num / den, returns quotient."""
    den = den[:]
    while len(den) > 1 and den[-1] % R == 0: den.pop()
    dd = len(den)-1
    work = [c % R for c in num]
    if len(work)-1 < dd: return [0]
    q = [0]*(len(work)-dd)
    lead_inv = pow(den[-1], -1, R)
    for i in range(len(work)-1, dd-1, -1):
        coef = (work[i]*lead_inv) % R
        q[i-dd] = coef
        for j in range(len(den)):
            work[i-dd+j] = (work[i-dd+j] - coef*den[j]) % R
    return q

def lagrange_nodes(xs, ys):
    """Interpolate through (xs[i], ys[i]) over arbitrary nodes."""
    M = [1]
    for x in xs: M = poly_mul(M, [(-x) % R, 1])
    phi = [0]*len(xs)
    for i in range(len(xs)):
        Mi, _ = poly_div_by_linear(M, xs[i])
        inv = pow(poly_eval(Mi, xs[i]), -1, R)
        sc = (ys[i]*inv) % R
        for k in range(len(Mi)): phi[k] = (phi[k] + sc*Mi[k]) % R
    return phi

def vanishing(xs):
    Z = [1]
    for x in xs: Z = poly_mul(Z, [(-x) % R, 1])
    return Z

# ===================== KZG group operations =====================
def g1_multiexp(scalars, points):
    acc = Z1
    for s, P in zip(scalars, points):
        s %= R
        if s: acc = add(acc, multiply(P, s))
    return acc

def g2_multiexp(scalars, points):
    acc = Z2
    for s, P in zip(scalars, points):
        s %= R
        if s: acc = add(acc, multiply(P, s))
    return acc

def setup(n):
    tau = secrets.randbelow(R-1)+1
    srs1, srs2, p = [], [], 1
    for _ in range(n):
        srs1.append(multiply(G1, p)); srs2.append(multiply(G2, p)); p = (p*tau) % R
    xs = list(range(n))
    basis = [None]*n
    M = [1]
    for x in xs: M = poly_mul(M, [(-x) % R, 1])
    for i in range(n):
        Mi, _ = poly_div_by_linear(M, xs[i])
        inv = pow(poly_eval(Mi, xs[i]), -1, R)
        basis[i] = [(c*inv) % R for c in Mi]
    L = [g1_multiexp(bi, srs1) for bi in basis]
    return {"n": n, "srs1": srs1, "srs2": srs2, "g2_one": G2,
            "g2_tau": multiply(G2, tau), "L": L, "xs": xs, "basis": basis}

def commit_vector(vec, pp):
    return g1_multiexp(vec, pp["L"])

def update(C, i, delta, pp):
    return add(C, multiply(pp["L"][i], delta % R))

def interpolate(vec, pp):
    n = pp["n"]; phi = [0]*n
    for vi, bi in zip(vec, pp["basis"]):
        for k in range(n): phi[k] = (phi[k] + vi*bi[k]) % R
    return phi

def open_position(vec, i, pp):
    phi = interpolate(vec, pp)
    z = pp["xs"][i]; y = poly_eval(phi, z)
    num = phi[:]; num[0] = (num[0]-y) % R
    q, rem = poly_div_by_linear(num, z)
    assert rem == 0
    return y, g1_multiexp(q, pp["srs1"])

def verify_open(C, i, y, proof, pp):
    z = pp["xs"][i]
    lhs = add(C, multiply(G1, (-y) % R))
    rhs = add(pp["g2_tau"], multiply(pp["g2_one"], (-z) % R))
    return pairing(pp["g2_one"], lhs) == pairing(rhs, proof)

def open_subvector(vec, idxs, pp):
    """ONE-group-element proof for the values at positions idxs."""
    phi = interpolate(vec, pp)
    xs_S = [pp["xs"][i] for i in idxs]; ys_S = [vec[i] for i in idxs]
    r = lagrange_nodes(xs_S, ys_S); Z = vanishing(xs_S)
    q = poly_div(poly_sub(phi, r), Z)
    return ys_S, g1_multiexp(q, pp["srs1"])

def verify_subvector(C, idxs, ys_S, proof, pp):
    xs_S = [pp["xs"][i] for i in idxs]
    r = lagrange_nodes(xs_S, ys_S); Z = vanishing(xs_S)
    rC = g1_multiexp(r, pp["srs1"]); ZC2 = g2_multiexp(Z, pp["srs2"])
    return pairing(pp["g2_one"], add(C, neg(rC))) == pairing(ZC2, proof)

def same_point(P, Q): return normalize(P) == normalize(Q)

# ===================== tamper-evident eviction history =====================
def history_init(): return hashlib.sha256(b"EvictVC-genesis").digest()
def _record(i, val, epoch): return f"{i}|{val}|{epoch}".encode()
def history_append(head, i, val, epoch):
    return hashlib.sha256(head + _record(i, val, epoch)).digest()

# ===================== eviction =====================
def evict(state_vec, score_vec, C_state, head, i, threshold, pp, epoch):
    s_val, s_proof = open_position(state_vec, i, pp)     # binds evicted value
    sc_val, sc_proof = open_position(score_vec, i, pp)   # binds policy score
    C_state_new = update(C_state, i, (-s_val) % R, pp)   # remove i (ONE group op)
    head_new = history_append(head, i, s_val, epoch)
    ev = {"value": s_val, "value_proof": s_proof, "score": sc_val, "score_proof": sc_proof}
    return C_state_new, head_new, ev

def verify_eviction(C_state_old, C_score, head_old, i, ev, threshold,
                    C_state_new, head_new, pp, epoch):
    v, vp, sc, scp = ev["value"], ev["value_proof"], ev["score"], ev["score_proof"]
    if not verify_open(C_state_old, i, v, vp, pp): return False, "bad value opening"
    if not verify_open(C_score, i, sc, scp, pp):   return False, "bad score opening"
    if not (sc < threshold):
        return False, f"policy violation: score {sc} >= threshold {threshold}"
    if not same_point(C_state_new, update(C_state_old, i, (-v) % R, pp)):
        return False, "live commitment not correctly updated"
    if head_new != history_append(head_old, i, v, epoch):
        return False, "history chain inconsistent"
    return True, "ok"

# ===================== Merkle baseline =====================
def _h(b): return hashlib.sha256(b).digest()
def _leaf(v): return _h(str(v).encode())

def merkle_build(values):
    level = [_leaf(v) for v in values]; levels = [level]
    while len(level) > 1:
        nxt = []
        for k in range(0, len(level), 2):
            right = level[k+1] if k+1 < len(level) else level[k]
            nxt.append(_h(level[k] + right))
        levels.append(nxt); level = nxt
    return levels

def merkle_multiproof_nodes(width, S):
    """Number of authentication-path nodes to prove leaves in set S."""
    frontier = set(S); nodes = 0
    while width > 1:
        nxt = set()
        for idx in frontier:
            if (idx ^ 1) not in frontier and (idx ^ 1) < width: nodes += 1
            nxt.add(idx // 2)
        frontier = nxt; width = (width + 1) // 2
    return nodes

# ===================== benchmark (Figure 1) =====================
def benchmark():
    n = 2**14                                   # fixed KV-cache size
    ks = [2**j for j in range(0, 15)]           # 1 ... 16384 audited positions
    kzg_bytes = [48]*len(ks)                    # one compressed G1, constant in k
    merkle_bytes = []
    for k in ks:
        trials = 5; tot = 0
        for _ in range(trials):
            S = set(secrets.randbelow(n) for _ in range(k))
            tot += merkle_multiproof_nodes(n, S)
        merkle_bytes.append(tot/trials * 32)    # 32 bytes per node
        print(f"k={k:6d}  KZG=48 B   Merkle={merkle_bytes[-1]:12.0f} B")

    plt.figure(figsize=(7,5))
    plt.loglog(ks, kzg_bytes, "o-", lw=2, label="EvictVC (KZG)  -- one group element, 48 B")
    plt.loglog(ks, merkle_bytes, "s-", lw=2, label="Merkle multiproof  -- k paths, O(k log n)")
    plt.xlabel("number of audited KV-cache positions  k")
    plt.ylabel("audit proof size  (bytes)")
    plt.title(f"EvictVC: batch-audit proof size  (n = {n} cache positions)")
    plt.legend(); plt.grid(True, which="both", ls=":"); plt.tight_layout()
    plt.savefig("evictvc_fig1.png", dpi=150)
    print("\nSaved figure -> evictvc_fig1.png")

# ===================== demo =====================
if __name__ == "__main__":
    print("="*64)
    print("PART 1  eviction: policy-bound + tamper-evident history")
    print("="*64)
    n = 8
    pp = setup(n)
    state_vec = [11, 23, 47, 5, 88, 200, 17, 99]
    score_vec = [90, 12, 75, 8, 60, 95, 40, 7]      # committed attention scores
    threshold, epoch = 20, 0                          # evict if score < 20
    C_state = commit_vector(state_vec, pp)
    C_score = commit_vector(score_vec, pp)
    head = history_init()

    i = 3                                             # score 8 < 20 -> eligible
    C2, head2, ev = evict(state_vec, score_vec, C_state, head, i, threshold, pp, epoch)
    t = time.perf_counter(); _ = update(C_state, i, (-ev["value"]) % R, pp)
    evict_us = (time.perf_counter()-t)*1e6
    ok, msg = verify_eviction(C_state, C_score, head, i, ev, threshold, C2, head2, pp, epoch)
    print(f"[evict {i}]  score={ev['score']} (<{threshold})   verify: {ok}  ({msg})")
    print(f"            live-commitment update = ONE group op, {evict_us:.0f} us (O(1))")

    j = 5                                             # score 95 -> ineligible
    Cb, hb, evb = evict(state_vec, score_vec, C_state, head, j, threshold, pp, epoch)
    okb, msgb = verify_eviction(C_state, C_score, head, j, evb, threshold, Cb, hb, pp, epoch)
    print(f"[evict {j}]  score={evb['score']} (>= {threshold})  verify: {okb}  ({msgb})")

    ev_t = dict(ev); ev_t["value"] = (ev["value"]+1) % R
    okt, msgt = verify_eviction(C_state, C_score, head, i, ev_t, threshold, C2, head2, pp, epoch)
    print(f"[tamper evicted value]              verify: {okt}  ({msgt})")

    print("\n" + "="*64)
    print("PART 2a  KZG subvector proof is ONE element regardless of k")
    print("="*64)
    pp2 = setup(64)
    vec = [secrets.randbelow(10**6) for _ in range(64)]
    C = commit_vector(vec, pp2)
    for k in (1, 2, 4, 8):
        idxs = list(range(k))
        ys, proof = open_subvector(vec, idxs, pp2)
        ok = verify_subvector(C, idxs, ys, proof, pp2)
        print(f"  k={k:2d}  proof = 1 G1 element (48 B compressed)   verify: {ok}")

    print("\n" + "="*64)
    print("PART 2b  benchmark (Figure 1)")
    print("="*64)
    benchmark()
