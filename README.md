# EvictVC: Evictable Vector Commitments

> **Verifiable LLM inference under KV-cache eviction.**

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![Paper: ePrint](https://img.shields.io/badge/paper-ePrint-b31b1b.svg)
![Backend: ark--bls12--381](https://img.shields.io/badge/Rust-ark--bls12--381-orange.svg)
![Reference: Python](https://img.shields.io/badge/Python-reference-3776AB.svg)

**EvictVC** is a vector-commitment primitive that supports **permanent, policy-checked
deletion** of committed positions while maintaining a **tamper-evident history** of
every deletion. It closes a gap in verifiable LLM inference: production KV-cache
compressors (StreamingLLM, H2O, SnapKV) *permanently evict* the key/value states of
tokens deemed unimportant, which breaks commitment schemes that assume the committed
state stays fully present and openable.

This repository contains the paper preprint, a **Python reference implementation**
(optimized for clarity) and a **high-performance Rust backend** (`evictvc-rs`, over
`ark-bls12-381`) with reproducible benchmarks.

---

## Core contributions

- **The problem.** A formal treatment of verifiable inference under *permanent* KV-cache
  eviction — a setting existing verifiable-inference commitments cannot express.
- **The primitive.** *Evictable vector commitments*, with three game-based security
  notions: **eviction soundness**, **policy-binding**, and **history tamper-evidence**.
- **The construction.** From KZG polynomial commitments over roots of unity: **O(1)
  eviction**, a **48-byte constant-size audit proof for an arbitrary set of positions**,
  and **O(n log n)** maintenance of *all* opening proofs via the Feist–Khovratovich
  technique.
- **Security.** Tight reductions to the **n-SDH** assumption and the **collision
  resistance** of a hash function.

## Performance highlights

Single thread, Apple M1, curve BLS12-381 (`ark-bls12-381`):

| Operation | Cost | Scaling |
| :-- | :-- | :-- |
| Eviction / commitment update | **≈ 135 µs** | **constant** (n = 2¹⁰ … 2¹⁶) |
| Verification | 2 pairings, **≈ 2.6 ms** | constant |
| Audit proof size | **48 bytes** | **constant** (any number of positions) |
| Trusted setup (n = 2¹⁶) | ≈ 9.95 s | O(n log n) |
| Commit / single open (n = 2¹⁶) | ≈ 0.68 s | O(n log n) |
| All n opening proofs (n = 2¹⁶) | **≈ 7 min (FK)** vs ≈ 12.5 h (naïve) | O(n log n) vs O(n²) |

The headline result: an audit proof for *any* number of audited positions is a single
48-byte group element, versus Merkle multiproofs that grow linearly with the number of
positions.

---

## Repository layout

```
evictvc/
├── paper/                 Preprint (source + PDF) and figure
│   ├── evictvc.pdf
│   ├── evictvc.tex
│   └── evictvc_fig1.png
├── python/                Reference implementation (clarity over speed)
│   ├── kzg_vc.py          KZG vector commitment, O(1) update
│   ├── evictvc.py         Eviction, policy check, hash-chain history, subvector openings
│   └── requirements.txt
├── rust/                  evictvc-rs: high-performance backend
│   ├── src/
│   │   ├── main.rs        Maintenance benchmark (commit / update / evict / verify)
│   │   └── bin/fk.rs      Feist–Khovratovich all-proofs benchmark
│   ├── Cargo.toml
│   ├── Cargo.lock
│   └── rust-toolchain.toml
└── benchmarks/
    └── results/           Raw timing CSVs reported in the paper
```

---

## Quickstart

### Python reference implementation

```bash
cd python
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python kzg_vc.py
python evictvc.py
```

`evictvc.py` runs the eviction soundness / policy-binding / history-tamper checks and
regenerates the audit-proof-size figure (`evictvc_fig1.png`).

### Rust backend (performance)

Requires a Rust toolchain (pinned in `rust/rust-toolchain.toml`).

```bash
cd rust
cargo build --release
cargo run --release -- 16
cargo run --release --bin fk -- 16
```

The trailing integer is the maximum `log2(n)` swept by the benchmark (default 16, i.e.
n up to 2¹⁶). The first binary measures commit / update / evict / verify; `fk` measures
amortized all-proof generation via Feist–Khovratovich against the naïve baseline.

To enable multi-threaded MSM/FFT, build with the `parallel` feature:

```bash
cargo run --release --features parallel -- 16
```

---

## Reproducing the paper

- **Figure (audit proof size, constant 48 B vs linear Merkle):** `python/evictvc.py`.
- **Maintenance table (update O(1), verify constant):** `rust/`, `cargo run --release -- 16`.
- **Amortized all-proofs table (FK vs naïve):** `rust/`, `cargo run --release --bin fk -- 16`.

Raw numbers from the reference machine are in `benchmarks/results/`. Absolute timings
depend on hardware; the *asymptotic* behavior (flat update cost, constant proof size,
O(n log n) all-proofs) is hardware-independent.

---

## Security and status

This is a **research prototype** accompanying an academic paper. The implementations
are written for clarity and benchmarking. **The code is not constant-time, has not been
independently audited, and must not be used in production.** The trusted setup here is
generated locally for benchmarking; a real deployment requires a properly run
trusted-setup ceremony (or a transparent/updatable alternative).

The scheme is proven secure under the **n-SDH** assumption and the **collision
resistance** of the history hash. See the paper for the formal model, the construction,
and the reductions.

---

## Citation

If you use EvictVC in your research, please cite the paper:

```bibtex
@misc{liu2026evictvc,
  author       = {Liu, Kefan},
  title        = {Evictable Vector Commitments: Verifiable {LLM} Inference under
                  {KV}-Cache Eviction},
  year         = {2026},
  howpublished = {Cryptology {ePrint} Archive, Paper 2026/XXXX},
  note         = {\url{https://eprint.iacr.org/2026/XXXX}}
}
```

(Update the ePrint identifier once the preprint is assigned one.)

---

## License

Released under the [MIT License](LICENSE).

## Acknowledgments

<!-- Optional: grants, institutions, people, and any tooling you wish to credit. -->
