//! Feist-Khovratovich amortized all-proofs benchmark (the `evictvc` library
//! crate): all n opening proofs in O(n log n) total vs the O(n) naive single
//! open, with a correctness check that FK proofs match naive openings and verify.
//!
//! Run:  cargo run --release --bin fk            (n up to 2^12)
//!       cargo run --release --bin fk -- 16      (n up to 2^16)

use ark_bls12_381::Fr;
use ark_ff::UniformRand;
use ark_std::rand::{rngs::StdRng, SeedableRng};
use std::time::Instant;

use evictvc::*;

fn main() {
    let mut rng = StdRng::seed_from_u64(7);

    {
        let n = 16;
        let pp = setup(n, &mut rng);
        let v: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let c = commit(&v, &pp);
        let coeffs = pp.ifft(&v);
        let fk = all_proofs_fk(&coeffs, &pp);
        let mut all_match = true;
        let mut all_verify = true;
        for i in 0..n {
            let (_, naive) = open(&v, i, &pp);
            if fk[i] != naive {
                all_match = false;
            }
            if !verify(c, i, v[i], fk[i], &pp) {
                all_verify = false;
            }
        }
        println!(
            "FK correctness (n={}): match_naive={}  all_verify={}",
            n, all_match, all_verify
        );
        assert!(all_match && all_verify);
        println!();
    }

    let max_log: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(12);
    println!(
        "{:<8}{:>15}{:>18}{:>16}{:>22}",
        "n", "FK all(ms)", "FK per-proof(us)", "naive open(ms)", "naive all est(s)"
    );
    for log_n in 8..=max_log {
        let n = 1usize << log_n;
        let pp = setup(n, &mut rng);
        let v: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let coeffs = pp.ifft(&v);

        let t = Instant::now();
        let fk = all_proofs_fk(&coeffs, &pp);
        let t_fk = t.elapsed().as_secs_f64() * 1e3;
        let per_proof_us = t_fk * 1e3 / n as f64;

        let c = commit(&v, &pp);
        assert!(verify(c, 7, v[7], fk[7], &pp));

        let t = Instant::now();
        let _ = open(&v, 7, &pp);
        let t_naive = t.elapsed().as_secs_f64() * 1e3;
        let naive_all_s = t_naive * n as f64 / 1e3;

        println!(
            "2^{:<6}{:>15.2}{:>18.3}{:>16.2}{:>22.1}",
            log_n, t_fk, per_proof_us, t_naive, naive_all_s
        );
    }
}
