//! Maintenance benchmark for the EvictVC core (the `evictvc` library crate):
//! commit / O(1) update / open / verify across n = 2^10 .. 2^max.
//!
//! Run:  cargo run --release            (n up to 2^16)
//!       cargo run --release -- 16      (set the maximum log2 n)

use ark_bls12_381::Fr;
use ark_ff::{One, UniformRand};
use ark_std::rand::{rngs::StdRng, SeedableRng};
use std::time::Instant;

use evictvc::*;

fn main() {
    let mut rng = StdRng::seed_from_u64(42);

    {
        let n = 16;
        let pp = setup(n, &mut rng);
        let v: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let c = commit(&v, &pp);
        let (y, proof) = open(&v, 3, &pp);
        assert!(verify(c, 3, y, proof, &pp), "honest opening must verify");
        assert!(!verify(c, 3, y + Fr::one(), proof, &pp), "tampered must fail");
        let l5 = lagrange_point(5, &pp);
        let delta = Fr::rand(&mut rng);
        let c2 = update(c, l5, delta);
        let mut v2 = v.clone();
        v2[5] += delta;
        assert_eq!(c2, commit(&v2, &pp), "O(1) update must equal fresh commit");
        println!("correctness: OK  (open verifies, tamper rejected, O(1) update consistent)");
        println!();
    }

    let max_log: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(16);
    println!(
        "{:<8}{:>13}{:>13}{:>13}{:>13}{:>13}",
        "n", "setup(ms)", "commit(ms)", "update(us)", "open(ms)", "verify(ms)"
    );
    for log_n in 10..=max_log {
        let n = 1usize << log_n;

        let t = Instant::now();
        let pp = setup(n, &mut rng);
        let t_setup = t.elapsed().as_secs_f64() * 1e3;

        let v: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();

        let t = Instant::now();
        let c = commit(&v, &pp);
        let t_commit = t.elapsed().as_secs_f64() * 1e3;

        let l0 = lagrange_point(0, &pp);
        let reps = 5000usize;
        let delta = Fr::rand(&mut rng);
        let t = Instant::now();
        let mut acc = c;
        for _ in 0..reps {
            acc = update(acc, l0, delta);
        }
        let _ = std::hint::black_box(acc);
        let t_update = t.elapsed().as_secs_f64() / reps as f64 * 1e6;

        let t = Instant::now();
        let (y, proof) = open(&v, 7, &pp);
        let t_open = t.elapsed().as_secs_f64() * 1e3;

        let t = Instant::now();
        let ok = verify(c, 7, y, proof, &pp);
        let t_verify = t.elapsed().as_secs_f64() * 1e3;
        assert!(ok);

        println!(
            "2^{:<6}{:>13.2}{:>13.2}{:>13.3}{:>13.2}{:>13.2}",
            log_n, t_setup, t_commit, t_update, t_open, t_verify
        );
    }
}
