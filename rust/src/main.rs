// EvictVC -- Rust/arkworks backend over BLS12-381 (parallel).
//
// commit(v): iFFT(v) over a roots-of-unity domain, then MSM against the
//   monomial SRS [tau^k]_1.  O(n log n) field FFT + O(n) MSM.
// update/evict: one group op against a precomputed Lagrange-basis commitment
//   L_i = [l_i(tau)]_1.  O(1), independent of n.  (The full Lagrange SRS is a
//   one-time setup artifact; we precompute a representative L_i to time it.)
// verify: two pairings, O(1).

use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective, G2Projective};
use ark_ec::{pairing::Pairing, CurveGroup, Group, VariableBaseMSM};
use ark_ff::{One, UniformRand, Zero};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_std::rand::{rngs::StdRng, SeedableRng};
use std::time::Instant;

struct Params {
    srs1_affine: Vec<G1Affine>,
    g2: G2Projective,
    g2_tau: G2Projective,
    domain: GeneralEvaluationDomain<Fr>,
}

fn div_by_linear(p: &[Fr], z: Fr) -> Vec<Fr> {
    let d = p.len() - 1;
    let mut q = vec![Fr::zero(); d];
    if d == 0 {
        return q;
    }
    q[d - 1] = p[d];
    let mut k = d as isize - 2;
    while k >= 0 {
        let ku = k as usize;
        q[ku] = p[ku + 1] + z * q[ku + 1];
        k -= 1;
    }
    q
}

fn setup(n: usize, rng: &mut StdRng) -> Params {
    let tau = Fr::rand(rng);
    let g1 = G1Projective::generator();
    let g2 = G2Projective::generator();
    let mut powers = Vec::with_capacity(n);
    let mut p = Fr::one();
    for _ in 0..n {
        powers.push(p);
        p *= tau;
    }
    let srs1: Vec<G1Projective> = powers.iter().map(|&pk| g1 * pk).collect();
    let srs1_affine = G1Projective::normalize_batch(&srs1);
    let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
    Params { srs1_affine, g2, g2_tau: g2 * tau, domain }
}

fn lagrange_point(i: usize, pp: &Params) -> G1Projective {
    let n = pp.srs1_affine.len();
    let mut e = vec![Fr::zero(); n];
    e[i] = Fr::one();
    let coeffs = pp.domain.ifft(&e);
    G1Projective::msm(&pp.srs1_affine, &coeffs).unwrap()
}

fn commit(v: &[Fr], pp: &Params) -> G1Projective {
    let coeffs = pp.domain.ifft(v);
    G1Projective::msm(&pp.srs1_affine, &coeffs).unwrap()
}

fn update(c: G1Projective, l_i: G1Projective, delta: Fr) -> G1Projective {
    c + l_i * delta
}

fn open(v: &[Fr], i: usize, pp: &Params) -> (Fr, G1Projective) {
    let coeffs = pp.domain.ifft(v);
    let y = v[i];
    let mut num = coeffs;
    num[0] -= y;
    let z = pp.domain.element(i);
    let q = div_by_linear(&num, z);
    let proof = G1Projective::msm(&pp.srs1_affine[..q.len()], &q).unwrap();
    (y, proof)
}

fn verify(c: G1Projective, i: usize, y: Fr, proof: G1Projective, pp: &Params) -> bool {
    let g1 = G1Projective::generator();
    let z = pp.domain.element(i);
    let lhs_g1 = c - g1 * y;
    let rhs_g2 = pp.g2_tau - pp.g2 * z;
    Bls12_381::pairing(lhs_g1.into_affine(), pp.g2.into_affine())
        == Bls12_381::pairing(proof.into_affine(), rhs_g2.into_affine())
}

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
