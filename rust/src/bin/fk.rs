// EvictVC -- Feist-Khovratovich (FK) amortized proofs.
//
// Computes ALL n opening proofs (one per domain position) in O(n log n) total,
// instead of O(n) per proof = O(n^2) naive. Identity used:
//     proof(z) = sum_k h_k z^k,   h_k = sum_j c_{k+j+1} s_j
// so the proofs are the DFT of a group-coefficient vector h, and h is a
// correlation of the polynomial coefficients with the SRS (computed by FFT).
//
// Run:  cargo run --release --bin fk            (default to 2^12)
//       cargo run --release --bin fk -- 16      (up to 2^16)

use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective, G2Projective};
use ark_ec::{pairing::Pairing, AffineRepr, CurveGroup, Group, VariableBaseMSM};
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

fn commit(v: &[Fr], pp: &Params) -> G1Projective {
    let coeffs = pp.domain.ifft(v);
    G1Projective::msm(&pp.srs1_affine, &coeffs).unwrap()
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

// All n proofs in O(n log n): proofs = DFT(h), h = correlation(coeffs, SRS).
fn all_proofs_fk(coeffs: &[Fr], pp: &Params) -> Vec<G1Projective> {
    let n = coeffs.len();
    let d = n - 1;
    let size2 = 2 * n;
    let domain2 = GeneralEvaluationDomain::<Fr>::new(size2).unwrap();
    let s2 = domain2.size();

    // reversed coefficient tail (field): a = [c_d, c_{d-1}, ..., c_1, 0, ...]
    let mut a = vec![Fr::zero(); s2];
    for j in 0..d {
        a[j] = coeffs[d - j];
    }
    // SRS (group): b = [s_0, s_1, ..., s_{d-1}, 0, ...]
    let mut b = vec![G1Projective::zero(); s2];
    for j in 0..d {
        b[j] = pp.srs1_affine[j].into_group();
    }

    let a_freq = domain2.fft(&a); // field FFT
    let b_freq = domain2.fft(&b); // group FFT (depends only on SRS -> precomputable)
    let prod: Vec<G1Projective> = a_freq
        .iter()
        .zip(b_freq.iter())
        .map(|(&af, &bf)| bf * af)
        .collect();
    let conv = domain2.ifft(&prod); // group iFFT -> linear convolution

    // h_m = conv[d-1-m]
    let mut h = vec![G1Projective::zero(); n];
    for mm in 0..d {
        h[mm] = conv[d - 1 - mm];
    }
    pp.domain.fft(&h) // group DFT over the domain -> all proofs
}

fn main() {
    let mut rng = StdRng::seed_from_u64(7);

    // ---- correctness: FK proofs match naive open AND all verify ----
    {
        let n = 16;
        let pp = setup(n, &mut rng);
        let v: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let c = commit(&v, &pp);
        let coeffs = pp.domain.ifft(&v);
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

    // ---- benchmark: FK all-proofs vs naive single open ----
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
        let coeffs = pp.domain.ifft(&v);

        let t = Instant::now();
        let fk = all_proofs_fk(&coeffs, &pp);
        let t_fk = t.elapsed().as_secs_f64() * 1e3;
        let per_proof_us = t_fk * 1e3 / n as f64;

        // spot-check one proof
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
