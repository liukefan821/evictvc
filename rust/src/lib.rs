//! EvictVC core library: a KZG vector commitment over a roots-of-unity domain,
//! with O(1) updates/evictions via Lagrange-basis points and Feist-Khovratovich
//! amortized all-proofs. The `evictvc` and `fk` binaries are thin benchmark
//! harnesses built on top of this crate.

use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective, G2Projective};
use ark_ec::{pairing::Pairing, AffineRepr, CurveGroup, Group, VariableBaseMSM};
use ark_ff::{One, UniformRand, Zero};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_std::rand::rngs::StdRng;

/// Public parameters: the monomial SRS [tau^k]_1 in G1, the G2 generator and
/// [tau]_2, and the evaluation domain.
pub struct Params {
    pub srs1_affine: Vec<G1Affine>,
    pub g2: G2Projective,
    pub g2_tau: G2Projective,
    pub domain: GeneralEvaluationDomain<Fr>,
}

impl Params {
    /// iFFT a value vector to coefficient form over the domain.
    pub fn ifft(&self, evals: &[Fr]) -> Vec<Fr> {
        self.domain.ifft(evals)
    }
}

/// Synthetic division of a polynomial by (X - z); returns the quotient.
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

/// Trusted setup: sample tau, build the monomial SRS and the domain.
pub fn setup(n: usize, rng: &mut StdRng) -> Params {
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

/// The Lagrange-basis commitment L_i = [l_i(tau)]_1 for position i.
pub fn lagrange_point(i: usize, pp: &Params) -> G1Projective {
    let n = pp.srs1_affine.len();
    let mut e = vec![Fr::zero(); n];
    e[i] = Fr::one();
    let coeffs = pp.domain.ifft(&e);
    G1Projective::msm(&pp.srs1_affine, &coeffs).unwrap()
}

/// Commit to a value vector: O(n log n) field iFFT followed by an O(n) MSM.
pub fn commit(v: &[Fr], pp: &Params) -> G1Projective {
    let coeffs = pp.domain.ifft(v);
    G1Projective::msm(&pp.srs1_affine, &coeffs).unwrap()
}

/// O(1) update / eviction: c + L_i * delta, independent of n.
pub fn update(c: G1Projective, l_i: G1Projective, delta: Fr) -> G1Projective {
    c + l_i * delta
}

/// Open position i: returns (value, KZG proof).
pub fn open(v: &[Fr], i: usize, pp: &Params) -> (Fr, G1Projective) {
    let coeffs = pp.domain.ifft(v);
    let y = v[i];
    let mut num = coeffs;
    num[0] -= y;
    let z = pp.domain.element(i);
    let q = div_by_linear(&num, z);
    let proof = G1Projective::msm(&pp.srs1_affine[..q.len()], &q).unwrap();
    (y, proof)
}

/// Verify an opening with two pairings. O(1).
pub fn verify(c: G1Projective, i: usize, y: Fr, proof: G1Projective, pp: &Params) -> bool {
    let g1 = G1Projective::generator();
    let z = pp.domain.element(i);
    let lhs_g1 = c - g1 * y;
    let rhs_g2 = pp.g2_tau - pp.g2 * z;
    Bls12_381::pairing(lhs_g1.into_affine(), pp.g2.into_affine())
        == Bls12_381::pairing(proof.into_affine(), rhs_g2.into_affine())
}

/// All n opening proofs in O(n log n) total via Feist-Khovratovich.
/// proofs = DFT(h), where h is a correlation of the polynomial coefficients
/// with the SRS (computed by a group FFT/iFFT).
pub fn all_proofs_fk(coeffs: &[Fr], pp: &Params) -> Vec<G1Projective> {
    let n = coeffs.len();
    let d = n - 1;
    let size2 = 2 * n;
    let domain2 = GeneralEvaluationDomain::<Fr>::new(size2).unwrap();
    let s2 = domain2.size();

    let mut a = vec![Fr::zero(); s2];
    for j in 0..d {
        a[j] = coeffs[d - j];
    }
    let mut b = vec![G1Projective::zero(); s2];
    for j in 0..d {
        b[j] = pp.srs1_affine[j].into_group();
    }

    let a_freq = domain2.fft(&a);
    let b_freq = domain2.fft(&b);
    let prod: Vec<G1Projective> = a_freq
        .iter()
        .zip(b_freq.iter())
        .map(|(&af, &bf)| bf * af)
        .collect();
    let conv = domain2.ifft(&prod);

    let mut h = vec![G1Projective::zero(); n];
    for mm in 0..d {
        h[mm] = conv[d - 1 - mm];
    }
    pp.domain.fft(&h)
}
