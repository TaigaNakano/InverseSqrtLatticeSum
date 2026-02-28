# InverseSqrtLatticeSum

## 1. Abstract

InverseSqrtLatticeSum evaluates the regularized 1D lattice sum

```math
S(a,b)=\sum_{n=n_0}^{\infty}\left(\frac{1}{\sqrt{n^2+a}}-\frac{1}{\sqrt{n^2+b}}\right),\qquad a,b>0.
```

with $`n_0\in\{0,1\}`$.

Although the difference is convergent, direct numerical summation can be inefficient at high accuracy:
each component behaves like $`1/n`$, and the final result relies on substantial cancellation.

This repository provides a small, embeddable C++ implementation with two methods:
a transparent direct summation baseline and a faster transform-based solver.
The goal is practical reliability: deterministic stopping criteria, explicit status codes,
and stable floating-point handling.

## 2. How to use

Build and run:

```bash
g++ -O3 main.cpp inverse_sqrt_lattice_sum.cpp -o isls
./isls
```

You can include `inverse_sqrt_lattice_sum.hpp` and call `InverseSqrtLatticeSum::get_value(...)`.

- `method=0` (`ISLS_FAST`): Poisson/transform accelerated evaluation.
- `method=1` (`ISLS_CONVENTIONAL`): direct Kahan-compensated partial sum.
- Inputs must satisfy $`a>0`$, $`b>0`$.
- Accuracy/termination are controlled by `atol`, `rtol`, and `max_evaluation`.

## 3. Theory & numerical strategy

### 3.1 Why acceleration is needed

Even though the difference cancels the leading $`O(1/n)`$ behavior and yields an $`O(n^{-3})`$ tail,
plain partial sums can still be slow and numerically delicate:

- Regime dependence: small $`a`$ or $`b`$, or strong imbalance $`a \ll b`$, can require many terms.
- Cancellation sensitivity: the final value may be much smaller than intermediate partial sums.
- Stopping criteria: absolute-only tolerances can fail when the scale of $`S(a,b)`$ varies widely.

The conventional path therefore uses compensated summation to reduce roundoff accumulation.

### 3.2 Poisson summation viewpoint (key idea)

A central observation is that $`S(a,b)`$ is a difference of two lattice sums of a smooth kernel.
Poisson summation relates sums over integers to sums over Fourier modes:

```math
\sum_{n\in\mathbb{Z}} f(n)=\sum_{k\in\mathbb{Z}} \widehat{f}(k).
```

So if the Fourier transform $`\widehat{f}`$ decays rapidly, the transformed series converges quickly.
Intuitively: more smoothness in the original function means faster decay in Fourier space, which
is precisely what makes the transformed representation efficient.
(For Schwartz-class functions this is the cleanest; for other kernels one uses regularization /
rearrangement to make the Poisson approach numerically well behaved.)

In this project, the “FAST” method uses a Poisson-summation-based reformulation that splits
$`S(a,b)`$ into (i) a closed-form main part and (ii) a rapidly decaying remainder.

### 3.3 The accelerated representation used here

The fast path uses the Poisson-summation reformulation (as implemented in this repo):

```math
S(a,b)=\frac{1}{2}\log\frac{b}{a}\;\pm\;\frac{1}{2}\left(\frac{1}{\sqrt{a}}-\frac{1}{\sqrt{b}}\right)+2\int_0^T\left[\frac{1}{e^{\alpha\cosh t}-1}-\frac{1}{e^{\beta\cosh t}-1}\right]dt,
```

with $`\alpha=2\pi\sqrt a`$, $`\beta=2\pi\sqrt b`$, and the sign determined by the index origin
($`n_0=0`$ vs $`n_0=1`$). The integrand decays exponentially in $`t`$, enabling efficient trapezoidal refinement.

Implementation notes:

- For large arguments we switch to numerically stable forms (e.g. using $`e^{-x}`$-based evaluation
  to avoid loss in $`1/(e^x-1)`$).
- Convergence is checked using a mixed criterion:

```math
|\Delta|\le \mathrm{rtol}\,|S_{\mathrm{old}}|+\mathrm{atol}.
```

### 3.4 Background & where this sum shows up (motivation)

This regularized inverse-sqrt lattice sum is a compact representative of a broader class of
(i) image-series / periodic Green-function constructions and (ii) regularized lattice sums
where naive summation is slow or unstable, motivating transform-based acceleration.

Fast evaluation of lattice sums is a standard theme in periodic scattering and related problems;
see Denlinger et al. for an open-access example of fast summation methods and regime-sensitive
behavior in lattice sums. [Denlinger]

Regularized lattice sums (zeta/theta/Epstein-type objects) provide a mathematical umbrella for
many “sum over lattice points with regularization” constructions. [Buchheit]

As a concrete applied example where infinite image arrays lead to correction factors (and where
careful evaluation of series matters), four-point probe sheet resistivity formulas are a classic
reference. [Smits]

This repository focuses on the numerical evaluation of $`S(a,b)`$ itself.

## 4. References (open access)

- [Poisson summation] Poisson summation formula (overview and standard assumptions).  
  https://en.wikipedia.org/wiki/Poisson_summation_formula

- [Fourier decay notes] Lecture notes discussing Fourier transforms and Poisson summation and the
  smoothness/decay heuristic (Fourier transform of smooth functions decays rapidly).  
  https://math.mit.edu/classes/18.785/2017fa/LectureNotes17.pdf

- [Denlinger] R. Denlinger, Z. Gimbutas, L. Greengard, V. Rokhlin,  
  A Fast Summation Method for Oscillatory Lattice Sums, arXiv:1609.08523.  
  https://arxiv.org/abs/1609.08523

- [Buchheit] A. A. Buchheit, C. Busse, J. Gutendorf,  
  Computation and properties of the Epstein zeta function ..., arXiv:2412.16317.  
  https://arxiv.org/abs/2412.16317

- [Smits] F. M. Smits,  
  Measurement of Sheet Resistivities with the Four-Point Probe (open PDF).  
  https://www.four-point-probes.com/smits.pdf
