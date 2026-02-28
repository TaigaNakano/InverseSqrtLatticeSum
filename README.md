# InverseSqrtLatticeSum

## 1. Abstract, Motivation
InverseSqrtLatticeSum evaluates
\[
S(a,b)=\sum_{n=n_0}^{\infty}\left(\frac{1}{\sqrt{a+n^2}}-\frac{1}{\sqrt{b+n^2}}\right),\qquad a,b>0,
\]
with \(n_0\in\{0,1\}\). This quantity appears when comparing closely related lattice kernels, removing shared singular behavior, or building regularized observables in periodic models. Although the difference is convergent, direct numerical summation can be inefficient at high accuracy because each component behaves like \(1/n\) and cancellation is essential.

This repository provides a small, embeddable C++ implementation with two methods: a transparent direct summation baseline and a faster transform-based solver. The goal is practical reliability: deterministic stopping criteria, explicit status codes, and stable floating-point handling.

## 2. How to use
Build and run:
```bash
g++ -O3 main.cpp inverse_sqrt_lattice_sum.cpp -o isls
./isls
```

You can also include `inverse_sqrt_lattice_sum.hpp` and call `InverseSqrtLatticeSum::get_value(...)`.

- `method=0` (`ISLS_FAST`): Poisson/integral accelerated evaluation.
- `method=1` (`ISLS_CONVENTIONAL`): direct Kahan-compensated partial sum.
- Inputs must satisfy `a>0`, `b>0`.
- Accuracy/termination are controlled by `atol`, `rtol`, and `max_evaluation`.

## 3. Theory
**Key acceleration idea (high level).**  
Instead of summing many slowly decaying terms directly, separate the problem into a closed-form main part plus a rapidly decaying remainder. The remainder is then integrated with a quadrature that converges much faster than naive term-by-term summation.

**Details (for experts).**  
The difference cancels the leading \(1/n\) asymptotic term, yielding an effective \(O(n^{-3})\) tail, but finite-precision accumulation is still sensitive to cancellation. The conventional path therefore uses Kahan compensation. The fast path applies a Poisson-summation reformulation:
\[
S(a,b)=\tfrac12\log\frac{b}{a}\;\pm\;\tfrac12\left(\frac{1}{\sqrt a}-\frac{1}{\sqrt b}\right)+2\int_0^T\left[\frac{1}{e^{\alpha\cosh t}-1}-\frac{1}{e^{\beta\cosh t}-1}\right]dt,
\]
with \(\alpha=2\pi\sqrt a\), \(\beta=2\pi\sqrt b\), and sign determined by index origin. The integrand decays exponentially in \(t\), enabling efficient trapezoidal refinement. Numerical stability is improved by switching from `1/expm1(x)` to `exp(-x)` in the large-\(x\) regime. Convergence is checked by
\[
|\Delta|\le \mathrm{rtol}\,|S_{\text{old}}|+\mathrm{atol}.
\]

## 4. References
- W. Kahan (1965), compensated summation.
- E. T. Whittaker & G. N. Watson, *A Course of Modern Analysis*.
- NIST Digital Library of Mathematical Functions: https://dlmf.nist.gov/
