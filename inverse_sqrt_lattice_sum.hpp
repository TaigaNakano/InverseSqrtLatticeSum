#pragma once

#include <cmath>

/** @brief Return code indicating invalid input arguments. */
#define ISLS_COMPUTATION_INVALID_INPUT -1
/** @brief Return code indicating successful convergence. */
#define ISLS_COMPUTATION_CONVERGED 0
/** @brief Return code indicating that the maximum evaluation count was reached. */
#define ISLS_COMPUTATION_REACHED_MAXEVAL 1
/** @brief Return code indicating that no numerical iteration was required. */
#define ISLS_COMPUTATION_NO_ITERATION_NEEDED 2

/** @brief Fast method based on Poisson summation with numerical integration. */
#define ISLS_FAST 0
/** @brief Conventional direct partial-sum method with Kahan summation. */
#define ISLS_CONVENTIONAL 1

/**
 * @brief Utility class for computing inverse square-root lattice sums.
 */
class InverseSqrtLatticeSum
{
    private:
        static constexpr double X_SWITCH_FOR_INVEXPM1 = 36.0436533891171560897;
        static constexpr double PI = 3.141592653589793238462643383279502884;

        InverseSqrtLatticeSum() = delete;
        ~InverseSqrtLatticeSum() = delete;

        /**
         * @brief Compute the lattice sum by direct accumulation with Kahan compensation.
         *
         * @param[out] result Computed sum value.
         * @param[out] eval_count Number of evaluated terms.
         * @param a Positive parameter in the summand.
         * @param b Positive parameter in the summand.
         * @param start_from_one If true, start summation from index 1; otherwise from 0.
         * @param atol Absolute convergence tolerance.
         * @param rtol Relative convergence tolerance.
         * @param max_evaluation Maximum number of term evaluations (negative means unlimited).
         * @return Status code defined by ISLS_COMPUTATION_* macros.
         */
        static int partialsum_with_kahans_summation(
            double* result,
            int* eval_count,
            double a,
            double b,
            bool start_from_one,
            double atol,
            double rtol,
            int max_evaluation);

        /**
         * @brief Compute the lattice sum using Poisson summation and trapezoidal integration.
         *
         * @param[out] result Computed sum value.
         * @param[out] eval_count Number of integrand evaluations.
         * @param a Positive parameter in the summand.
         * @param b Positive parameter in the summand.
         * @param start_from_one If true, apply correction for index start at 1.
         * @param atol Absolute convergence tolerance.
         * @param rtol Relative convergence tolerance.
         * @param max_evaluation Maximum number of integrand evaluations (negative means unlimited).
         * @return Status code defined by ISLS_COMPUTATION_* macros.
         */
        static int poissons_summation_with_integral(
            double* result,
            int* eval_count,
            double a,
            double b,
            bool start_from_one,
            double atol,
            double rtol,
            int max_evaluation);

        /**
         * @brief Evaluate 1/(exp(x)-1) for positive x with a stable large-x branch.
         * @param x Positive real value.
         * @return Stable evaluation of 1/(exp(x)-1).
         */
        static inline double inv_expm1_pos(double x)
        {
            return (x >= InverseSqrtLatticeSum::X_SWITCH_FOR_INVEXPM1) ? std::exp(-x) : 1.0 / std::expm1(x);
        }

        /**
         * @brief Evaluate the integrand expressed by hyperbolic cosine.
         * @param coshT Value of cosh(t).
         * @param alpha First scaled parameter.
         * @param beta Second scaled parameter.
         * @return Integrand value at t.
         */
        static inline double integrand_from_cosh(double coshT, double alpha, double beta)
        {
            return inv_expm1_pos(alpha * coshT) - inv_expm1_pos(beta * coshT);
        }

    public:
        static constexpr double DEFAULT_EPS = 2.220446049250313e-16;
        static constexpr int DEFAULT_MAXEVAL = 2147483647;

        /**
         * @brief Public API for computing the inverse square-root lattice sum.
         *
         * @param[out] result Computed sum value.
         * @param[out] eval_count Number of evaluations consumed by the selected method.
         * @param a Positive parameter in the summand.
         * @param b Positive parameter in the summand.
         * @param method Computation method (ISLS_FAST or ISLS_CONVENTIONAL).
         * @param start_from_one If true, start index at 1; otherwise at 0.
         * @param atol Absolute convergence tolerance.
         * @param rtol Relative convergence tolerance.
         * @param max_evaluation Maximum evaluations allowed (negative means unlimited).
         * @return Status code defined by ISLS_COMPUTATION_* macros.
         */
        static int get_value(
            double* result,
            int* eval_count,
            double a,
            double b,
            int method = 0,
            bool start_from_one = false,
            double atol = 0.0,
            double rtol = InverseSqrtLatticeSum::DEFAULT_EPS,
            int max_evaluation = InverseSqrtLatticeSum::DEFAULT_MAXEVAL);
};
