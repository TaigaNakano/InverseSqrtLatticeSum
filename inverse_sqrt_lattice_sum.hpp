#pragma once

#include <cmath>

// return codes
#define ISLS_COMPUTATION_INVALID_INPUT -1
#define ISLS_COMPUTATION_CONVERGED 0
#define ISLS_COMPUTATION_REACHED_MAXEVAL 1
#define ISLS_COMPUTATION_NO_ITERATION_NEEDED 2

// computation method
#define ISLS_FAST 0
#define ISLS_CONVENTIONAL 1

class InverseSqrtLatticeSum
{
    private:
        static constexpr double X_SWITCH_FOR_INVEXPM1 = 36.0436533891171560897;
        static constexpr double PI = 3.141592653589793238462643383279502884;

        InverseSqrtLatticeSum() = delete;
        ~InverseSqrtLatticeSum() = delete;

        static int partialsum_with_kahans_summation(
            double* result,
            int* eval_count,
            double a,
            double b,
            bool start_from_one,
            double atol,
            double rtol,
            int max_evaluation);

        static int poissons_summation_with_integral(
            double* result,
            int* eval_count,
            double a,
            double b,
            bool start_from_one,
            double atol,
            double rtol,
            int max_evaluation);

        // inline helpers (header-only)
        static inline double inv_expm1_pos(double x)
        {
            return (x >= InverseSqrtLatticeSum::X_SWITCH_FOR_INVEXPM1) ? std::exp(-x) : 1.0 / std::expm1(x);
        }

        static inline double integrant_from_cosh(double coshT, double alpha, double beta)
        {
            return inv_expm1_pos(alpha * coshT) - inv_expm1_pos(beta * coshT);
        }

    public:
        static constexpr double DEFAULT_EPS = 2.220446049250313e-16;
        static constexpr int DEFAULT_MAXEVAL = 2147483647;

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