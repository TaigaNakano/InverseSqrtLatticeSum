#include "inverse_sqrt_lattice_sum.hpp"

#if defined(__clang__)
  #define NOOPT [[clang::optnone]]
#elif defined(__GNUC__)
  #define NOOPT [[gnu::optimize("O0")]]
#else
  #define NOOPT
#endif

NOOPT
/**
 * @brief Directly evaluates the truncated lattice sum with Kahan compensated summation.
 */
int InverseSqrtLatticeSum::partialsum_with_kahans_summation(
    double* result,
    int* eval_count,
    double a,
    double b,
    bool start_from_one,
    double atol,
    double rtol,
    int max_evaluation)
{
    if (!result || !eval_count) return ISLS_COMPUTATION_INVALID_INPUT;
    if (a <= 0.0 || b <= 0.0) return ISLS_COMPUTATION_INVALID_INPUT;

    double s = 0.0;
    double c = 0.0; // Kahan compensation

    int k = start_from_one ? 1 : 0;

    int eval_count_local = 0;
    double old_s, nn, y, t, delta;

    while (max_evaluation < 0 || eval_count_local < max_evaluation)
    {
        old_s = s;

        nn = (double)k;
        nn *= nn;

        y = (1.0 / std::sqrt(a + nn) - 1.0 / std::sqrt(b + nn)) - c;
        t = s + y;
        c = (t - s) - y;
        s = t;

        eval_count_local++;

        delta = s - old_s;
        if (std::fabs(delta) <= rtol * std::fabs(old_s) + atol)
        {
            *eval_count = eval_count_local;
            *result = s;
            return ISLS_COMPUTATION_CONVERGED;
        }

        k++;
    }

    *eval_count = eval_count_local;
    *result = s;
    return ISLS_COMPUTATION_REACHED_MAXEVAL;
}

/**
 * @brief Evaluates the transformed integral form derived from Poisson summation.
 */
int InverseSqrtLatticeSum::poissons_summation_with_integral(
    double* result,
    int* eval_count,
    double a,
    double b,
    bool start_from_one,
    double atol,
    double rtol,
    int max_evaluation)
{
    if (!result || !eval_count) return ISLS_COMPUTATION_INVALID_INPUT;
    if (a <= 0.0 || b <= 0.0) return ISLS_COMPUTATION_INVALID_INPUT;

    double sa = std::sqrt(a);
    double sb = std::sqrt(b);

    double alpha = 2.0 * InverseSqrtLatticeSum::PI * sa;
    double beta  = 2.0 * InverseSqrtLatticeSum::PI * sb;

    double pref = 0.5 * std::log(b / a);
    pref += start_from_one
        ? -0.5 * (1.0 / sa - 1.0 / sb)
        : +0.5 * (1.0 / sa - 1.0 / sb);

    double m = std::fmin(alpha, beta);
    double r = InverseSqrtLatticeSum::X_SWITCH_FOR_INVEXPM1 / m;
    double T = (r <= 1.0) ? 0.0 : std::acosh(r);

    /** @brief If the integration interval is effectively zero, return only the prefactor term. */
    if (T < InverseSqrtLatticeSum::DEFAULT_EPS)
    {
        *eval_count = 0;
        *result = pref;
        return ISLS_COMPUTATION_NO_ITERATION_NEEDED;
    }

    /** @brief At least two evaluations are required for endpoint sampling. */
    if (max_evaluation >= 0 && max_evaluation < 2)
    {
        *eval_count = 0;
        *result = pref;
        return ISLS_COMPUTATION_REACHED_MAXEVAL;
    }

    double g0 = integrant_from_cosh(1.0, alpha, beta);          // cosh(0)=1
    double gT = integrant_from_cosh(std::cosh(T), alpha, beta); // cosh(T)
    int eval_count_local = 2;

    double Tk = 0.5 * (g0 + gT) * T;
    int n = 1;

    while (true)
    {
        n <<= 1;
        int half = n >> 1;

        /** @brief Never exceed max_evaluation strictly. */
        if (max_evaluation >= 0 && eval_count_local + half > max_evaluation)
            break;

        double h = T / n;

        /** @brief Midpoint sequence: t = (2j-1)h. */
        double delta2h = 2.0 * h;
        double coshDelta = std::cosh(delta2h);
        double sinhDelta = std::sinh(delta2h);

        double cosh_t = std::cosh(h);
        double sinh_t = std::sinh(h);

        double sum_mid = 0.0;
        for (int j = 0; j < half; j++)
        {
            sum_mid += integrant_from_cosh(cosh_t, alpha, beta);

            /** @brief Increment t by 2h using hyperbolic-angle addition identities. */
            double newCosh = cosh_t * coshDelta + sinh_t * sinhDelta;
            double newSinh = sinh_t * coshDelta + cosh_t * sinhDelta;
            cosh_t = newCosh;
            sinh_t = newSinh;
        }

        eval_count_local += half;

        double Tk2 = 0.5 * Tk + h * sum_mid;
        double delta = Tk2 - Tk;

        if (std::fabs(delta) <= rtol * std::fabs(Tk) + atol)
        {
            *eval_count = eval_count_local;
            *result = pref + 2.0 * Tk2;
            return ISLS_COMPUTATION_CONVERGED;
        }

        Tk = Tk2;

        if (max_evaluation >= 0 && eval_count_local == max_evaluation) break;
    }

    *eval_count = eval_count_local;
    *result = pref + 2.0 * Tk;
    return ISLS_COMPUTATION_REACHED_MAXEVAL;
}

/**
 * @brief Dispatches computation to the selected method implementation.
 */
int InverseSqrtLatticeSum::get_value(
    double* result,
    int* eval_count,
    double a,
    double b,
    int method,
    bool start_from_one,
    double atol,
    double rtol,
    int max_evaluation)
{
    if (!result || !eval_count) return ISLS_COMPUTATION_INVALID_INPUT;

    switch (method)
    {
    case ISLS_FAST:
        return poissons_summation_with_integral(result, eval_count, a, b, start_from_one, atol, rtol, max_evaluation);
    case ISLS_CONVENTIONAL:
        return partialsum_with_kahans_summation(result, eval_count, a, b, start_from_one, atol, rtol, max_evaluation);
    default:
        return poissons_summation_with_integral(result, eval_count, a, b, start_from_one, atol, rtol, max_evaluation);
    }
}
