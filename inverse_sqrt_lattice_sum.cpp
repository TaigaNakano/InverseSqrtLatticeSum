#include "inverse_sqrt_lattice_sum.hpp"

#if defined(__clang__)
  #define NOOPT [[clang::optnone]]
#elif defined(__GNUC__)
  #define NOOPT [[gnu::optimize("O0")]]
#else
  #define NOOPT
#endif

NOOPT
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

    int _eval_count = 0;
    double old_s, nn, y, t, delta;

    while (max_evaluation < 0 || _eval_count < max_evaluation)
    {
        old_s = s;

        nn = (double)k;
        nn *= nn;

        y = (1.0 / std::sqrt(a + nn) - 1.0 / std::sqrt(b + nn)) - c;
        t = s + y;
        c = (t - s) - y;
        s = t;

        _eval_count++;

        delta = s - old_s;
        if (std::fabs(delta) <= rtol * std::fabs(old_s) + atol)
        {
            *eval_count = _eval_count;
            *result = s;
            return ISLS_COMPUTATION_CONVERGED;
        }

        k++;
    }

    *eval_count = _eval_count;
    *result = s;
    return ISLS_COMPUTATION_REACHED_MAXEVAL;
}

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

    // 積分区間が実質ゼロなら pref のみ
    if (T < InverseSqrtLatticeSum::DEFAULT_EPS)
    {
        *eval_count = 0;
        *result = pref;
        return ISLS_COMPUTATION_NO_ITERATION_NEEDED;
    }

    // 端点評価に最低2回必要
    if (max_evaluation >= 0 && max_evaluation < 2)
    {
        *eval_count = 0;
        *result = pref;
        return ISLS_COMPUTATION_REACHED_MAXEVAL;
    }

    double g0 = integrant_from_cosh(1.0, alpha, beta);          // cosh(0)=1
    double gT = integrant_from_cosh(std::cosh(T), alpha, beta); // cosh(T)
    int _eval_count = 2;

    double Tk = 0.5 * (g0 + gT) * T;
    int n = 1;

    while (true)
    {
        n <<= 1;
        int half = n >> 1;

        // max_evaluation を厳密に超えない
        if (max_evaluation >= 0 && _eval_count + half > max_evaluation)
            break;

        double h = T / n;

        // 内点列 t = (2j-1)h
        double delta2h = 2.0 * h;
        double coshDelta = std::cosh(delta2h);
        double sinhDelta = std::sinh(delta2h);

        double cosh_t = std::cosh(h);
        double sinh_t = std::sinh(h);

        double sum_mid = 0.0;
        for (int j = 0; j < half; j++)
        {
            sum_mid += integrant_from_cosh(cosh_t, alpha, beta);

            // t <- t + 2h
            double newCosh = cosh_t * coshDelta + sinh_t * sinhDelta;
            double newSinh = sinh_t * coshDelta + cosh_t * sinhDelta;
            cosh_t = newCosh;
            sinh_t = newSinh;
        }

        _eval_count += half;

        double Tk2 = 0.5 * Tk + h * sum_mid;
        double delta = Tk2 - Tk;

        if (std::fabs(delta) <= rtol * std::fabs(Tk) + atol)
        {
            *eval_count = _eval_count;
            *result = pref + 2.0 * Tk2;
            return ISLS_COMPUTATION_CONVERGED;
        }

        Tk = Tk2;

        if (max_evaluation >= 0 && _eval_count == max_evaluation) break;
    }

    *eval_count = _eval_count;
    *result = pref + 2.0 * Tk;
    return ISLS_COMPUTATION_REACHED_MAXEVAL;
}

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