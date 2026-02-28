#include <iostream>
#include <iomanip>
#include "inverse_sqrt_lattice_sum.hpp"

int main()
{
    std::cout << std::scientific;
    std::cout << std::setprecision(15);
    int eval_count = 0.0;
    double value = 0.0;

    double a = 0.0, b = 0.0;
    int method = 0;    
    int use_input_start_from_one; bool start_from_one = false;
    double atol, rtol;
    int max_eval;
    std::cout << "Please input value of a" << std::endl;
    std::cout << "> ";
    std::cin >> a;
    std::cout << "Please input value of b" << std::endl;
    std::cout << "> ";
    std::cin >> b;
    std::cout << "Please input computation method(0: Poisson's sum, 1: Direct computation)" << std::endl;
    std::cout << "> ";
    std::cin >> method;
    std::cout << "Please input index start from zero or one(1: One, other wise 0)" << std::endl;
    std::cout << "> ";
    std::cin >> use_input_start_from_one;
    std::cout << "Please input absolute tolerance" << std::endl;
    std::cout << "> ";
    std::cin >> atol;
    std::cout << "Please input relative tolerance" << std::endl;
    std::cout << "> ";
    std::cin >> rtol;
    std::cout << "Please input maximum iteration" << std::endl;
    std::cout << "> ";
    std::cin >> max_eval;

    start_from_one = use_input_start_from_one == 1;
    int result = InverseSqrtLatticeSum::get_value(&value, &eval_count, a, b, method, start_from_one, atol, rtol, max_eval);
    
    std::cout << "--- Computation result ---" << std::endl;
    switch (result)
    {
        case ISLS_COMPUTATION_INVALID_INPUT:
            std::cout << "Reuslt: Invalid input (a, b) = (" << a << ", " << b << ") " << std::endl;
            break;
        
        case ISLS_COMPUTATION_CONVERGED:
            std::cout << " -Reuslt: Computation was converged" << std::endl;
            break;

        case ISLS_COMPUTATION_REACHED_MAXEVAL:
            std::cout << " -Reuslt: Computation reached max iteration" << std::endl;
            break;

        case ISLS_COMPUTATION_NO_ITERATION_NEEDED:
            std::cout << " -Reuslt: No iteration needed" << std::endl;
            break;

        default:
            break;
    }
    std::cout << " -Value of ISLS = " << value << std::endl;
    std::cout << " -Counts of evaluation = " << eval_count << std::endl;
    std::cout << "--- Computation result END ---" << std::endl;
    return 0;
}