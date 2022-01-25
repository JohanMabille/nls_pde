#include "PDE_solver.hpp"
#include "BS_solver.h"
#include "Payoff.h"
#include "cmath"
#include<math.h>

int main(int argc, char* argv[])
{


    // INPUT PARAMETERS

    size_t n= 10000;
    double K = 48;
    double vol = 0.16;
    double r = 0.001;
    double S0 = 50;
    double dt = 1.0 / 365;
    double T = 30.0 / 365;
    //double lw_bound = (log(S0) - 5 *vol * pow( T, 0.5));
    //double up_bound = (log(S0) + 5 * vol *pow( T, 0.5));

    // --- Payoff Construction
    dauphine::Call* my_payoff=new dauphine::Call(K, S0); // Call payoff
    //dauphine::Put* my_payoff=new dauphine::Put(K, S0); // Put payoff
    //dauphine::Payoff* my_payoff=new dauphine::Payoff(K,S0,vectorDePayoff) //Custom payoff - possibility to implement any payoff

    // --- Boundary Conditions construction

    // The boundary conditions is an unordered map, you can leave it empty ( i.e not specify boundary conditions)
    // It can take 6 different key possible
            /*= {{'upperValue'      , _double },
                {'upperDerivative'  , _double },
                {'upperDerivative2' , _double },
                {'lowerValue'       , _double },
                {'lowerDerivative'  , _double },
                {'lowerDerivative2' , _double }}
            */
    // If you specify an upper and/or lower value, it will override the upper/lower Derivatives

    // - Exemple /different Boundary conditions
    //std::unordered_map< std::string, double > boundary_conditions ={{"lowerValue", 0 },{"upperValue",25}}; // Exemple
    std::unordered_map< std::string, double > boundary_conditions = { {"lowerValue", 0 },{"upperDerivative2",0} }; //Best Boundary conditions (saw in class)
    //std::unordered_map< std::string, double > boundary_conditions = {};

    // Visualise the Boundary conditions if any
    //std::cout <<"Boundary Conditions  : " << std::endl;
    //dauphine::print_map(boundary_conditions);
    //std::cout <<"\n "<< std::endl;

    // --- Solver BS

    std::cout <<" ----------------- SOLVER BLACK-SCHOLES SIMPLE ----------------- \n " << std::endl;


    dauphine::BS_solver my_BSsolver(vol, r, dt, T, my_payoff, S0,n + 1 , boundary_conditions); // if no boundary_conditions remove the argument in the solver, it will take the default ones
    my_BSsolver.compute_bis();
    my_BSsolver.print_greeks();

    // --- Solver PDE

    //std::cout <<" ----------------- SOLVER PDE ----------------- \n "  << std::endl;

    // INPUT PARAMETERS
      //double a = -0.5 * pow(vol, 2);
    //double b = 0.5 * pow(vol, 2) - r;
    //double c = r;
    //double d = 0;

    //dauphine::PDE_solver my_solver{ a, b, c, d, dt, T, n + 1,  my_payoff, up_bound , lw_bound, boundary_conditions };
    //my_solver.compute_bis();
    //double price = my_solver.get_value(log(S0));
    //std::cout <<"Price PDE :  " << price << std::endl;

    return 0;








}
