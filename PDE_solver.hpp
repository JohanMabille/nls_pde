#ifndef PDE_SOLVER_HPP_INCLUDED
#define PDE_SOLVER_HPP_INCLUDED

#include "matrix.hpp"
#include <iostream>
#include <unordered_map>
#include <string>
#include <utility>
#include "TridiagExtended.h"
#include "Payoff.h"

namespace dauphine
{
    class PDE_solver
    {
    public:
        PDE_solver(double a, double b, double c, double d, double dt, double T, std::size_t n,
            Payoff* payoff, double upper_bound, double lower_bound,
        std::unordered_map< std::string, double > boundary_conditions = {} , double theta = 0.5);

        ~PDE_solver();

        double compute_bottom_alpha() const;
        double compute_bottom_beta() const;
        double compute_bottom_gamma() const;
        double compute_alpha() const;
        double compute_beta() const;
        double compute_gamma() const;
        double compute_top_alpha() const;
        double compute_top_beta() const;
        double compute_top_gamma() const;

        TridiagExtended compute_A_bis(double alpha, double beta, double gamma) const;
        TridiagExtended compute_B_bis(const TridiagExtended& A) const;
        std::vector<double> compute_D_bis(double alpha, double gamma) const;
        void compute_bis();

        double get_value(double S_l) const;

    protected:
        double m_a;
        double m_b;
        double m_c;
        double m_d;
        double m_dt;
        double m_T;
        double m_theta;
        std::size_t m_N;
        int m_system_size;
        matrix m_end_values;
        double m_upper_bound;
        double m_lower_bound;
        double m_dx;
        matrix m_values;
        matrix m_space_grid;
        bool m_grid_computed;
        bool m_UV_present;
        bool m_LV_present;
        int m_LD1_nset;
        int m_LD2_nset;
        int m_UD1_nset;
        int m_UD2_nset;
        std::unordered_map< std::string,  double> m_boundary_conditions;
        Payoff* m_payoffPtr;



    };
};

#endif // PDE_SOLVER_HPP_INCLUDED
