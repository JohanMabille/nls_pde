#include "PDE_solver.hpp"
#include "cmath"
#include <iostream>
#include <ostream>
#include <iterator>
#include <tuple>
#include<math.h>


namespace dauphine
{
    PDE_solver::PDE_solver(double a, double b, double c, double d, double dt, double T, std::size_t n, Payoff* payoff,
    double upper_bound, double lower_bound,
    std::unordered_map< std::string,  double > boundary_conditions, double theta ):
        m_a(a),
        m_b(b),
        m_c(c),
        m_d(d),
        m_dt(dt),
        m_T(T),
        m_theta(theta),
        m_N(n),
        m_system_size(m_N),
        m_upper_bound(upper_bound),
        m_lower_bound(lower_bound),
        m_dx((upper_bound - lower_bound) / (m_N-1)),
        m_values(matrix(n, std::size_t(round(T/dt)+1))),
        m_space_grid(linspace(lower_bound, upper_bound, n)),
        m_grid_computed(false),
        m_UV_present(false),
        m_LV_present(false),
        m_LD1_nset(1),
        m_LD2_nset(1),
        m_UD1_nset(1),
        m_UD2_nset(1),
        m_boundary_conditions(boundary_conditions),
        m_payoffPtr(payoff),
        m_end_values(m_N,1)
    {
        m_end_values = m_payoffPtr->get_payoff(m_space_grid);
        m_values.fill(0, m_N-1, round(m_T / m_dt), round(m_T / m_dt), m_end_values);

        if (m_boundary_conditions.find("lowerValue") != m_boundary_conditions.end())
        {
            double LV = m_boundary_conditions["lowerValue"];
            for (std::size_t i=0; i < round(T/dt) ; i++)
            {
                m_values(0,i)= LV;
            }
            m_system_size -=1;
            m_LV_present = true;
        }

        if (m_boundary_conditions.find("upperValue") != m_boundary_conditions.end())
        {
            double UV = m_boundary_conditions["upperValue"];

            for (std::size_t i=0; i < round(T/dt) ; i++)
            {
                m_values(m_N-1,i)= UV;
            }
            m_system_size -=1;
            m_UV_present = true;
        }
        if (m_boundary_conditions.find("lowerDerivative") != m_boundary_conditions.end())
        {
            m_LD1_nset = 0;
        }
        if (m_boundary_conditions.find("lowerDerivative2") != m_boundary_conditions.end())
        {
            m_LD2_nset = 0;
        }
        if (m_boundary_conditions.find("upperDerivative") != m_boundary_conditions.end())
        {
            m_UD1_nset = 0;
        }
        if (m_boundary_conditions.find("upperDerivative2") != m_boundary_conditions.end())
        {
            m_UD2_nset = 0;
        }

    }

    PDE_solver::~PDE_solver()
    {
        delete m_payoffPtr;
    }

    double PDE_solver::compute_bottom_alpha() const
    {
        return m_a/pow(m_dx, 2)*m_LD2_nset - m_b/m_dx*m_LD1_nset + m_c;
    }

    double PDE_solver::compute_bottom_beta() const
    {
        return -2*m_a/pow(m_dx, 2) * m_LD2_nset + m_b/m_dx * m_LD1_nset;
    }

    double PDE_solver::compute_bottom_gamma() const
    {
        return m_a/pow(m_dx, 2) * m_LD2_nset;
    }

    double PDE_solver::compute_alpha() const
    {
        return m_a/pow(m_dx, 2) - m_b/(2*m_dx);
    }

    double PDE_solver::compute_beta() const
    {
        return m_c - 2*m_a/pow(m_dx, 2);
    }

    double PDE_solver::compute_gamma() const
    {
        return m_a/pow(m_dx,2) + m_b/(2*m_dx);
    }

    double PDE_solver::compute_top_alpha() const
    {
        return m_a/pow(m_dx, 2) * m_UD2_nset;
    }

    double PDE_solver::compute_top_beta() const
    {
        return -2*m_a/pow(m_dx, 2) * m_UD2_nset - m_b/m_dx * m_UD1_nset;
    }

    double PDE_solver::compute_top_gamma() const
    {
        return m_a/pow(m_dx, 2) * m_UD2_nset + m_b/m_dx * m_UD1_nset + m_c;
    }



    double PDE_solver::get_value(double S_l) const

    {

        if (!m_grid_computed)
        {
            throw std::logic_error("Error: Grid not computed yet");
        }
        if (S_l < m_lower_bound)
        {
            throw std::logic_error("Error: Value outside the grid (too low)");
        }
        else if (S_l > m_upper_bound)
        {
            throw std::logic_error("Error: Value outside the grid (too high)");
        }


        double closest_index = 0.;
        double min_gap = 1000.;
        for (size_t i = 0; i < m_N; i++)
        {
            if (abs(S_l - m_space_grid(i, 0)) < min_gap)
            {
                min_gap = abs(S_l - m_space_grid(i, 0));
                closest_index = i;
            }
        }

        return m_values(closest_index, 0);

    }




    TridiagExtended PDE_solver::compute_A_bis(double alpha, double beta, double gamma) const
    {
        std::vector<double> l(m_system_size - 1, alpha);
        std::vector<double> d(m_system_size, beta);
        std::vector<double> u(m_system_size - 1, gamma);
        double e;
        double f;

        if (!m_LV_present) {
            d[0] = compute_bottom_alpha();
            u[0] = compute_bottom_beta();
            e = compute_bottom_gamma();
        }
        else {
            e = 0;
        };

        if (!m_LV_present) {
            f = compute_top_alpha();
            l[m_system_size - 2] = compute_top_beta();
            d[m_system_size - 1] = compute_top_gamma();
        }
        else {
            f = 0;
        };

        TridiagExtended A = TridiagExtended(d, l, u, e, f);

        return A;

    }

    TridiagExtended PDE_solver::compute_B_bis(const TridiagExtended& A) const
    {
        TridiagExtended B(A);
        B = B * (-(1 - m_theta) * m_dt);
        std::vector<double> d = B.get_d();
        std::transform(d.begin(), d.end(), d.begin(), [](double arg) {return 1+ arg;});
        B.get_d() = d;
        return B;
    }

    std::vector<double> PDE_solver::compute_D_bis(double alpha, double gamma) const
    {
        std::vector<double> D(m_system_size, m_d);

        if (m_LV_present) {
            D[0] += alpha * m_boundary_conditions.at("lowerValue");
        }
        else {
            if (m_boundary_conditions.find("lowerDerivative") != m_boundary_conditions.end()) {
                D[0] += m_b * m_boundary_conditions.at("lowerDerivative");
            }
            if (m_boundary_conditions.find("lowerDerivative2") != m_boundary_conditions.end()) {
                D[0] += m_a * m_boundary_conditions.at("lowerDerivative2");
            }
        }

        if (m_UV_present) {
            D[m_system_size - 1] += gamma * m_boundary_conditions.at("upperValue");
        }
        else {
            if (m_boundary_conditions.find("upperDerivative") != m_boundary_conditions.end()) {
                D[m_system_size - 1] += m_b * m_boundary_conditions.at("upperDerivative");
            }
            if (m_boundary_conditions.find("upperDerivative2") != m_boundary_conditions.end()) {
                D[m_system_size - 1] += m_a * m_boundary_conditions.at("upperDerivative2");
            }
        }

        return D;
    }

    void PDE_solver::compute_bis()
    {
        double alpha = compute_alpha();
        double beta = compute_beta();
        double gamma = compute_gamma();

        TridiagExtended A = compute_A_bis(alpha, beta, gamma);
        TridiagExtended B = compute_B_bis(A);
        std::vector<double> D = compute_D_bis(alpha, gamma);
        std::vector<double> E(D.size());
        std::transform(D.begin(), D.end(), E.begin(), [this](double arg) {return -m_dt*arg;});
        A = A * (m_theta * m_dt);
        std::transform(A.get_d().begin(), A.get_d().end(), A.get_d().begin(), [](double arg) {return 1 + arg;});

        double lower_index;
        double upper_index;
        if (m_LV_present) {
            lower_index = 1;
        }
        else {
            lower_index = 0;
        }
        if (m_UV_present) {
            upper_index = m_N - 2;
        }
        else {
            upper_index = m_N - 1;
        }

        for (std::size_t i = round(m_T / m_dt) - 1; i < round(m_T / m_dt); i--)
        {
            std::vector <double> F = B.prod_withVector(m_values(lower_index, upper_index, i + 1));
            std::transform(F.begin(), F.end(), E.begin(), F.begin(), std::plus<double>());
            std::vector <double> X = A.Solve(F);
            m_values.fill(lower_index, upper_index, i, X);
        }

        m_grid_computed = true;


    }


};
