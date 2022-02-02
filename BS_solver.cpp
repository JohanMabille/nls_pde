#include "BS_solver.h"
#include "cmath"
#include<math.h>
#include <iostream>

namespace dauphine
{
	BS_solver::BS_solver(double vol, double r, double dt, double T, Payoff* payoff, double S0, std::size_t n,
		double upper_bound, double lower_bound, std::unordered_map< std::string, double > boundary_conditions, double theta) :
		m_S0(S0),
		m_vol(vol),
		m_r(r),
		PDE_solver(-0.5 * pow(vol, 2), 0.5 * pow(vol, 2) - r, r, 0, dt, T, n,payoff, upper_bound, lower_bound, boundary_conditions, theta)
	{

	};

	BS_solver::BS_solver(double vol, double r, double dt, double T, Payoff* payoff, double S0, std::size_t n,
		std::unordered_map< std::string, double > boundary_conditions, double theta) :
		m_S0(S0),
		m_vol(vol),
		m_r(r),
		PDE_solver(-0.5 * pow(vol, 2), 0.5 * pow(vol, 2) - r, r, 0, dt, T,n, payoff, log(S0) + 5 * vol * pow( T, 0.5), log(S0) - 5 *vol * pow( T, 0.5), boundary_conditions, theta)
	{

	};

	double BS_solver::get_price() const
	{
		return get_value(log(m_S0));
	}

	double BS_solver::get_price(double S) const
	{
	    if (log(S) > log(m_S0) + 5 * m_vol * pow( m_T, 0.5) || log(S) < log(m_S0) - 5 *m_vol * pow( m_T, 0.5) )
        {
		throw std::logic_error("S out of range - please renter an S within grid bound");
		}
		return get_value(log(S));
	}

	double BS_solver::closest_ind(double S) const
	{

		double closest_index = 0;
		double min_gap = 100;

		for (size_t i = 0; i < m_N; i++)
		{
			if (abs(log(S) - m_space_grid(i, 0)) < min_gap)
			{
				min_gap = abs(log(S) - m_space_grid(i, 0));
				closest_index = i;
			}
		}
		return closest_index;
	}

	double BS_solver::get_delta(double S) const
	{
	    if (log(S) > log(m_S0) + 5 * m_vol * pow( m_T, 0.5) || log(S) < log(m_S0) - 5 *m_vol * pow( m_T, 0.5) )
        {
		throw std::logic_error("S out of range - please renter an S within grid bound");
		}
		double closest_index = closest_ind(S);

		double price_up = m_values(closest_index + 1, 0);
		double price_down = m_values(closest_index - 1, 0);
		double deltaS = exp(m_space_grid(closest_index + 1, 0)) - exp(m_space_grid(closest_index - 1, 0));
		double delta = (price_up - price_down) / deltaS;

		return delta;
	}

	double BS_solver::get_delta() const
	{
		return get_delta(m_S0);
	}

	double BS_solver::get_gamma(double S) const
	{
	    if (log(S) > log(m_S0) + 5 * m_vol * pow( m_T, 0.5) || log(S) < log(m_S0) - 5 *m_vol * pow( m_T, 0.5) )
        {
		throw std::logic_error("S out of range - please renter an S within grid bound");
		}
		double closest_index = closest_ind(S);

		double price_up2 = m_values(closest_index + 2, 0);
		double price_curr = m_values(closest_index, 0);
		double price_down2 = m_values(closest_index - 2, 0);
		double deltaS = exp(m_space_grid(closest_index + 1, 0)) - exp(m_space_grid(closest_index - 1, 0));
		double gamma = (price_up2 + price_down2 - 2 * price_curr) / pow(deltaS, 2);
        return gamma;
	}

	double BS_solver::get_gamma() const
	{
		return get_gamma(m_S0);
	}

	double BS_solver::get_theta(double S) const
	{
	    if (log(S) > log(m_S0) + 5 * m_vol * pow( m_T, 0.5) || log(S) < log(m_S0) - 5 *m_vol * pow( m_T, 0.5) )
        {
		throw std::logic_error("S out of range - please renter an S within grid bound");
		}

		double closest_index = closest_ind(S);

		double theta = (m_values(closest_index, 1) - m_values(closest_index, 0)) / m_dt;
		return theta;
	}

	double BS_solver::get_theta() const
	{
		return get_theta(m_S0);
	}

	double BS_solver::get_vega(double S, double vol_bump) const
	//Ressource heavy -> only check the impact of a positive vol bump
	{
	    if (log(S) > log(m_S0) + 5 * m_vol * pow( m_T, 0.5) || log(S) < log(m_S0) - 5 *m_vol * pow( m_T, 0.5) )
        {
		throw std::logic_error("S out of range - please renter an S within grid bound");
		}
		double price = get_price(S);
                // Here you pass m_payoffPtr to ta new BS_solver object. YOu copy the pointer, not the pointer object.
                // Since PDE_Solver, the mother class of BS_solver, is reponsible for the deletion of this pointer
                // you will have a double deletion (one from BSsolverBumped, and one form this).
                // A better solution would be to allocate the payoff on the stack in the main function, pass it
                // by pointer (it's fine), but not delete it in the solver. When exiting the main function, the payoff
                // will be automagically deleted.
		BS_solver BSsolverBumped(m_vol + vol_bump, m_r, m_dt, m_T, m_payoffPtr, m_S0, m_N, m_upper_bound, m_lower_bound, m_boundary_conditions);
		BSsolverBumped.compute_bis();
		double price_bumped = BSsolverBumped.get_price(S);

		double vega = (price_bumped - price) / vol_bump;
        return vega;
	}

	double BS_solver::get_vega(double vol_bump) const
	{
		return get_vega(m_S0, vol_bump);
	}

	void BS_solver::print_greeks(double S, double vol_bump)
	{
	    if (log(S) > log(m_S0) + 5 * m_vol * pow( m_T, 0.5) || log(S) < log(m_S0) - 5 *m_vol * pow( m_T, 0.5) )
        {
		throw std::logic_error("S out of range - please renter an S within grid bound");
		}

		double price = get_price(S);
		std::cout << "Price BS :" << price << std::endl;
		double delta = get_delta(S);
        std::cout << "Delta : " << delta * 100 << "%" << std::endl;
        double gamma = get_gamma(S);
		std::cout << "Gamma : " << gamma * 100 << "%" << std::endl;;
		double theta = get_theta(S);
		std::cout << "Theta :" << theta << std::endl;
		double thetaNormalized = theta / 365;
		std::cout << "Theta normalized :" << thetaNormalized * 100 << "%" << std::endl;
		double vega = get_vega(S, vol_bump);
		std::cout << "Vega :" << vega << std::endl;

	}

	void BS_solver::print_greeks( double vol_bump)
	{

		double price = get_price(m_S0);
		std::cout << "Price BS :" << price << std::endl;
		double delta = get_delta(m_S0);
		std::cout << "Delta : " << delta * 100 << "%" << std::endl;
		double gamma = get_gamma(m_S0);
		std::cout << "Gamma : " << gamma * 100 << "%" << std::endl;
		double theta = get_theta(m_S0);
		std::cout << "Theta :" << theta << std::endl;
		double thetaNormalized = theta / 365;
		std::cout << "Theta normalized :" << thetaNormalized * 100 << "%" << std::endl;
		double vega = get_vega(m_S0, vol_bump);
		std::cout << "Vega :" << vega << std::endl;

	}

};
