#ifndef BS_SOLVER_HPP_INCLUDED
#define BS_SOLVER_HPP_INCLUDED

#include "PDE_solver.hpp"


namespace dauphine
{
	class BS_solver : public PDE_solver
	{
	public:
		BS_solver(double vol, double r, double dt, double T, Payoff* payoff, double S0, std::size_t n,
			double upper_bound, double lower_bound, std::unordered_map< std::string, double > boundary_conditions = {}, double theta = 0.5);
		BS_solver(double vol, double r, double dt, double T, Payoff* payoffs, double S0, std::size_t n,
			std::unordered_map< std::string, double > boundary_conditions = {}, double theta = 0.5);

		double get_price(double S) const;
		double get_price() const;
		double closest_ind(double S) const;
		double get_delta(double S) const;
		double get_delta() const;
		double get_gamma(double S) const;
		double get_gamma() const;
		double get_theta(double S) const;
		double get_theta() const;
		double get_vega(double S, double vol_bump = 0.0001) const;
		double get_vega(double vol_bump = 0.0001) const;
		void print_greeks(double S, double vol_bump = 0.0001);
		void print_greeks(double vol_bump = 0.0001);


	protected:
		double m_S0;
		double m_vol;
		double m_r;


	};



};


#endif // BS_SOLVER_HPP_INCLUDED



