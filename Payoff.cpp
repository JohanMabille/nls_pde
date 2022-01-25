#include "Payoff.h"
#include "cmath"
#include <iostream>
#include <ostream>
#include <iterator>
#include <tuple>
#include "math.h"

namespace dauphine
{
    Payoff::Payoff(const double& K, const double& S0, matrix end_values)
        :m_K(K), m_S0(S0),m_end_values(end_values)
    {
    }

    Payoff::Payoff(const double& K, const double& S0)
        : m_K(K), m_S0(S0)
    {
    }

    matrix Payoff::get_payoff(matrix space_grid) const
    {


        return m_end_values;
    }

    Call::Call(const double& K, const double& S0)
        :Payoff(K,S0)
    {
    }

    matrix Call::get_payoff(matrix space_grid) const
	{
        matrix call = maxi(expon(space_grid) - m_K, 0.0);

        return call;
	}

	Put::Put(const double& K, const double& S0)
        :Payoff(K, S0)
        {
        }

    matrix Put::get_payoff(matrix space_grid) const
	{
        matrix put = maxi( m_K - expon(space_grid)  , 0.0);
        return put;
	}

}
