#include "TridiagExtended.h"
#include <algorithm>
#include<math.h>

namespace dauphine
{
	TridiagExtended::TridiagExtended(std::vector<double> d, std::vector<double> l, std::vector<double> u, double e, double f)
		: m_d(d), m_l(l), m_u(u), m_e(e), m_f(f), m_size(d.size())
	{
	}

	std::vector<double>& TridiagExtended::get_d()
	{
		return m_d;
	}

	std::vector<double> TridiagExtended::Solve(std::vector<double> B)
	{
		std::size_t N = m_size;
		std::vector<double> h(N - 1, 0);
		std::vector<double> i(N, 0);
		std::vector<double> j(N - 1, 0);

		// Find LU coefficients
		i[0] = m_d[0];
		j[0] = m_u[0];
		h[0] = m_l[0] / i[0];
		i[1] = m_d[1] - h[0] * j[0];
		j[1] = m_u[1] - h[0] * m_e;

		for (std::size_t m = 1; m < (N - 2); m++)
		{
			h[m] = m_l[m] / i[m];
			i[m + 1] = m_d[m] - h[m] * j[m];
			j[m + 1] = m_u[m];
		}

		double w = m_f / i[N - 3];
		h[N - 2] = (m_l[N - 2] - w * j[N - 3]) / i[N - 2];
		i[N - 1] = m_d[N - 1] - h[N - 2] * j[N - 2];

		//Solve LY=B
		std::vector<double> y(N, 0);
		y[0] = B[0];
		for (std::size_t m = 1; m < (N - 1); m++)
		{
			y[m] = B[m] - y[m - 1] * h[m - 1];
		}
		y[N - 1] = B[N - 1] - h[N - 2] * y[N - 2] - w * y[N - 3];

		//Solve UX=Y
		std::vector<double> x(N, 0);
		x[N - 1] = y[N - 1] / i[N - 1];
		for (std::size_t m = (N - 2); m < (N - 1); m--)
		{
			x[m] = (y[m] - j[m] * x[m + 1]) / i[m];
			x[0] = (y[0] - j[0] * x[1] - m_e * x[2]) / i[0];
		}

		return x;
	}

	std::vector<double> TridiagExtended::prod_withVector(std::vector<double> x)
	{
		std::vector<double> res(m_size, 0);
		res[0] = m_d[0] * x[0] + m_u[0] * x[1] + m_e * x[2];
		for (std::size_t m = 1; m < (m_size - 1); m++)
		{
			res[m] = m_l[m - 1] * x[m - 1] + m_d[m] * x[m] + m_u[m] * x[m + 1];
		}
		res[m_size - 1] = m_f * x[m_size - 3] + m_l[m_size - 2] * x[m_size - 2] + m_d[m_size - 1] * x[m_size - 1];

		return res;

	}

	TridiagExtended operator*(const TridiagExtended& lhs, double factor)
	{
		TridiagExtended res(lhs);
		std::transform(lhs.m_d.begin(), lhs.m_d.end(), res.m_d.begin(), [factor](double arg) {return factor * arg;});
		std::transform(lhs.m_l.begin(), lhs.m_l.end(), res.m_l.begin(), [factor](double arg) {return factor * arg;});
		std::transform(lhs.m_u.begin(), lhs.m_u.end(), res.m_u.begin(), [factor](double arg) {return factor * arg;});
		res.m_e *= factor;
		res.m_f *= factor;

		//std::cout << res.m_d[1] << std::endl;
		//std::cout << res.m_e << std::endl;

		return res;

	}



}
