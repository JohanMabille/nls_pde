#pragma once
#include <vector>
#include <iostream>
#include <ostream>

namespace dauphine
{
    class TridiagExtended
    {
    public:
        TridiagExtended(std::vector<double> d, std::vector<double> l, std::vector<double> u, double e, double f);
        std::vector<double> Solve(std::vector<double> B);
        std::vector<double> prod_withVector(std::vector<double> x);
        std::vector<double>& get_d();
        

    private:
        std::vector<double> m_d;
        std::vector<double> m_l;
        std::vector<double> m_u;
        double m_e;
        double m_f;
        std::size_t m_size;

        friend TridiagExtended operator*(const TridiagExtended& lhs, double factor);
    };

    


};