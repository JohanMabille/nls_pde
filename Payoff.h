#ifndef PAYOFF_H_INCLUDED
#define PAYOFF_H_INCLUDED

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include "matrix.hpp"

namespace dauphine
{
    class Payoff
    {
        public:
            Payoff(const double& K, const double& S0,matrix end_values);
            virtual ~Payoff() = default;
            Payoff& operator=(const Payoff&) = delete;
            Payoff(const Payoff&&) = delete;
            Payoff& operator=(Payoff&) = delete;
            virtual matrix get_payoff(matrix space_grid) const;


        protected:
            matrix m_end_values;
            const double& m_K;
            const double& m_S0;
            Payoff(const double& K, const double& S0);


        };

        class Call : public Payoff
        {
        public:
            Call( const double& K, const double& S0) ;
            ~Call() = default;

            matrix get_payoff(matrix space_grid) const;

        };

        class Put : public Payoff
        {
        public:
            Put( const double& K, const double& S0) ;
            ~Put() = default;

            matrix get_payoff(matrix space_grid) const;

        };

};


#endif // PAYOFF_H_INCLUDED
