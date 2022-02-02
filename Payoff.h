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
            //missing copy constructor, move constructor and move assign
            //signature are incorrect
            Payoff(const Payoff&) = delete;
            Payoff& operator=(const Payoff&) = delete;
            Payoff(/*const*/ Payoff&&) = delete;
            Payoff& operator=(Payoff&&) = delete;
            // Why passing the space_gfrid argument by copy instead of const ref?
            // a matrix object can be costly to copy
            virtual matrix get_payoff(matrix space_grid) const;
            

        protected:
            matrix m_end_values;
            // Why storing these by const ref instead of values?
            // a double is not a big object, copy is cheap
            const double& m_K;
            const double& m_S0;
            Payoff(const double& K, const double& S0);

            
        };

        class Call : public Payoff
        {
        public:
            Call( const double& K, const double& S0) ;
            // No need to declare a default destructor here since
            // the destructor of the base class has been declared
            // as virtual
            ~Call() = default;

            // You could use the 'override' keyword to explicitly
            // state that this method is an overload of a method
            // declared in the mother class
            matrix get_payoff(matrix space_grid) const;

        };

        class Put : public Payoff
        {
        public:
            // Same remarks here
            Put( const double& K, const double& S0) ;
            ~Put() = default;

            matrix get_payoff(matrix space_grid) const;

        };

};


#endif // PAYOFF_H_INCLUDED
