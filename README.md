# Project C++ - Implementation of a PDE solver and Application to the Black-Scholes model
C++ Project - M2 203 - Sofia BELARBI - Nicolas DE COSTER - Louis MABILEAU


Our project consists of the following files (code sources and/or headers):

* **main.cpp** : Enter the different input parameters to BS solver and launch it.
* **matrix.cpp** & **matrix.hpp**  : Class  that builds the custom matrix container and all the function needed to acces, manipulate and change the data in it (fill, operator, transpose ...).
* **PDE_solver.cpp** & **PDE_solver.hpp** : Class that builds and computes the whole grid (every points in it) in order to return a single solution (price) in the end. We can enter any paramter we want (a, b, c, d) and choose to specify or not the boundary conditions. The pricer can be recovered for any point in time t and for any price - S0 in particular.
* **BS_solver.cpp** & **BS_solver.h** : Class that inherits from PDE _solver but that return a price within the simple BS framework with financial input parameters ( volatility, rate) as well as the associated greeks.
* **TridiagExtended.cpp** & **TridiagExtended.h** : Class that holds the matrix of our system. It is used to comput a "custom tridiaginal matrix" (with two elements added) The mathematical theory is specified in the associated report.
* **Payoff.cpp** & **Payoff.h** : Class that is used to generate the payoffs of standard and/or custom financial products such as options. We implemented the payoff generator for Calls and Puts. If you wish to generate a custom payoff, you can enter in argument of the Payoff class constructor your vector of payoff. 


