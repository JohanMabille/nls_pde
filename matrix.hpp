#ifndef MATRIX_HPP_INCLUDED
#define MATRIX_HPP_INCLUDED


#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <ostream>
#include <cstddef>
#include <unordered_map>

namespace dauphine
{

    class matrix
    {
    public:

        matrix(std::size_t nb_rows, std::size_t nb_cols, double value = 0);
        matrix();

        std::size_t nb_rows() const;
        std::size_t nb_cols() const;

        void resize(std::size_t nb_rows, std::size_t nb_cols);

        matrix get_row(std::size_t i);
        matrix get_col(std::size_t j);

        double& operator()(std::size_t i, std::size_t j);
        const double& operator()(std::size_t i, std::size_t j) const;
        matrix operator()(std::size_t i_1 ,std::size_t i_2 , std::size_t j_1 , std::size_t j_2) const;
        std::vector<double> operator()(std::size_t i_1, std::size_t i_2, std::size_t j_1) const; // return a vector

        void fill(std::size_t i_1, std::size_t i_2, std::size_t j_1, std::size_t j_2, const matrix& src);
        void fill(std::size_t i_1, std::size_t i_2, std::size_t j_1, const std::vector<double>& src);

        matrix& operator+=(const matrix& rhs);
        matrix& operator-=(const matrix& rhs);

        matrix& operator+=(double rhs);
        matrix& operator-=(double rhs);
        matrix& operator/=(const matrix& rhs);

        matrix transpose() const;
        matrix inverse() const;


    private:

        std::size_t m_nb_rows;
        std::size_t m_nb_cols;
        std::vector<double> m_data;
    };

    std::ostream& operator<<(std::ostream& out, const matrix& m);

    // Add
    matrix operator+(const matrix& lhs, const matrix& rhs);
    matrix operator+(const matrix& lhs, double rhs);
    matrix operator+(double lhs, const matrix& rhs);

    // Substract
    matrix operator-(const matrix& lhs, const matrix& rhs);
    matrix operator-(const matrix& lhs, double rhs);
    matrix operator-(double lhs, const matrix& rhs);


    // Multiply

    matrix operator*(const matrix& lhs, const matrix& rhs);
    matrix operator*(const matrix& lhs, double rhs);


    // Matrix construction
    matrix linspace(double start, double end, std::size_t nb_step);
    matrix eye(std::size_t size);
    matrix expon( const matrix& m);
    matrix maxi( const matrix& m, double tresh);
    bool is_mat_equal(const matrix& lhs, const matrix& rhs);
    void print_map(std::unordered_map< std::string,  double > &m);

}

#endif // MATRIX_HPP_INCLUDED
