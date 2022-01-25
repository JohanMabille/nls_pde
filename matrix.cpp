#include "matrix.hpp"
#include <algorithm>
#include <vector>
#include <exception>
#include <ostream>
#include<math.h>

namespace dauphine
{
    matrix::matrix(std::size_t nb_rows, std::size_t nb_cols, double value)
        : m_nb_rows(nb_rows),
          m_nb_cols(nb_cols),
          m_data(nb_rows * nb_cols,value)
    {
    }

    matrix::matrix()
        : m_nb_rows(1),
        m_nb_cols(1),
        m_data(1* 1, 0)
    {
    }

    std::size_t matrix::nb_rows() const
    {
        return m_nb_rows;
    }

    std::size_t matrix::nb_cols() const
    {
        return m_nb_cols;
    }

    void matrix::resize(std::size_t nb_rows, std::size_t nb_cols)
    {
        m_nb_rows = nb_rows;
        m_nb_cols = nb_cols;
        m_data.resize(m_nb_rows * m_nb_cols);
    }

     matrix matrix::get_row(std::size_t i)
    {
        if (i > m_nb_rows - 1)
        {
        throw std::logic_error("Error: row not included in matrix");
        }
        matrix res(1, m_nb_cols);
        for (std::size_t j = 0;j < m_nb_cols;j++)
        {
            res(0,j) = m_data[i * m_nb_cols + j];
        }
        return res;
    }
    matrix matrix::get_col(std::size_t j)
    {
        if (j > m_nb_cols - 1)
        {
		throw std::logic_error("Error: column not included in matrix");
		}

        matrix res(m_nb_rows, 1);
        for (std::size_t i = 0;i < m_nb_rows;i++)
        {
            res(i, 0) = m_data[i * m_nb_cols + j];
        }
        return res;
    }

    double& matrix::operator()(std::size_t i, std::size_t j)
    {
        return m_data[i * m_nb_cols + j];
    }

    const double& matrix::operator()(std::size_t i, std::size_t j) const
    {
        return m_data[i * m_nb_cols + j];
    }

    matrix matrix::operator()(std::size_t i_1 ,std::size_t i_2 , std::size_t j_1 , std::size_t j_2) const
    {

		matrix submat {i_2 - i_1 + 1 , j_2 - j_1 +1};

            int subi = 0;
            for (size_t i = i_1; i <= i_2; i++)
            {
                int subj = 0;
                for (size_t j = j_1; j <= j_2; j++)
                {
                   submat(subi,subj) = m_data[i * m_nb_cols + j];

                    subj++;
                }
                subi++;
            }

            return submat;
    }

    std::vector<double> matrix::operator()(std::size_t i_1, std::size_t i_2, std::size_t j_1) const
    {

        std::vector<double> vector(i_2 - i_1 + 1);

        int subi = 0;
        for (size_t i = i_1; i <= i_2; i++)
        {
            vector[subi] = m_data[i * m_nb_cols + j_1];
            subi++;
        }

        return vector;
    }

    void matrix::fill(std::size_t i_1, std::size_t i_2, std::size_t j_1, std::size_t j_2, const matrix& src)
    {
        int subi = 0;
        for (size_t i = i_1; i <= i_2; i++)
        {
            int subj = 0;
            for (size_t j = j_1; j <= j_2; j++)
            {
                m_data[i * m_nb_cols + j]= src(subi, subj) ;

                subj++;
            }
            subi++;
        }
    }

    void matrix::fill(std::size_t i_1, std::size_t i_2, std::size_t j_1, const std::vector<double>& src)
    {
        int subi = 0;
        for (size_t i = i_1; i <= i_2; i++)
        {
            m_data[i * m_nb_cols + j_1] = src[subi];
            subi++;
        }
    }


    std::ostream& operator<<(std::ostream& out, const matrix& m)
    {
        for(std::size_t i = 0; i < m.nb_rows(); ++i)
        {
            for(std::size_t j = 0; j < m.nb_cols(); ++j)
            {
                out << m(i, j) << ", ";
            }
            out << std::endl;
        }
        return out;
    }


    // Add

    matrix& matrix::operator+=(const matrix& rhs)
    {
        if (m_nb_rows != rhs.nb_rows() || m_nb_cols != rhs.nb_cols())
        {
		throw std::logic_error("You can only add matrix of same dimension");
		}
        std::transform(m_data.begin(), m_data.end(), rhs.m_data.begin(),
                       m_data.begin(), std::plus<double>());
        return *this;
    }

    matrix& matrix::operator+=(double rhs)
    {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                       [rhs](double arg) { return arg + rhs; });
        return *this;
    }

    matrix operator+(const matrix& lhs, const matrix& rhs)
    {
        if (lhs.nb_rows() != rhs.nb_rows() && lhs.nb_cols() != rhs.nb_cols())
        {
		throw std::logic_error("You can only add matrix of same dimension");
		}
        matrix tmp(lhs);
        tmp += rhs;
        return tmp;
    }

    matrix operator+(double lhs, const matrix& rhs)
    {
        matrix tmp(rhs);
        tmp += lhs;
        return tmp;
    }

    matrix operator+ (const matrix& rhs,double lhs)
    {
        return lhs + rhs;
    }

    // Substract

    matrix& matrix::operator-=(const matrix& rhs)
    {
        if (m_nb_rows != rhs.nb_rows() && m_nb_cols != rhs.nb_cols())
        {
		throw std::logic_error("You can only substract matrix of same dimension");
		}

		std::transform(m_data.begin(), m_data.end(), rhs.m_data.begin(),
                       m_data.begin(), std::minus<double>());
        return *this;
    }

    matrix& matrix::operator-=(double rhs)
    {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                       [rhs](double arg) { return arg - rhs; });
        return *this;
    }

    matrix operator-(const matrix& lhs, const matrix& rhs)
    {
        if (lhs.nb_rows() != rhs.nb_rows() || lhs.nb_cols() != rhs.nb_cols())
        {
		throw std::logic_error("You can only substract matrix of same dimension");
		}
		matrix tmp(lhs);
        tmp -= rhs;
        return tmp;
    }


    matrix operator- (const matrix& rhs,double lhs)
    {
        matrix tmp(rhs);
        tmp -= lhs;
        return tmp;
    }
    matrix operator- (double lhs, const matrix& rhs)
    {
        matrix tmp(rhs * (-1));
        tmp += lhs;
        return tmp;
    }
    // Multiply

    matrix operator*(const matrix& lhs, const matrix& rhs)
    {
        if (lhs.nb_cols() != rhs.nb_rows() )
        {
		throw std::logic_error("Problem in the dimensions of the matrix you are trying to multiply");
		}
        const size_t c = rhs.nb_cols();
        const size_t r = lhs.nb_rows();

        matrix result(r, c);

        for (size_t i = 0; i < r; i++)
            {
                for (size_t j = 0; j < c; j++)
                {

                for (size_t  k = 0; k < lhs.nb_cols(); k++)
                {
                    result(i,j) +=  lhs(i,k) * rhs(k,j);
                }
            }
        }
        return result;
    }

    matrix operator*(const matrix& lhs, double rhs)
    {
        const size_t c = lhs.nb_cols();
        const size_t r = lhs.nb_rows();

        matrix result(r, c);
        for (size_t i = 0; i < r; i++)
            {
                for (size_t j = 0; j < c; j++)
                {
                    result(i,j) =  lhs(i,j) * rhs;
                }
            }
        return result;
    }

    //Divide
    matrix& matrix::operator/=(const matrix& rhs)
    {
    matrix m=*this;
	for (std::size_t i = 0; i < m_nb_rows; i++) {
		matrix row = m.get_row(i);
		for (std::size_t j = 0; j < m_nb_cols; j++) {
			row.get_col(j) /= rhs;
		}
	}
	return *this;
    }

    //Transpose
    matrix matrix::transpose() const

    {
        matrix m(*this);
        const size_t c = m_nb_cols;
        const size_t r = m_nb_rows;

        matrix result(c,r);
        for (size_t i = 0; i < r; i++)
            {
                for (size_t j = 0; j < c; j++)
                {
                    result(j,i) = m(i,j);
                }
            }


        return result;
    }


    // Matrix construstuction
    matrix linspace(double start, double end, std::size_t nb_step)
    {
        double step = (end - start) / (nb_step - 1);
        matrix res(nb_step, 1);
        res(0, 0) = start;
        for (std::size_t i = 1;i < nb_step;i++)
        {
            res(i, 0) = res(i - 1, 0) + step;
        }
        return res;
    }

    // Identity matrix
    matrix eye(std::size_t size)
    {
        matrix res(size, size);
        for (std::size_t i = 0;i < size;i++)
        {
            res(i, i) = 1;
        }
        return res;
    }



    //Inverse

   matrix matrix::inverse() const
    {

    if ( m_nb_rows != m_nb_cols)
        {
		throw std::logic_error("Can't invert non-square matrix");
        }

	std::size_t n = m_nb_rows;
	matrix m(*this);
	matrix id = eye(n);
	float d = 0.0;
	float **mat = NULL;
	mat = new float*[2*n];

    for (std::size_t i = 0; i < 2*n; ++i)
    {
        mat[i] = new float[2*n]();
    }

	matrix im{n,n};


	for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {

             mat[i][j]= m(i,j);

        }
        for (std::size_t j = n; j < n*2; ++j)
        {
            mat[i][j]= id(i,j-n);
        }
    }


    // Partial pivoting
    for(std::size_t i = n; i > 1; --i)
    {
        if(mat[i-1][1] < mat[i][1])
        {
            for(std::size_t j = 0; j < 2*n; ++j)
            {
                d = mat[i][j];

                mat[i][j] = mat[i-1][j];

                mat[i-1][j] = d;
            }
        }
    }

    // Reducing To Diagonal Matrix
    for(std::size_t  i = 0; i < n; ++i)
    {
        for( std::size_t j = 0; j < n*2; ++j)
        {
            if(j != i)
            {
                if ( mat[i][i] == 0)
                {
                throw std::logic_error("Matrix Can't be inverted - division by 0 ");
                }

                d = mat[j][i] / mat[i][i];
                for(std::size_t k = 0; k < n*2; ++k)
                {
                    mat[j][k] -= mat[i][k]*d;
                }
            }
        }
    }


    // Reducing To Unit Matrix
    for(std::size_t i = 0; i < n; ++i)
    {
        d = mat[i][i];
        for(std::size_t j = 0; j < 2*n; ++j)
        {
            mat[i][j] = mat[i][j]/d;
        }
    }


    for(std::size_t i=0; i < n; ++i)
    {
        for(std::size_t j = 0; j < n; ++j)
        {
          im(i,j)=mat[i][j+n];
        }

    }

    // Deleting the memory allocated
    for (std::size_t i = 0; i < n; ++i)
    {
        delete[] mat[i];
    }
    delete[] mat;

    return im;
    }

    matrix expon(const matrix& m)
    {
        std::size_t rows = m.nb_rows();
        std::size_t cols = m.nb_cols();
        matrix mat (rows, cols);

        for(std::size_t i=0; i < rows; ++i)
        {
            for(std::size_t j = 0; j < cols; ++j)
            {
              mat(i,j)=exp(m(i,j));
            }

        }
        return mat;
    }

    matrix maxi( const matrix& m, double tresh)

    {
        std::size_t rows = m.nb_rows();
        std::size_t cols = m.nb_cols();
        matrix mat(rows, cols);

        for(std::size_t i=0; i < rows; ++i)
        {
            for(std::size_t j = 0; j < cols; ++j)
            {

              if (m(i,j) > tresh )
              {
                  mat(i,j) = m(i,j);
              }
              else
              {
                  mat(i,j) = tresh;
              }
            }
        }
        return mat;

    }

    // Equality

    bool is_mat_equal(const matrix& lhs, const matrix& rhs)
    {
        bool res = true;
        if (lhs.nb_rows() == rhs.nb_rows() && lhs.nb_cols() == rhs.nb_cols())
        {
            for (std::size_t i = 0; i < lhs.nb_rows(); ++i)
            {
                for (std::size_t j = 0; j < lhs.nb_cols(); ++j)
                {
                    if (lhs(i, j) != rhs(i, j))
                    {
                        res = false;
                    }
                }
            }
        }
        else
        {
            res = false;
        }
        return res;
    }

    void print_map(std::unordered_map< std::string,  double > &m)
    {
        for (auto const &pair: m) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
        }
    }


}
