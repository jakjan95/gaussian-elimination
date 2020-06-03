/*
A program to solving linear equations using gaussian elimination.

classical_gaussian_elimination() is solving linear equations in simple 
way by just zeroes elements under pivot, but when it meet zero on pivot it crashes.

gaussian_elimination_with_partial_pivot() is solving linear equations 
by swapping rows of matix and zeroes elements under pivot.

both of this functions using back_substitution() which directly 
calculate variable value by back substitutions

*/

#include <iostream>
#include <numeric>
#include <vector>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

//---------------------------------------------- Classical gaussian elimination
void classical_elimination(Matrix& A, Vector& b) {
    const int row_size = A[0].size();

    for (int j = 0; j < row_size - 1; ++j) {
        const double pivot = A[j][j];

        if (pivot == 0) {
            throw std::runtime_error{"Elim failure"};
        }

        for (int i = j + 1; i < row_size; ++i) {
            const double mult = A[i][j] / pivot;

            for (int k = j; k < A[j].size(); ++k) {
                A[i][k] += A[j][k] * -mult;
            }

            b[i] -= mult * b[j];
        }
    }
}

Vector back_substitution(const Matrix& A, const Vector& b) {
    const int row_size = A[0].size();
    Vector x(row_size);

    for (int i = row_size - 1; i >= 0; --i) {
        //I tried std::inner_product here but it wasn't working correctly
        double dot_product = 0;
        for (int j = i + 1; j < A[i].size(); ++j) {
            dot_product += A[i][j] * x[j];
        }

        double s = b[i] - dot_product;

        if (double m = A[i][i]) {
            x[i] = s / m;
        } else {
            throw std::runtime_error{"Back substitution failure"};
        }
    }
    return x;
}

Vector classical_gaussian_elimination(Matrix A, Vector b) {
    classical_elimination(A, b);
    return back_substitution(A, b);
}

//-------------------------- Gaussian elimination using patrial pivot
void elim_with_partial_pivot(Matrix& A, Vector& b) {
    const int row_size = A[0].size();

    for (int j = 0; j < row_size; ++j) {
        int pivot_row = j;

        for (int k = j + 1; k < row_size; ++k) {
            if (std::abs(A[k][j]) > std::abs(A[pivot_row][j])) {
                pivot_row = k;
            }
        }

        if (pivot_row != j) {
            std::swap(A[j], A[pivot_row]);
            std::swap(b[j], b[pivot_row]);
        }

        for (int i = j + 1; i < row_size; ++i) {
            const double pivot = A[j][j];

            if (pivot == 0) {
                throw std::runtime_error{"can't calculate: pivot == 0"};
            }

            const double mult = A[i][j] / pivot;
            for (int k = j; k < A[j].size(); ++k) {
                A[i][k] += A[j][k] * -mult;
            }

            b[i] -= mult * b[j];
        }
    }
}

Vector gaussian_elimination_with_partial_pivot(Matrix A, Vector b) {
    elim_with_partial_pivot(A, b);
    return back_substitution(A, b);
}

//----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    os << "{\n";
    for (const auto& el : m) {
        os << "{ ";
        for (const auto& item : el) {
            os << item << ' ';
        }
        os << "}\n";
    }
    os << '}';
    return os;
}

std::ostream& operator<<(std::ostream& os, const Vector& v) {
    os << "{ ";
    for (const auto& el : v) {
        os << el << ' ';
    }
    os << '}';
    return os;
}

//----------------------------------------------------------------------------
Vector operator*(const Matrix& m, const Vector& u) {
    const int row_size = m[0].size();
    Vector v(row_size);
    for (int i = 0; i < row_size; ++i) {
        double dot_product = 0;
        for (int j = 0; j < row_size; ++j) {
            dot_product += m[i][j] * u[j];
        }
        v[i] = dot_product;
    }
    return v;
}

//----------------------------------------------------------------------------
void solve_matrix(Matrix A, Vector b) {
    std::cout << "***Classical gausian elimination***\n";
    std::cout << A << '\n';
    std::cout << b << '\n';

    try {
        Vector x = classical_gaussian_elimination(A, b);
        std::cout << "Solution of classic gaussian elimination is:\nx = " << x << '\n';
        Vector v = A * x;
        std::cout << "A * x = " << v << '\n';
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
    }
}

//----------------------------------------------------------------------------
void solve_matrix_with_central(Matrix A, Vector b) {
    std::cout << "***Classical gausian elimination using patrial pivot***\n";
    std::cout << A << '\n';
    std::cout << b << '\n';

    try {
        Vector x = gaussian_elimination_with_partial_pivot(A, b);
        std::cout << "Solution of classic gaussian elimination with partial pivot is:\nx = " << x << '\n';
        Vector v = A * x;
        std::cout << "A * x = " << v << '\n';
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
    }
}

int main() {
    Matrix A1 = {{11.9, 99.1, 0.01}, {15.9, 19.2, 2}, {11.9, 25.1, 13}};
    Vector b1 = {31, 0.05, 91};

    solve_matrix(A1, b1);
    solve_matrix_with_central(A1, b1);

    Matrix A2 = {{0, 1}, {1, 0}};
    Vector b2 = {2, 3};

    solve_matrix(A2, b2);
    solve_matrix_with_central(A2, b2);
}
