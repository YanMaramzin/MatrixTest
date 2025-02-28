#include "Matrix.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>


Matrix::Matrix(const size_t rows, const size_t cols) : m_rows(rows), m_cols(cols)
{
    if (m_rows == 0 || m_cols == 0)
        throw std::invalid_argument("Matrix dimensions cannot be zero");

    m_data.resize(m_rows, std::vector<double>(m_rows, 0.0));
}

Matrix::Matrix(const std::vector<std::vector<double> > &data)
{
    if (m_data.empty() || m_data[0].empty())
        throw std::invalid_argument("Matrix m_data cannot be empty");

    m_data = data;
    m_rows = data.size();
    m_cols = data[0].size();
}

size_t Matrix::rows() const
{
    return m_rows;
}

size_t Matrix::cols() const
{
    return m_cols;
}

double Matrix::get(const size_t row, const size_t col) const
{
    if (row >= m_rows || col >= m_cols)
        throw std::out_of_range("Matrix indices out of range");

    return m_data[row][col];
}

void Matrix::set(const size_t row, const size_t col, const double value)
{
    if (row >= m_rows || col >= m_cols)
        throw std::out_of_range("Matrix indices out of range");

    m_data[row][col] = value;
}

Matrix Matrix::operator+(const Matrix &other) const
{
    if (m_rows != other.m_rows || m_cols != other.m_cols)
        throw std::invalid_argument("Matrices dimensions must match for addition");

    Matrix result(m_rows, m_cols);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_cols; ++j)
            result.m_data[i][j] = m_data[i][j] + other.m_data[i][j];
    }

    return result;
}

Matrix Matrix::operator-(const Matrix &other) const
{
    if (m_rows != other.m_rows || m_cols != other.m_cols)
        throw std::invalid_argument("Matrices dimensions must match for subtraction");

    Matrix result(m_rows, m_cols);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_cols; ++j)
            result.m_data[i][j] = m_data[i][j] - other.m_data[i][j];
    }
    return result;
}

Matrix Matrix::operator*(const Matrix &other) const
{
    if (m_cols != other.m_rows)
        throw std::invalid_argument("Matrices dimensions must match for multiplication");

    Matrix result(m_rows, other.m_cols);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < other.m_cols; ++j) {
            for (size_t k = 0; k < m_cols; ++k)
                result.m_data[i][j] += m_data[i][k] * other.m_data[k][j];
        }
    }

    return result;
}

Matrix Matrix::transpose() const
{
    Matrix result(m_cols, m_rows);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_cols; ++j)
            result.set(j, i, m_data[i][j]);
    }

    return result;
}

double Matrix::determinant() const
{
    if (m_rows != m_cols)
        throw std::invalid_argument("Matrix must be square to calculate determinant");

    return calculateDeterminant(*this);
}

Matrix Matrix::inverse() const
{
    if (m_rows != m_cols)
        throw std::invalid_argument("Matrix must be square to calculate inverse");

    const auto det = determinant();
    if (std::fabs(det) < 1e-9)
        throw std::invalid_argument("Matrix is singular, cannot compute inverse");

    Matrix adjugate(m_rows, m_cols);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_cols; ++j) {
            Matrix minor = getMinor(i, j);
            const auto cofactor = ((i + j) % 2 == 0 ? 1 : -1) * minor.determinant();
            adjugate.set(j, i, cofactor / det);
        }
    }
    return adjugate;
}

bool Matrix::isSymmetric() const
{
    if (m_rows != m_cols)
        return false;

    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_cols; ++j) {
            if (m_data[i][j] != m_data[j][i])
                return false;
        }
    }

    return true;
}

double Matrix::trace() const
{
    if (m_rows != m_cols)
        throw std::invalid_argument("Matrix must be square to calculate trace");

    auto trace = 0.0;
    for (size_t i = 0; i < m_rows; ++i)
        trace += m_data[i][i];

    return trace;
}

bool Matrix::isDiagonal() const
{
    if (m_rows != m_cols)
        return false;

    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_cols; ++j) {
            if (i != j && std::fabs(m_data[i][j]) > 1e-9)
                return false;
        }
    }

    return true;
}

std::vector<double> Matrix::eigenvalues() const
{
    if (m_rows != 2 || m_cols != 2)
        throw std::invalid_argument("Eigenvalues can only be calculated for 2x2 matrices");

    const auto a = m_data[0][0];
    const auto b = m_data[0][1];
    const auto c = m_data[1][0];
    const auto d = m_data[1][1];
    const auto trace = a + d;
    const auto det = a * d - b * c;
    const auto discriminant = trace * trace - 4 * det;
    if (discriminant < 0)
        throw std::invalid_argument("Complex eigenvalues are not supported");

    auto sqrtDiscriminant = std::sqrt(discriminant);
    return {(trace + sqrtDiscriminant) / 2, (trace - sqrtDiscriminant) / 2};
}

bool Matrix::isOrthogonal() const
{
    if (m_rows != m_cols)
        return false;

    const auto product = (*this) * this->transpose();
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_cols; ++j) {
            if (i == j && std::fabs(product.m_data[i][j] - 1.0) > 1e-9)
                return false;

            if (i != j && std::fabs(product.m_data[i][j]) > 1e-9)
                return false;
        }
    }
    return true;
}

double Matrix::norm() const
{
    auto sum = 0.0;
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_cols; ++j)
            sum += m_data[i][j] * m_data[i][j];
    }
    return std::sqrt(sum);
}

bool Matrix::isUpperTriangular() const
{
    if (m_rows != m_cols)
        return false;

    for (size_t i = 1; i < m_rows; ++i) {
        for (size_t j = 0; j < i; ++j) {
            if (std::fabs(m_data[i][j]) > 1e-9)
                return false;
        }
    }
    return true;
}

bool Matrix::isLowerTriangular() const
{
    if (m_rows != m_cols)
        return false;

    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = i + 1; j < m_cols; ++j) {
            if (std::fabs(m_data[i][j]) > 1e-9)
                return false;
        }
    }

    return true;
}

size_t Matrix::rank() const
{
    std::vector<std::vector<double> > mat = m_data;
    reduceToRowEchelonForm(mat);
    size_t rank = 0;
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_cols; ++j) {
            if (std::fabs(mat[i][j]) > 1e-9) {
                rank++;
                break;
            }
        }
    }
    return rank;
}

Matrix Matrix::getMinor(size_t row, size_t col) const
{
    Matrix minor(m_rows - 1, m_cols - 1);
    for (size_t i = 0, mi = 0; i < m_rows; ++i) {
        if (i == row)
            continue;

        for (size_t j = 0, mj = 0; j < m_cols; ++j) {
            if (j == col) continue;
            minor.m_data[mi][mj] = m_data[i][j];
            ++mj;
        }
        ++mi;
    }
    return minor;
}

double Matrix::calculateDeterminant(const Matrix &mat) const
{
    if (mat.m_rows == 1)
        return mat.m_data[0][0];

    if (mat.m_rows == 2)
        return mat.m_data[0][0] * mat.m_data[1][1] - mat.m_data[0][1] * mat.m_data[1][0];

    auto det = 0.0;
    for (size_t j = 0; j < mat.m_cols; ++j) {
        Matrix minor = mat.getMinor(0, j);
        const auto cofactor = ((0 + j) % 2 == 0 ? 1 : -1) * minor.determinant();
        det += mat.m_data[0][j] * cofactor;
    }
    return det;
}

void Matrix::reduceToRowEchelonForm(std::vector<std::vector<double> > &mat) const
{
    size_t lead = 0;
    for (size_t r = 0; r < m_rows; ++r) {
        if (lead >= m_cols)
            return;

        size_t i = r;
        while (std::fabs(mat[i][lead]) < 1e-9) {
            ++i;
            if (i == m_rows) {
                i = r;
                ++lead;
                if (lead == m_cols)
                    return;
            }
        }
        std::swap(mat[i], mat[r]);
        const auto lv = mat[r][lead];
        for (size_t j = 0; j < m_cols; ++j)
            mat[r][j] /= lv;

        for (size_t i = 0; i < m_rows; ++i) {
            if (i != r) {
                const auto lv = mat[i][lead];
                for (size_t j = 0; j < m_cols; ++j) {
                    mat[i][j] -= lv * mat[r][j];
                }
            }
        }
        ++lead;
    }
}

std::ostream &operator<<(std::ostream &out, Matrix const &matrix)
{
    for (size_t i = 0; i < matrix.rows(); ++i) {
        out << "[";
        for (size_t j = 0; j < matrix.cols(); ++j) {
            if (j != matrix.cols() - 1)
                out << matrix.m_data[i][j] << " ";
            else
                out << matrix.m_data[i][j] << "]";
        }
        if (i != matrix.rows() - 1)
            out << "\n";
    }
    return out;
}
