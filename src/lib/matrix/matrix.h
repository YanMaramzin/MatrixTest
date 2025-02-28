#pragma once

#include <vector>
#include <cmath>
#include <ostream>

class Matrix {
public:
    // Конструктор для создания матрицы заданного размера
    Matrix(size_t rows, size_t cols);

    // Конструктор для создания матрицы из двумерного вектора
    explicit Matrix(const std::vector<std::vector<double>> &data);

    // Получение количества строк
    size_t rows() const;

    // Получение количества столбцов
    size_t cols() const;

    // Получение элемента матрицы
    double get(size_t row, size_t col) const;

    // Установка элемента матрицы
    void set(size_t row, size_t col, double value);

    // Сложение матриц
    Matrix operator+(const Matrix &other) const;

    // Вычитание матриц
    Matrix operator-(const Matrix &other) const;

    // Умножение матриц
    Matrix operator*(const Matrix &other) const;

    // Транспонирование матрицы
    Matrix transpose() const;

    // Вычисление определителя матрицы
    double determinant() const;

    // Нахождение обратной матрицы
    Matrix inverse() const;

    // Проверка на симметричность
    bool isSymmetric() const;

    // Вычисление следа матрицы
    double trace() const;

    // Проверка на диагональность
    bool isDiagonal() const;

    // Нахождения собственных значений (только для 2x2 матриц)
    std::vector<double> eigenvalues() const;

    // Проверка на ортогональность
    bool isOrthogonal() const;

    // Вычисление нормы матрицы
    double norm() const;

    // Проверка на верхнюю треугольность
    bool isUpperTriangular() const;

    // Проверка на нижнюю треугольность
    bool isLowerTriangular() const;

    // Нахождение ранга матрицы
    size_t rank() const;

    friend std::ostream &operator<<(std::ostream &out, Matrix const &matrix);

private:
    std::vector<std::vector<double>> m_data;
    size_t m_rows {0};
    size_t m_cols {0};

    // Вспомогательная функция для вычисления минора
    Matrix getMinor(size_t row, size_t col) const;

    // Вспомогательная функция для рекурсивного вычисления определителя
    double calculateDeterminant(const Matrix &mat) const;

    // Вспомогательная функция для приведения матрицы к ступенчатому виду
    void reduceToRowEchelonForm(std::vector<std::vector<double>> &mat) const;
};
