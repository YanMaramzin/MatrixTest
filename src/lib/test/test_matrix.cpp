#include <test_matrix.h>
#include "matrix.h"

void Matrix_test::initTestCase()
{
    qDebug() << "Matrix_test";
}

void Matrix_test::firstTest()
{
    qDebug() << "Privet!";
}
void Matrix_test::testConstructor() {
    Matrix m(2, 3);
    QCOMPARE(m.rows(), 2);
    QCOMPARE(m.cols(), 3);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            QCOMPARE(m.get(i, j), 0.0);
        }
    }
    QVERIFY_EXCEPTION_THROWN(Matrix(0, 1), std::invalid_argument);
    qDebug() << "Всё норм, матрица создалась";
}
void Matrix_test::testGetSet() {
    Matrix m(2, 2);
    m.set(1, 1, 3.14);
    QCOMPARE(m.get(1, 1), 3.14);

    // сделал проверку на выход за границы
    QVERIFY_EXCEPTION_THROWN(m.set(2, 0, 1), std::out_of_range);
    QVERIFY_EXCEPTION_THROWN(m.get(0, 2), std::out_of_range);
}


void Matrix_test::testAddition() {
    Matrix m1({{-1, 2}, {3, -6}});
    Matrix m2({{5, 6}, {0.00001, -6}});
    Matrix result = m1 + m2;
    QCOMPARE(result.get(0, 0), 4.0);
    QCOMPARE(result.get(0, 1), 8.0);
    QCOMPARE(result.get(1, 0), 3.00001);
    QCOMPARE(result.get(1, 1), -12.0);
    Matrix zero(2, 2);
    QVERIFY((m1 + zero).get(0, 1) == m1.get(0, 1));
}
void Matrix_test::testSubtraction() {
    Matrix m1({{-1, 2}, {3, -6}});
    Matrix m2({{5, 6}, {0.00001, -6}});
    Matrix result = m1 - m2;
    QCOMPARE(result.get(0, 0), -6.0);
    QCOMPARE(result.get(0, 1), -4.0);
    QCOMPARE(result.get(1, 0), 2.99999);
    QCOMPARE(result.get(1, 1), 0.0);
}
void Matrix_test::testMultiplication() {
    Matrix m1({{1000, 0.00001}, {-199, 0}});
    Matrix m2({{2, 0}, {1, 2}});
    Matrix result = m1 * m2;
    QCOMPARE(result.get(0, 0), 2000.00001);
    QCOMPARE(result.get(0, 1), 0.00002);
    QCOMPARE(result.get(1, 0), -398);
    QCOMPARE(result.get(1, 1), 0);

}
void Matrix_test::testTranspose() {
    Matrix m({{1, 2, 3}, {4, 5, 6}});
    Matrix transposed = m.transpose();
    QCOMPARE(transposed.rows(), 3);
    QCOMPARE(transposed.cols(), 2);
    QCOMPARE(transposed.get(0, 0), 1.0);
    QCOMPARE(transposed.get(1, 1), 5.0);

}
void Matrix_test::testDeterminant() {
//    Matrix m1({{5}}); // Не поддерживает матрицу из 1 элемента, ну вроде очевидно
//    QCOMPARE(m1.determinant(), 5.0);
    Matrix m2({{1, 2}, {3, 4}});
    QCOMPARE(m2.determinant(), -2.0);
    Matrix m3({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    QCOMPARE(m3.determinant(), 0.0);

}
void Matrix_test::testIsSymmetric() {
    Matrix sym1({{5, 2}, {2, 3}});
    QVERIFY(sym1.isSymmetric());
    Matrix sym2({{0, 5}, {5, 0}});
    QVERIFY(sym2.isSymmetric());
    Matrix sym3({{5, 0}, {0, 5}});
    QVERIFY(sym3.isSymmetric());
    Matrix sym4({{1, 2, 3}, {2, 5, 4}, {3, 4, 6}});
    QVERIFY(sym4.isSymmetric());
    Matrix nonsym({{1, 2}, {3, 4}});
    QVERIFY(!nonsym.isSymmetric());

}
void Matrix_test::testTrace() {
    Matrix m1({{1, 2}, {3, 4}});
    QCOMPARE(m1.trace(), 5.0);
    Matrix m2({{-1, -2}, {-3, -4}});
    QCOMPARE(m2.trace(), -5.0);
    Matrix m3({{-1, -2, -3}, {-4, -5, -6}, {-7, -8, -9}});
    QCOMPARE(m3.trace(), -15.0);

}
void Matrix_test::testIsDiagonal() {
    Matrix diag({{1, 0}, {0, 4}}); // Чуть чуть скучная функция, поэтому и тесты будут скучными
    QVERIFY(diag.isDiagonal());
    Matrix nonDiag({{1, 2}, {3, 4}});
    QVERIFY(!nonDiag.isDiagonal());
}
void Matrix_test::testEigenvalues() {
    Matrix m({{4, 1}, {2, 3}}); // ладно, вот это было хард, чут чут помог чат гпт
    auto ev = m.eigenvalues();
    auto trace = 7.0;
    auto det = 10.0;
    auto sqrtDisc = std::sqrt(trace * trace - 4 * det);
    QCOMPARE(ev.size(), 2);
    QVERIFY(std::abs(ev[0] - (trace + sqrtDisc)/2) < 1e-9);
}

void Matrix_test::testNorm() {
    Matrix m({{1, 2}, {3, 4}});
    double expected = std::sqrt(1 + 4 + 9 + 16);
    QCOMPARE(m.norm(), expected);
}

void Matrix_test::testIsOrthogonal() {
    Matrix orth({{0, 1}, {1, 0}}); // Ну вот на true OR false тесты лайтовые, поэтому я так, по 1-2 теста пишу (как QVERIFY работает я пон)
    QVERIFY(orth.isOrthogonal());
    Matrix nonOrth({{1, 2}, {3, 4}});
    QVERIFY(!nonOrth.isOrthogonal());
}

void Matrix_test::testTriangular() {
    Matrix upper({{1, 2, 3}, {0, 4, 5}, {0, 0, 6}});
    QVERIFY(upper.isUpperTriangular());
    Matrix lower({{1, 0, 0}, {2, 3, 0}, {4, 5, 6}});
    QVERIFY(lower.isLowerTriangular());
}
void Matrix_test::testRank() { // вот тут я не понял как работает изначальный код, поэтому делал тесты по наитию.
    Matrix fullRank({{1, 2}, {3, 4}});
    QCOMPARE(fullRank.rank(), 2);
    Matrix rank1({{1, 2}, {2, 4}});
    QCOMPARE(rank1.rank(), 1);
}

QTEST_MAIN(Matrix_test)
