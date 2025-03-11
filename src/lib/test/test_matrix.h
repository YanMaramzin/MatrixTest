#pragma once
#include <QObject>
#include <QTest>

class Matrix_test: public QObject
{
    Q_OBJECT

private slots:
    void initTestCase();
    void firstTest();
    void testConstructor();
    void testAddition();
    void testSubtraction();
    void testMultiplication();
    void testTranspose();
    void testDeterminant();
    void testIsSymmetric();
    void testNorm();
    void testTrace();
    void testIsDiagonal();
    void testEigenvalues();
    void testIsOrthogonal();
    void testTriangular();
    void testRank();
    void testGetSet();
    void testInverse();
};
