#include "Common.hpp"
#include "gtest/gtest.h"

struct AbelianBgfTest : public CommonBgfTest, public ::testing::Test {

    AbelianBgfTest() {
        std::array<double, 3 * 4> random;
        std::normal_distribution<double> distrib(0, 1);
        for (auto &item : random)
            item = distrib(gen);

        for (auto i{0}; i < 3; ++i) {
            su3B[i * 4] = Cplx(random[4 * i], random[4 * i + 1]);
            bgfB[i] = su3B[i * 4];
            su3C[i * 4] = Cplx(random[4 * i + 2], random[4 * i + 3]);
            bgfC[i] = su3C[i * 4];
            su3One[i * 4] = Cplx(1., 0.);
        }
    };
    SU3 su3B;
    SU3 su3C;
    SU3 su3One;
    bgf::AbelianBgf bgfB;
    bgf::AbelianBgf bgfC;
};

TEST_F(AbelianBgfTest, CopyConstructor) {
    bgf::AbelianBgf Bcopy(bgfB);
    bgf::AbelianBgf Ccopy = bgfC;
    ASSERT_TRUE(Bcopy == bgfB);
    ASSERT_TRUE(Ccopy == bgfC);
    ASSERT_FALSE(Bcopy == Ccopy);
}
TEST_F(AbelianBgfTest, ApplyFromRight) {
    ASSERT_TRUE(SU3Cmp(bgfB.ApplyFromLeft(A), su3B * A)());
}

TEST_F(AbelianBgfTest, ApplyFromLeft) {
    ASSERT_TRUE(SU3Cmp(bgfB.ApplyFromRight(A), A * su3B)());
}

TEST_F(AbelianBgfTest, BgfProduct) {
    SU3 D = (bgfB * bgfC).ApplyFromLeft(su3One);
    ASSERT_TRUE(SU3Cmp(D, su3B * su3C)());
}

TEST_F(AbelianBgfTest, VectorProduct) {
    CVector b;
    CVector c;
    CVector Bb;
    CVector cC;
    std::array<double, 3 * 4> random;
    std::normal_distribution<double> distrib(0, 1);
    for (auto &item : random)
        item = distrib(gen);

    for (int i = 0; i < 3; ++i) {
        c[i] = Cplx(random[4 * i], random[4 * i + 1]);
        b[i] = Cplx(random[4 * i + 2], random[4 * i + 3]);
    }
    cC = bgfC.ApplyFromRight(c);
    Bb = bgfB.ApplyFromLeft(b);
    for (int i = 0; i < 3; ++i) {
        ASSERT_DOUBLE_EQ(cC[i].real(), (bgfC[i] * c[i]).real());
        ASSERT_DOUBLE_EQ(cC[i].imag(), (bgfC[i] * c[i]).imag());
        ASSERT_DOUBLE_EQ(Bb[i].real(), (bgfB[i] * b[i]).real());
        ASSERT_DOUBLE_EQ(Bb[i].imag(), (bgfB[i] * b[i]).imag());
    }
}

TEST_F(AbelianBgfTest, CplxProduct) {
    std::array<double, 2> random;
    std::normal_distribution<double> distrib(0, 1);
    for (auto &item : random)
        item = distrib(gen);

    Cplx alpha;
    (random[0], random[1]);
    SU3 D = (bgfB * alpha).ApplyFromLeft(su3One);
    ASSERT_TRUE(SU3Cmp(D, su3B * alpha)());
}

TEST(AbelianBgf, Arithmetic) {
    bgf::AbelianBgf A{{1.2, 2.3, 3.4}};
    bgf::AbelianBgf B{{0.2, 1.3, 2.4}};

    bgf::AbelianBgf C = A + B;
    bgf::AbelianBgf D = A - B;

    ASSERT_TRUE(Cmp(C[0], 1.4)());
    ASSERT_TRUE(Cmp(C[1], 3.6)());
    ASSERT_TRUE(Cmp(C[2], 5.8)());
    ASSERT_TRUE(Cmp(D[0], 1.)());
    ASSERT_TRUE(Cmp(D[1], 1.)());
    ASSERT_TRUE(Cmp(D[2], 1.)());
}

SU3 makeDiagonal(const Cplx &d0, const Cplx &d1, const Cplx &d2) {
    SU3 matrix;
    matrix[0] = d0;
    matrix[4] = d1;
    matrix[8] = d2;
    return matrix;
}

TEST(AbelianBgf, DISABLED_KnownValues) {
    //  This tests the class
    //     bgf::AbelianBgfFactory
    //  using the formula
    //    V(x0) = V(x0 = 0)* exp{ i a E x0 }
    // Unit 3x3 matrix
    SU3 su3One = makeDiagonal({1., 0.}, {1., 0.}, {1., 0.});
    // pi / 3
    double pio3 = std::atan(1.) * 4. / 3;
    // do 1000 checks

    for (int _n = 0; _n < 1000; ++_n) {
        // random T, L (< 100), eta, nu
        int L = std::uniform_int_distribution<>(0, 32)(gen) + 4;
        int T = L;
        // int T = std::rand() % 100 + 10;
        //  initialize the background field
        bgf::AbelianBgfFactory factory(T, L);
        // random x0
        int x0 = std::uniform_int_distribution<>(0, T)(gen);
        // now calculate
        //       exp (i E t) = dV
        //           = exp ( - i gamma * x0 * diag(2,-1,-1) )
        double gamma = 1. / L / T * pio3;
        SU3 su3dV = makeDiagonal(exp({0, -2. * gamma * x0}), exp({0, gamma * x0}), exp({0, gamma * x0}));
        // get V(x0) and V(0)
        bgf::AbelianBgf V(factory.get(x0));
        bgf::AbelianBgf V0(factory.get(0));

        // std::cout << V0 << "\n"<< V0.ApplyFromLeft(su3dV) << "\n"
        // << V.ApplyFromLeft(su3One) << "\n"
        // << "\n";

        EXPECT_TRUE(SU3Cmp(V0.ApplyFromLeft(su3dV),
                           V.ApplyFromLeft(su3One))());
    }
}

TEST(AbelianBgf, KnownValuesOnAFixedLattice) {
    //  This tests the class
    //     bgf::AbelianBgfFactory
    //  using the formula
    //    V(x0) = V(x0 = 0)* exp{ i a E x0 }
    // Unit 3x3 matrix
    SU3 su3One = makeDiagonal({1., 0.}, {1., 0.}, {1., 0.});
    // pi / 3
    double pio3 = std::atan(1.) * 4. / 3;

    int L = 64;
    int T = L;

    bgf::AbelianBgfFactory factory(T, L);
    // random x0
    for (auto x0 = 0; x0 < L; ++x0) {
        // now calculate
        //       exp (i E t) = dV
        //           = exp ( - i gamma * x0 * diag(2,-1,-1) )
        double gamma = pio3 / (L * T);
        SU3 su3dV = makeDiagonal(exp({0, -2. * gamma * x0}), exp({0, gamma * x0}), exp({0, gamma * x0}));
        // get V(x0) and V(0)
        bgf::AbelianBgf V(factory.get(x0));
        bgf::AbelianBgf V0(factory.get(0));

        EXPECT_TRUE(SU3Cmp(V0.ApplyFromLeft(su3dV),
                           V.ApplyFromLeft(su3One))());
    }
}

// Test tree level formula for E,

TEST(AbelianBgf, E) {
    int L = 8;
    int T = 8;
    // initialize
    bgf::get_abelian_bgf(0, 0, T, L, 0);
    SU3 C = makeDiagonal({0, 1. / L}, {0, -.5 / L}, {0, -.5 / L});
    Cplx E = 0;
    for (int k = 1; k < 4; ++k) {
        E -= (bgf::get_abelian_bgf(0, k) *
              bgf::get_abelian_bgf(0, 0) *
              bgf::get_abelian_bgf(1, k).dag() *
              bgf::get_abelian_bgf(0, 0).dag())
                 .ApplyFromRight(C)
                 .tr();
        E -= (bgf::get_abelian_bgf(T, k).dag() *
              bgf::get_abelian_bgf(T - 1, 0).dag() *
              bgf::get_abelian_bgf(T - 1, k) *
              bgf::get_abelian_bgf(T - 1, 0).dag())
                 .ApplyFromRight(C)
                 .tr();
    }
    E *= L * L * L * 2;
    double pio3 = std::atan(1.) * 4. / 3;
    double k = 12 * L * L * (std::sin(pio3 / L / L) + std::sin(pio3 * 2 / L / L));
    EXPECT_DOUBLE_EQ(E.real(), k);

    L = 8;
    int s = -1;
    T = L - s;
    C = makeDiagonal({0, 1. / L}, {0, -.5 / L}, {0, -.5 / L});
    bgf::AbelianBgfFactory factory(T, L, s);
    E = 0;
    E -= (bgf::AbelianBgf(factory.get(0)) *
          bgf::AbelianBgf(factory.get(1)).dag())
             .ApplyFromRight(C)
             .tr();
    E -= (bgf::AbelianBgf(factory.get(T)).dag() *
          bgf::AbelianBgf(factory.get(T - 1)))
             .ApplyFromRight(C)
             .tr();
    E *= L * L * L * 2 * 3;
    EXPECT_DOUBLE_EQ(E.real(), 37.6945384607827);

    s = 1;
    T = L - s;
    C = makeDiagonal({0, 1. / L}, {0, -.5 / L}, {0, -.5 / L});
    factory = bgf::AbelianBgfFactory(T, L, s);
    E = 0;
    E -= (bgf::AbelianBgf(factory.get(0)) *
          bgf::AbelianBgf(factory.get(1)).dag())
             .ApplyFromRight(C)
             .tr();
    E -= (bgf::AbelianBgf(factory.get(T)).dag() *
          bgf::AbelianBgf(factory.get(T - 1)))
             .ApplyFromRight(C)
             .tr();
    E *= L * L * L * 2 * 3;
    ASSERT_DOUBLE_EQ(E.real(), 37.69169953984173);
}
