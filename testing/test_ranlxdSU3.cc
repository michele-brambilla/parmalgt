#include "gtest/gtest.h"
#include <MyMath.h>
#include <MyRand.h>
#include <ranlxd.hpp>

const int NN = 1000000;

TEST(Norm, Old){
  MyRand mr;
  mr.init(123);
  SU3 A, AA;
  for (int i = 0; i < NN; ++i){
    SU3 tmp = SU3rand(mr);
    A += tmp;
    AA += tmp*tmp;
  }
  std::cout << A / NN << std::endl;
  std::cout << AA / NN << std::endl;
}
TEST(Norm, New){
  ranlxd::Rand rl(1,123);
  SU3  B, BB;// = SU3rand(mr), B = SU3rand(rl);
  for (int i = 0; i < NN; ++i){
    SU3 tmp = SU3rand(rl);
    B += tmp;
    BB += tmp*tmp;
  }
  std::cout << B / NN << std::endl;
  std::cout << BB / NN << std::endl;
}
