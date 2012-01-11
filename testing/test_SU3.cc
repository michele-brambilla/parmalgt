int PTORD = 6;

#include "gtest/gtest.h"
#include "Helper.h"
#include <MyMath.h>

MyRand r(23797);

TEST(SU3, Multiplication){
  SU3 A, Acpy, alphaA;
  Cplx alpha(r.Rand(), r.Rand());
  for (SU3::iterator ait = A.begin(), 
         aait = alphaA.begin(),
         acpit = Acpy.begin(); 
       ait != A.end(); ++ait, ++aait, ++acpit){
    *ait = Cplx(r.Rand(), r.Rand());
    *acpit = *ait;
    *aait = *ait*alpha;
  }
  Acpy *= alpha;
  EXPECT_TRUE ( SU3Cmp( alphaA, alpha*A)() );
  EXPECT_TRUE ( SU3Cmp( alphaA, A*alpha)() );
  EXPECT_TRUE ( SU3Cmp( alphaA, Acpy)() );
}
