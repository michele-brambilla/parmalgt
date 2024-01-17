#include "gtest/gtest.h"

#include <kernels/Kernels.hpp>

#include <fstream>

static constexpr int DIM = 4;
static constexpr int L = 4;

TEST(Kernels, trivial_test){
// EXPECT(false);
}



std::ostream& operator<<(std::ostream& os, const SU3& v) {
    for(int row=0;row<3;++row) {
        for(int col=0;col<3;++col) {
            std::cout << "(" << v(row, col).real() << "," << v(row, col).imag() << ")\t";
        }
        std::cout << "\n";   
    }
}

TEST(Kernels, test_Staple_square_kernel){
    using bgf_t = bgf::ScalarBgf;
    using field_t = fields::LocalField<BGptGluon<bgf_t, 0, DIM>, DIM>;
    geometry::Geometry<DIM>::extents_t extents;
    std::fill(extents.begin(), extents.end(), L);
    field_t field{extents,1,0, {}};
    pt::Direction<DIM> direction;
    geometry::Geometry<DIM>::raw_pt_t coords{0,0,0,0};

    auto point = field.mk_point(coords);

    // std::iota(U.begin(), U.end(),1);

    kernels::StapleSqKernel<field_t> kernel{direction};
    kernel(field, point);
    auto value = kernel.val;
    std::cout << value << "\n";
}