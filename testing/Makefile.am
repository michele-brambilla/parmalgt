# global includes

AM_CPPFLAGS = -I${GTEST_DIR}/include -I${GTEST_DIR} -I..\
	-I${top_srcdir} -DHAVE_TR1 -DDSFMT_MEXP=19937 

# gtest stuff

lib_LTLIBRARIES = lib/libgtest.la lib/libgtest_main.la
lib_libgtest_la_SOURCES = ${GTEST_DIR}/src/gtest-all.cc
lib_libgtest_main_la_SOURCES = ${GTEST_DIR}/src/gtest_main.cc
lib_libgtest_main_la_LIBADD = lib/libgtest.la

# test programs we want to build

# NOTE: The point test does not work with the new lattice class at the
# moment

bin_PROGRAMS = test_Geometry test_Bgfield test_Helper test_SU3 test_newQCDpt \
	test_PtTypes #test_Point

test_Geometry_SOURCES = test_Geometry.cc
test_Geometry_LDADD = lib/libgtest_main.la -lpthread


test_Bgfield_SOURCES = test_Bgfield.cc ../MyRand.cc \
		       ../MyMath.cc ../Background.cc
test_Bgfield_LDADD = lib/libgtest_main.la -lpthread

test_Helper_SOURCES = test_Helper.cc ../MyRand.cc
test_Helper_LDADD = lib/libgtest_main.la -lpthread

test_SU3_SOURCES = test_SU3.cc ../MyMath.cc ../MyRand.cc
test_SU3_LDADD = lib/libgtest_main.la -lpthread

test_newQCDpt_SOURCES = test_newQCDpt.cc ../MyMath.cc ../MyRand.cc \
	../newQCDpt.cc ../Background.cc
test_newQCDpt_LDADD = lib/libgtest.la -lpthread

test_PtTypes_SOURCES = test_PtTypes.cc ../MyMath.cc ../MyRand.cc
test_PtTypes_LDADD = lib/libgtest_main.la -lpthread

TESTS = ${bin_PROGRAMS}