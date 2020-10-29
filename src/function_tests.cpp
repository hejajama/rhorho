#include "functions.hpp"

//#define BOOST_TEST_MODULE Helpers
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(HelperFunctions)
BOOST_AUTO_TEST_CASE(B0tests)
{
    BOOST_CHECK_CLOSE(B0(1, 2), 0.462098, 0.0001f);
    BOOST_CHECK_CLOSE(B0(10, 50), 0.0193589, 0.0001f);
    BOOST_CHECK_CLOSE(B0(0.01,0.001), 418.839640629419, 0.0001f);
}

BOOST_AUTO_TEST_CASE(Ftest)
{
    Vec l(2,1);
    Vec l1(-1,3);
    double alpha=0.1; double m2=0.1*0.1;
    BOOST_CHECK_CLOSE(F_int_B0(l, l1, alpha, m2), -0.33946139702, 0.0001);
    
    l1=Vec(0,0,0);
    l=Vec(1,0);
    alpha=0.01;
    m2=0.05*0.05;
    BOOST_CHECK_CLOSE(F_int_B0(l, l1, alpha, m2), -0.584170078, 0.0001);
}
BOOST_AUTO_TEST_SUITE_END()
