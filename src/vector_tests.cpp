#include "vector.hpp"

//#define BOOST_TEST_MODULE Vector
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(VectorTest)
BOOST_AUTO_TEST_CASE(VecTest)
{
    Vec v1(1,2,3);
    Vec v2(3,4,5);
    Vec v3 = v1+v2;
    BOOST_CHECK_CLOSE(v3.GetX(), 4, 0.0001);
    BOOST_CHECK_CLOSE(v3.GetY(), 6, 0.0001);
    BOOST_CHECK_CLOSE((v2*2).GetZ(), 10, 0.0001);
    
    // Dot product
    BOOST_CHECK_CLOSE(v1*v2, 26, 0.0001);
    
    Vec v4 = v1+v2*2-v3*3;
    BOOST_CHECK_CLOSE(v4.GetY(), -8, 0.0001);
    
    Vec q1(1,2);
    Vec q2(3,4);
    Vec K = q1*(-1)-q2;
    BOOST_CHECK_CLOSE(K.GetY(), -6, 0.0001);
    BOOST_CHECK_CLOSE(K.GetX(), -4, 0.0001);
}

BOOST_AUTO_TEST_SUITE_END()
