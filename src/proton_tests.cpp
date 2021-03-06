
#include "proton.hpp"

//#define BOOST_TEST_MODULE Helpers
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(ProtonTest)
BOOST_AUTO_TEST_CASE(ProtonWF_uunorm)
{
    Proton p;
    p.SetM(0.26);
    p.SetBeta(0.55);
    Vec k1(0.3,0);
    Vec k2(0.4,0);
    double x1=0.2; double x2=0.4;
    p.SetWFNormalizationCoefficient(1);
    BOOST_CHECK_CLOSE(p.WaveFunction(k1, k2, x1, x2), 0.00189579, 0.01);
}
BOOST_AUTO_TEST_CASE(ProtonWFNorm)
{
    Proton p;
    p.SetM(0.26);
    p.SetBeta(0.55);
    // Numerical MC integral, can't require too good precision...
    BOOST_CHECK_CLOSE(p.ComputeWFNormalizationCoefficient(), 43389.7232, 0.1);
}

BOOST_AUTO_TEST_CASE(ProtonPowerLawWFNorm)
{
    Proton p;
    p.SetWaveFunction(Power);
    p.SetM(0.26);
    p.SetBeta(0.55);
    p.SetP(3.5);
    // Numerical MC integral, can't require too good precision...
    BOOST_CHECK_CLOSE(p.ComputeWFNormalizationCoefficient(), 1.50206095e6, 1);
}

BOOST_AUTO_TEST_SUITE_END()
