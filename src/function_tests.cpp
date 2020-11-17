#include "functions.hpp"
#include "diagram_integrator.hpp"
//#define BOOST_TEST_MODULE Helpers
#include <boost/test/unit_test.hpp>

double inthelperf_mc_lo(double *vec, size_t dim, void* p);

BOOST_AUTO_TEST_SUITE(HelperFunctions)
BOOST_AUTO_TEST_CASE(B0tests)
{
    BOOST_CHECK_CLOSE(B0(1, 2), 0.462098, 0.0001f);
    BOOST_CHECK_CLOSE(B0(10, 50), 0.0193589, 0.0001f);
    BOOST_CHECK_CLOSE(B0(0.01,0.001), 418.839640629419, 0.0001f);
}

BOOST_AUTO_TEST_CASE(Ftest)
{
    F_worker F(20,0.00001);
    
    Vec l(2,1);
    Vec l1(-1,3);
    double alpha=0.1; double m2=0.1*0.1;
    BOOST_CHECK_CLOSE(F.F_int_B0(l, l1, alpha, m2), -0.33946139702, 0.0001);
    
    l1=Vec(0,0,0);
    l=Vec(1,0);
    alpha=0.01;
    m2=0.05*0.05;
    BOOST_CHECK_CLOSE(F.F_int_B0(l, l1, alpha, m2), -0.584170078, 0.0001);
}


BOOST_AUTO_TEST_CASE(LO_DIAG)
{
    double x1=0.3; double x2=0.4;
    double intvec[6] = {0.1, 0, 0.2,0, x1,x2}; // k1x  k1y k2x k2y x1 x2
    
    DiagramIntegrator *integrator = new DiagramIntegrator;
    integrator->GetProton().SetBeta(0.55);
    integrator->GetProton().SetM(0.26);
    integrator->GetProton().ComputeWFNormalizationCoefficient();
    
    Diagram diag = integrator->DiagramType("LO");
    //integrator->SetMCIntPoints(1e5);
    inthelper_diagint helper;
    helper.q1=Vec(1,0); helper.q2=Vec(0.5,0);
    helper.integrator=integrator;
    helper.diag = diag;
    
    double integrand =inthelperf_mc_lo(intvec, 6, &helper);
    
    // Jacobian and normalizatoin factors should be added to the mathematica result
    double mathematica =-850030;
    mathematica *= 16.*std::pow(M_PI,3.)  / (8.0*x1*x2*(1.-x1-x2)*std::pow(2.0*M_PI,6.0));
    
    BOOST_CHECK_CLOSE(integrand, mathematica, 0.001);
    
    
    
    delete integrator;
}
BOOST_AUTO_TEST_SUITE_END()
