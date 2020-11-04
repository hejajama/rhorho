#ifndef _FUNCTIONS_HPP
#define _FUNCTIONS_HPP

#include "vector.hpp"
#include <gsl/gsl_integration.h>

double B0(double hsqr, double Delta);

/*
 * Class to handle evaluation of F(l,l1,alpha,m)
 * Note that THIS IS NOT THREAD SAFE!
 */
class F_worker
{
public:
    F_worker(const int int_divisions_, const double intaccuracy=0.001);
    double F_int_B0(Vec l, Vec l1, double alpha, double msqr);
    ~F_worker();
   
private:
    gsl_integration_workspace *ws;
    double int_accuracy;
    int int_divisions;
};



#endif
