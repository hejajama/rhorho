
#include <iostream>
#include <string>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <cmath>
#include "vector.hpp"
#include "gitsha1.h"

#include "dipoleamplitude.hpp"
using namespace std;

double StrToReal(std::string str)
{
    std::stringstream buff(str);
    double tmp;
    buff >> tmp;
    return tmp;
}
long int StrToInt(std::string str)
{
    std::stringstream buff(str);
    double  tmp;
    buff >> tmp;
    return static_cast<long int>(tmp);
}


void handler (const char * reason,
const char * file,
int line,
int gsl_errno)
{
    cerr << "# Error " << gsl_errno << " line " << line << " file " << file << " reason " << reason << endl;

}

int main(int argc, char* argv[])
{
    gsl_set_error_handler(handler);

    cout << "# Command: ";
    for (int i=1; i<argc; i++)
        cout << argv[i] << " ";
    cout << endl;
    cout << "# Git version " << g_GIT_SHA1 << " local repo " << g_GIT_LOCAL_CHANGES << " main build " << __DATE__  << " " << __TIME__ << endl;


    LCPT_Dipole dipole("/Users/hejajama/Nextcloud/projects/rhorho/dipole_2d_data/x_0.01/fixed_ir_nlo_mc_5e7_mq_0.2_as_0.2_large.dat");
    dipole.Set_out_of_range_warnings(false);
    
    cout << "N(r=1,b=2)=" << dipole.Evaluate(1,2) << endl;
    cout << "N(r=1.2,b=4.1)=" << dipole.Evaluate(1.2,4.1) << endl;
    cout << "N(r=7,b=2)=" << dipole.Evaluate(7,4.1) << endl;
}
