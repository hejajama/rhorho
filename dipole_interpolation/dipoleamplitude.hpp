
/*
 * Dipole amplitude from LCPT calculation
 *
 * H. Mäntysaari, R. Paatelainen, A. Dumitru, 2021
 * arXiv:2103.11682 [hep-ph]
 */

#ifndef LCPT_DIPOLE_H_
#define LCPT_DIPOLE_H_

#include "interpolation2d.hpp"
#include <string>

class LCPT_Dipole
{
public:
    LCPT_Dipole(std::string avgdipolefile, std::string v2file="");
    
    // Evaluate dipole in untis [1/GeV, 1/GeV]
    double Evaluate(double r, double b);
    double Evaluate(double r, double b, double phi_r_b);
    
    double v2(double r, double b);
    
    ~LCPT_Dipole();
    void Set_out_of_range_warnings(bool s) { out_of_range_warnings = s;}
    
    
private:
    DipoleInterpolator2D* interpolator2d;
    DipoleInterpolator2D* v2_interpolator2d;
    double minr;
    double maxr;
    double minb;
    double maxb;
    bool out_of_range_warnings;
    bool angledep;
};


#endif
