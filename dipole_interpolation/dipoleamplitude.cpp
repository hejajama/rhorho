
#include "dipoleamplitude.hpp"
#include "interpolation2d.hpp"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

DipoleInterpolator2D* InterpolatorFromFile(string fname)
{
    // Note, we want to initailize a 2d interpolator
    // It needs xgrid and ygrid containing ONLY the
    // distinct grid points
    // z[x_i,y_j] = zgrid[j*xgridsize+i]
    
    // So we say that b is y coordinate, r is x coord
    
    ifstream file(fname.c_str());
    
    bool initialized = false;
    bool r_grid_ready = false;
    vector<double> bvals;
    vector<double> rvals;
    vector<double> datavals;
    
    while(!file.eof())
    {
        string line;
        getline(file, line);
        if (line.substr(0,1)=="#")
            continue;
        if (line.length() < 3)
            continue;
        stringstream ss(line);
        double b,r,d;
        ss >> b;
        ss >> r;
        ss >> d;
        
        if (!initialized) // First line
        {
            bvals.push_back(b);
            rvals.push_back(r);
            initialized=true;
        }
        else
        {
            // Did we find new b value
            if (std::abs(b - bvals[bvals.size()-1]) > 1e-5)
            {
                bvals.push_back(b);
                r_grid_ready=true;
            }
            
            if (r_grid_ready != true)
            {
                rvals.push_back(r);
            }
            
        }
        datavals.push_back(d);
        
    }
    
    cout << "# Interpolation constructed from the data file " << fname  << endl;
    cout <<  "# In total " << datavals.size() << " points (" << rvals.size() << " rpoints, " << bvals.size() << " bpoints)" << endl;
	cout << "# Minr: " << rvals[0] << ", maxr: " << rvals[rvals.size()-1] << ", maxb: " << bvals[bvals.size()-1] << " GeV^-1" << endl;
    
    return  new DipoleInterpolator2D(rvals, bvals, datavals);
    
    
    
}

LCPT_Dipole::LCPT_Dipole(string fname, string v2fname)
{
    // Datafile structure is
    // b r N(r,b)
    interpolator2d = InterpolatorFromFile(fname);
    angledep=false;
    if (v2fname != "")
    {
        angledep=true;
        v2_interpolator2d = InterpolatorFromFile(v2fname);
    }
        
    
    minr=interpolator2d->MinX();
    maxr=interpolator2d->MaxX();
    minb=interpolator2d->MinY();
    maxb=interpolator2d->MaxY();
    
    
    out_of_range_warnings=true;
    
    
    
}

double LCPT_Dipole::Evaluate(double r, double b)
{
    // Too small r or too large b: assume N(r,b)=0
    if (r < minr or b>maxb)
        return 0;
    
    // Other out of range errors
    if (r > maxr or b < minb )
    {
        if (out_of_range_warnings)
            cerr << "r,b=" << r <<", " << b << " is out of range! returning 0..." << endl;
        return 0;
    }
    
    return interpolator2d->Evaluate(r,b);
    
}

double LCPT_Dipole::v2(double r, double b)
{
    if (angledep == false)
        return 0;
    
    // Too small r or too large b: assume N(r,b)=0
    if (r < minr or b>maxb)
        return 0;
    
    // Other out of range errors
    if (r > maxr or b < minb )
        return 0;
    
    return v2_interpolator2d->Evaluate(r,b);
}

double  LCPT_Dipole::Evaluate(double r, double b, double phirb)
{
    if (angledep == false)
        return Evaluate(r,b);
    
    return Evaluate(r,b)*(1.0 + 2.0*v2(r,b)*std::cos(2.0*phirb));
    
}
LCPT_Dipole::~LCPT_Dipole()
{
    delete interpolator2d;
    if (angledep)
        delete v2_interpolator2d;
}
