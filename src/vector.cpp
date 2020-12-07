/*
 * Simple class for 2D/3D vectors
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "vector.hpp"
#include <cmath>
#include <cstdlib>

inline double SQR(double x){ return x*x; }

// **********
// Vec Class

Vec::Vec() { x=0; y=0; z=0; }
Vec::Vec(REAL x_, REAL y_) { x=x_; y=y_; z=0; }
Vec::Vec(REAL x_, REAL y_, REAL z_) { x=x_; y=y_, z=z_;}
Vec::Vec(const Vec& v) { x=v.GetX(); y=v.GetY(); z=v.GetZ(); }

void Vec::SetX(REAL x_) { x=x_; }
void Vec::SetY(REAL y_) { y=y_; }
void Vec::SetZ(REAL z_) { z=z_; }

REAL Vec::GetX() const { return x; }
REAL Vec::GetY() const { return y; }
REAL Vec::GetZ() const { return z; }

Vec& Vec::operator+=(Vec& v)
{
    x+=v.GetX();
    y+=v.GetY();
    z+=v.GetZ();
    
    return *this;
}

Vec& Vec::operator-=(Vec& v)
{
    x-=v.GetX();
    y-=v.GetY();
    z-=v.GetZ();
    
    return *this;
} 

Vec& Vec::operator=(const Vec& v)
{
    x=v.x;
    y=v.y;
    z=v.z;
    return *this;
} 

Vec  Vec::operator+(const Vec& v)
{
    Vec sum;
    sum.SetX(x+v.GetX());
    sum.SetY(y+v.GetY());
    sum.SetZ(z+v.GetZ());
    return sum; 
}

Vec Vec::operator-(const Vec& v)
{
    Vec sum;
    sum.SetX(x-v.GetX());
    sum.SetY(y-v.GetY());
    sum.SetZ(z-v.GetZ());
    return sum; 
}

Vec& Vec::operator*=(REAL c)
{
    x*=c; y*=c; z*=c;
    
    return *this;
}

double Vec::operator*(Vec& v)
{
    return x*v.GetX() + y*v.GetY() + z*v.GetZ();
}

Vec Vec::operator*(REAL c)
{
    Vec tmp(x*c,y*c,z*c);
    return tmp;
}

REAL Vec::LenSqr()
{
    return SQR(x)+SQR(y)+SQR(z);
}

REAL Vec::Len()
{
    return sqrt(LenSqr());
}

std::ostream& operator<<(std::ostream& os, Vec& ic)
{
    return os << "(" << ic.GetX() << ", " << ic.GetY() <<
        ", " << ic.GetZ() << "), |vec| = " << ic.Len() << " ";
}

void Vec::Rotate2D(double angle)
{
    // Rotate counterclokwise by angle [rad]
    double s = std::sin(angle);
    double c = std::cos(angle);
    
    double oldx = x;
    double oldy = y;
    
    x = c*oldx + s*oldy;
    y = -s*oldx + c*oldy;
    
}
