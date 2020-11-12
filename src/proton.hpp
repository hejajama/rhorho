#ifndef _PROTON_H
#define _PROTON_H

#include "vector.hpp"


enum ProtonWaveFunction
{
    HarmoinicOscillator,
};

class Proton
{
public:
    Proton();
    double WaveFunction(Vec k1, Vec k2, double x1, double x2);
    double ComputeWFNormalizationCoefficient();
    
    void SetBeta(double b) { beta=b; }
    void SetM(double m) { mq=m; }
    void SetWFNormalizationCoefficient(double n) { wf_normalization=n; }
    double GetM() { return mq; }
    double GetBeta() { return beta; }
private:
    ProtonWaveFunction wave_function;
    double beta;   // In HO wave function
    double mq;      // Quark mass can be different from other parts of the code
    double wf_normalization;   // normalization coefficient
    
   
};

#endif

