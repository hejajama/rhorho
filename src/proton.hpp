#ifndef _PROTON_H
#define _PROTON_H

#include "vector.hpp"
#include <string>

enum ProtonWaveFunction
{
    HarmoinicOscillator, // (2) in hep-ph/9402214
    Power //(3) in hep-ph/9402214
};

std::string WaveFunctionString(ProtonWaveFunction wf);

class Proton
{
public:
    Proton();
    double WaveFunction(Vec k1, Vec k2, double x1, double x2);
    double ComputeWFNormalizationCoefficient();
    
    void SetBeta(double b) { beta=b; }
    void SetM(double m) { mq=m; }
    void SetWFNormalizationCoefficient(double n) { wf_normalization=n; }
    void SetP(double p) { wf_power = p; }
    double GetM() { return mq; }
    double GetBeta() { return beta; }
    double GetP() { return wf_power; }
    void SetWaveFunction(ProtonWaveFunction wf) { wave_function = wf; }
    ProtonWaveFunction GetWaveFunction() { return wave_function; }
private:
    ProtonWaveFunction wave_function;
    double beta;   // In HO wave function
    double mq;      // Quark mass can be different from other parts of the code
    double wf_normalization;   // normalization coefficient
    double wf_power;    // Exponent in the powerlaw wave function
    
   
};

#endif

