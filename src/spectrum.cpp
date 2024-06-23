#include "spectrum.h"

void Spectrum::append(double mz, double inten)
{
    mzs.push_back(mz);
    intensities.push_back(inten);
}
void Spectrum::append(const Spectrum &sp)
{
    mzs.insert(mzs.end(), sp.getMZs().begin(), sp.getMZs().end());
    intensities.insert(intensities.end(), sp.getIntensities().begin(), sp.getIntensities().end());
}
