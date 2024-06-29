#include "spectrum-raw.h"

void SpectrumRaw::append(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten)
{
    mzs.push_back(mz);
    intensities.push_back(inten);
}
void SpectrumRaw::append(const SpectrumRaw &sp)
{
    mzs.insert(mzs.end(), sp.getMZs().begin(), sp.getMZs().end());
    intensities.insert(intensities.end(), sp.getIntensities().begin(), sp.getIntensities().end());
}
