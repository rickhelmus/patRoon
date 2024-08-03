#include "spectrum-raw.h"

void SpectrumRaw::append(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten)
{
    mzs.push_back(mz);
    intensities.push_back(inten);
}

void SpectrumRaw::append(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten, SpectrumRawTypes::Mobility mob)
{
    append(mz, inten);
    mobilities.push_back(mob);
}

void SpectrumRaw::append(const SpectrumRaw &sp)
{
    mzs.insert(mzs.end(), sp.getMZs().begin(), sp.getMZs().end());
    intensities.insert(intensities.end(), sp.getIntensities().begin(), sp.getIntensities().end());
    mobilities.insert(mobilities.end(), sp.getMobilities().begin(), sp.getMobilities().end());
}

void SpectrumRaw::insert(size_t i, SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten)
{
    mzs.insert(mzs.begin() + i, mz);
    intensities.insert(intensities.begin() + i, inten);
}

void SpectrumRaw::insert(size_t i, SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten,
                         SpectrumRawTypes::Mobility mob)
{
    mzs.insert(mzs.begin() + i, mz);
    intensities.insert(intensities.begin() + i, inten);
    mobilities.insert(mobilities.begin() + i, mob);
}
