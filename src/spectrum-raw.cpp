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

void SpectrumRaw::sort(SpectrumRawTypes::MSSortType stype)
{
    if (empty() || stype == SpectrumRawTypes::MSSortType::NONE)
        return;
    
    const auto mobMZLT = [this](size_t i, size_t j)
    {
        return std::tie(mobilities[i], mzs[i]) < std::tie(mobilities[j], mzs[j]);
    };
    const auto sortedInds = (stype == SpectrumRawTypes::MSSortType::MZ) ? getSortedInds(mzs) : getSortedInds(mzs, mobMZLT);
    
    std::vector<SpectrumRawTypes::Mass> sortedMZs(mzs.size());
    std::vector<SpectrumRawTypes::Intensity> sortedInts(intensities.size());
    for (size_t i=0; i<sortedInds.size(); ++i)
    {
        const auto j = sortedInds[i];
        sortedMZs[i] = mzs[j];
        sortedInts[i] = intensities[j];
    }
    mzs = std::move(sortedMZs);
    intensities = std::move(sortedInts);
    
    if (hasMobilities())
    {
        std::vector<SpectrumRawTypes::Mobility> sortedMobs(mobilities.size());
        for (size_t i=0; i<sortedInds.size(); ++i)
            sortedMobs[i] = mobilities[sortedInds[i]];
        mobilities = std::move(sortedMobs);
    }
}

void SpectrumRawAveraged::append(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten,
                                 SpectrumRawTypes::PeakAbundance abr, SpectrumRawTypes::PeakAbundance aba)
{
    SpectrumRaw::append(mz, inten);
    abundancesRel.push_back(abr); abundancesAbs.push_back(aba);
}

void SpectrumRawAveraged::append(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten,
                                 SpectrumRawTypes::PeakAbundance abr, SpectrumRawTypes::PeakAbundance aba,
                                 SpectrumRawTypes::PeakAbundance avgPrvAbR, SpectrumRawTypes::PeakAbundance avgPrvAbA)
{
    append(mz, inten, abr, aba);
    averagedPrevAbundancesRel.push_back(avgPrvAbR); averagedPrevAbundancesAbs.push_back(avgPrvAbA);
}

void SpectrumRawAveraged::append(const SpectrumRawAveraged &sp)
{
    SpectrumRaw::append(sp);
    abundancesRel.insert(abundancesRel.end(), sp.getAbundancesRel().begin(), sp.getAbundancesRel().end());
    abundancesAbs.insert(abundancesAbs.end(), sp.getAbundancesAbs().begin(), sp.getAbundancesAbs().end());
}

void SpectrumRawAveraged::setAvgPrevAbundance(size_t i, SpectrumRawTypes::PeakAbundance abr,
                                              SpectrumRawTypes::PeakAbundance aba)
{
    averagedPrevAbundancesRel[i] = abr; averagedPrevAbundancesAbs[i] = aba;
}
