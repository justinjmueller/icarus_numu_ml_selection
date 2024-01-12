/**
 * @file container.h
 * @brief Header file defining a container for CAFAna Spectrum
 * @author justin.mueller@colostate.edu
*/
#ifndef CONTAINER_H
#define CONTAINER_H

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

/**
 * Container class for CAFAna Spectrum objects. Allows for easier
 * configuration of a set of CAFAna Spectrum, and handles the output of
 * the resulting histograms to a ROOT file.
*/
struct SpecContainer
{
    ana::SpectrumLoader loader;
    std::vector<const char *> names;
    std::vector<ana::Spectrum*> spectra;
    TFile output_file;
    float override_pot;
    float target_pot;

    /**
     * Constructor for SpecContainer.
     * @param in_name is the name of the input CAF file.
     * @param out_name is the name of the output ROOT file.
    */
    SpecContainer(const char * in_name, const char * out_name, float opot=-1, float tpot=-1)
    : loader(in_name),
      output_file(out_name, "recreate"),
      override_pot(opot),
      target_pot(tpot) { }
    
    /**
     * Adds a new CAFAna Spectrum (1D) object to the container.
     * @param n is the name of the spectrum.
     * @param b is the Binning of the spectrum.
     * @param v is the variable defining the spectrum
     * @return none.
    */
    void add_spectrum1d(const char * n, const ana::Binning b, const ana::SpillMultiVar v)
    {
        names.push_back(n);
        spectra.push_back(new ana::Spectrum(n, b, loader, v, ana::kNoSpillCut));
        if(override_pot != -1) spectra.back()->OverridePOT(override_pot);
    }

    /**
     * Adds a new CAFAna Spectrum (2D) object to the container.
     * @param n is the name of the spectrum.
     * @param b0 is the first set of Binnings.
     * @param b1 is the second set of Binnings.
     * @param v0 is the first variable.
     * @param v1 is the second variable.
     * @return none.
    */
    void add_spectrum2d(const char * n, const ana::Binning b0, const ana::Binning b1,
                        const ana::SpillMultiVar v0, const ana::SpillMultiVar v1)
    {
        names.push_back(n);
        spectra.push_back(new ana::Spectrum(n, loader, b0, v0, b1, v1, ana::kNoSpillCut));
        if(override_pot != -1) spectra.back()->OverridePOT(override_pot);
    }

    /**
     * Runs the selection and fills each spectrum in the container.
     * @return none.
    */
    void run()
    {
        loader.Go();
        for(size_t i(0); i < spectra.size(); ++i)
            output_file.WriteObject(spectra[i]->ToTHX(target_pot != -1 ? target_pot : 1), names[i]);
        output_file.Close();
    }

    /**
     * Destructor for the container. Releases dynamically allocated memory
     * for the CAFAna Spectrum objects.
    */
    ~SpecContainer()
    {
        for(ana::Spectrum * s : spectra)
            delete s;
    }
};
#endif