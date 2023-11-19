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

    /**
     * Constructor for SpecContainer.
     * @param in_name is the name of the input CAF file.
     * @param out_name is the name of the output ROOT file.
    */
    SpecContainer(const char * in_name, const char * out_name)
    : loader(in_name),
      output_file(out_name, "recreate") { }
    
    /**
     * Adds a new CAFAna Spectrum object to the container.
     * @param n is the name of the spectrum.
     * @param b is the Binning of the spectrum.
     * @param v is the variable defining the spectrum
     * @return none.
    */
    void add_spectrum(const char * n, const ana::Binning b, const ana::SpillMultiVar v)
    {
        names.push_back(n);
        spectra.push_back(new ana::Spectrum(n, b, loader, v, ana::kNoSpillCut));
    }

    /**
     * Runs the selection and fills each spectrum in the container.
     * @return none.
    */
    void run()
    {
        loader.Go();
        for(size_t i(0); i < spectra.size(); ++i)
            output_file.WriteObject(spectra[i]->ToTH1(1), names[i]);
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