/**
 * @file reweight.h
 * @brief Header file defining the calculation of systematics implemented as
 * reweightable interactions.
 * @author justin.mueller@colostate.edu
*/

#ifndef REWEIGHT_H
#define REWEIGHT_H

#include <string>
#include <map>
#include <vector>
#include <chrono>
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH2D.h"
#include "TH1D.h"
#include "types.h"
#include "vars.h"
#include "utilities.h"

/**
 * Calculates the TH2D histograms for the reconstructed quantities and the
 * systematic universe weights. The histograms are stored in a map that maps a
 * string (systematic name) to a TH2D* object (X = reconstructed quantity, Y =
 * systematic universe). The systematic parameters are specified in vars.h.
 * @param input_file_name The name of the input file (TFile).
 * @param reco_map The map that stores the selected interactions.
 * @param weights The map that will store the TH2D histograms.
 * @return the POT of the input file.
*/
double calc_reweight_systematics(std::string input_file_name, std::map<index_t, std::vector<double>> & reco_map, weights_t & weights)
{
    /**
     * Open the input file (TFile) and attach a TTreeReader to the "recTree".
     * The reader serves as an interface to the TTree, allowing us to access the
     * event metadata and the true interactions (weights). It also saves some
     * overhead on the I/O operations as not all event information is loaded.
     * TTreeReaderValues and TTreeReaderArrays are used to access the event
     * metadata and the true interactions, respectively.
    */
    TFile * file = new TFile(input_file_name.c_str(), "READ");
    if(file->IsZombie() || !file->GetListOfKeys()->Contains("recTree"))
    {
        std::cerr << "Error: File " << input_file_name << " does not exist." << std::endl;
        return 0;
    }
    TTreeReader reader("recTree", file);
    TTreeReaderValue<uint32_t> run(reader, "rec.hdr.run");
    TTreeReaderValue<uint32_t> subrun(reader, "rec.hdr.subrun");
    TTreeReaderValue<uint32_t> evt(reader, "rec.hdr.evt");
    TTreeReaderArray<caf::SRTrueInteraction> mc(reader, "rec.mc.nu");

    /**
     * Begin main loop over the events in the TTree. For each event, we will
     * loop over the true interactions (possibly more than one neutrino per
     * event) and check that the interaction has indeed been selected. If so,
     * we will loop over the systematic parameters and then over the
     * reconstructed quantities. For each reconstructed quantity, we will fill
     * the corresponding TH2D histogram with the reconstructed value and the
     * systematic universe weight.
    */
    while(reader.Next())
    {
        // Loop over the true interactions (neutrinos) in the event.
        for(const caf::SRTrueInteraction & nu : mc)
        {
            // Check if the interaction has been selected.
            index_t index(*run, *subrun, *evt, nu.index);
            if(reco_map.find(index) != reco_map.end())
            {
                // Loop over the systematic parameters.
                for(std::pair<std::string, size_t> syst : systs)
                {
                    // Loop over the reconstructed quantities.
                    for(size_t ri(0); ri < reco_vars.size(); ++ri)
                    {
                        RecoVar & r = reco_vars[ri];
                        std::string name = syst.first + "_" + r.name;
                        // Check if the TH2D histogram has been created.
                        if(weights.find(name) == weights.end())
                        {
                            size_t nuniv = nu.wgt[syst.second].univ.size();
                            weights.insert(std::make_pair(name, new TH2D(name.c_str(), name.c_str(), r.nbins, r.xmin, r.xmax, nuniv, 0, nuniv)));
                            weights[name]->SetDirectory(nullptr);
                            weights.insert(std::make_pair(name+"_cv", new TH1D((name+"_cv").c_str(), (name+"_cv").c_str(), r.nbins, r.xmin, r.xmax)));
                            weights[name+"_cv"]->SetDirectory(nullptr);
                        }
                        // Fill the TH2D histogram with the reconstructed value and the systematic universe weight.
                        for(size_t i(0); i < nu.wgt[syst.second].univ.size(); ++i)
                            static_cast<TH2D*>(weights[name])->Fill(reco_map[index][ri], i, nu.wgt[syst.second].univ[i]);
                        // Fill the TH1D histogram with the reconstructed value.
                        weights[name+"_cv"]->Fill(reco_map[index][ri]);
                    } // End loop over the reconstructed quantities.
                } // End loop over the systematic parameters.
            } // End check if the interaction has been selected.
        } // End loop over the true interactions.
    } // End main loop over the events.

    /**
     * For each systematic and each reconstructed quantity, we want to use
     * the created histograms to calculate the covariance matrix across the
     * bins of the reconstructed quantity.
    */
    /*
    for(std::pair<std::string, size_t> syst : systs)
    {
        for(const RecoVar & r : reco_vars)
        {
            std::string name = syst.first + "_" + r.name;
            calc_covariance(weights, name);
        }
    }*/

    double POT = (static_cast<TH1D*>(file->Get("TotalPOT")))->GetArray()[1];

    // Close the file and release the allocated memory.
    file->Close();
    delete file;

    return POT;
}

#endif