/**
 * @file utilities.h
 * @brief This file contains various utility functions for the systematics code.
 * @author justin.mueller@colostate.edu
*/

#ifndef UTILITIES_H
#define UTILITIES_H

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH2D.h"
#include "TH1D.h"
#include "vars.h"

/**
 * Read the TTree containing metadata of all events from the input file and
 * store them in a vector of meta_t objects. The meta_t object contains the
 * run, subrun, and event number of all events in the sample.
 * @param events The vector to store the event metadata.
 * @param file_name The name of the input file.
*/
void read_event_metadata(std::vector<meta_t> & events, const std::string & file_name)
{
    /**
     * Attach to the input file and check that it has been opened successfully.
    */
    TFile * file = new TFile(file_name.c_str(), "READ");
    if(!file->IsOpen())
    {
        std::cerr << "Error: selected file not found" << std::endl;
        return;
    }

    /**
     * Attach a TTreeReader to the "selected" TTree and create TTreeReaderValues
     * for the event metadata.
    */
    TTreeReader reader("events", file);
    TTreeReaderValue<double> run(reader, "run");
    TTreeReaderValue<double> subrun(reader, "subrun");
    TTreeReaderValue<double> event(reader, "event");

    /**
     * Loop over the events and store the event metadata.
    */
    while(reader.Next())
    {
        meta_t meta(*run, *subrun, *event);
        events.push_back(meta);
    }
    file->Close();
}

/**
 * Read the TTree containing selected events from the input file and store the
 * reconstructed quantities in a map. The map is indexed by an index_t object
 * (run, subrun, event, nu_id) and contains a vector of reconstructed quantities
 * (as configured by the reco_vars object in vars.h).
 * @param reco_map The map to store the reconstructed quantities.
 * @param file_name The name of the input file.
*/
void read_selected(std::map<index_t, std::vector<double>> & reco_map, const std::string & file_name)
{
    /**
     * Attach to the input file and check that it has been opened successfully.
    */
    TFile * file = new TFile(file_name.c_str(), "READ");
    if(!file->IsOpen())
    {
        std::cerr << "Error: selected file not found" << std::endl;
        return;
    }

    /**
     * Attach a TTreeReader to the "selected" TTree and create TTreeReaderValues
     * for the event metadata and the reconstructed quantities (as configured
     * in vars.h).
    */
    TTreeReader reader("selected_1mu1p", file);
    TTreeReaderValue<double> run(reader, "run");
    TTreeReaderValue<double> subrun(reader, "subrun");
    TTreeReaderValue<double> event(reader, "event");
    TTreeReaderValue<double> nu_id(reader, "nu_id");
    std::vector<TTreeReaderValue<double>> vars;
    for(size_t ri(0); ri < reco_vars.size(); ++ri)
        vars.push_back(TTreeReaderValue<double>(reader, reco_vars[ri].name.c_str()));

    /**
     * Loop over the selected interactions and store the reconstructed quantities.
    */
    while(reader.Next())
    {
        index_t index(*run, *subrun, *event, *nu_id);
        reco_map.insert(std::make_pair(index, std::vector<double>()));
        for(size_t ri(0); ri < reco_vars.size(); ++ri)
            reco_map[index].push_back(*vars[ri]);
    }
    file->Close();
}

/**
 * Calculate the covariance matrix for a given TH2D histogram across the
 * reconstructed quantity bins (X-axis). The covariance matrix is stored in a
 * TH2D object with the same binning (in X) as the input histogram.
*/
void calc_covariance(weights_t weights, std::string systname)
{
    /**
     * Retrieve the TH2D histogram for the given systematic parameter
     * and the corresponding central value (TH1D).
    */
    TH2D * hist = static_cast<TH2D*>(weights[systname]);
    TH1D * central = static_cast<TH1D*>(weights[systname + "_cv"]);

    /**
     * Create the covariance matrix TH2D object with the same binning as the
     * input histogram. The name of the covariance matrix follows the pattern
     * <hist_name>_cov.
    */
    std::string name = std::string(hist->GetName()) + "_cov";
    TH2D * cov = new TH2D(name.c_str(), name.c_str(), hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    /**
     * Loop over the reconstructed quantity bins and calculate the covariance
     * between each pair of bins.
    */
    for(int xi(1); xi <= hist->GetNbinsX(); ++xi)
    {
        for(int xj(1); xj <= hist->GetNbinsX(); ++xj)
        {
            double cov_ij(0);
            for(int yi(1); yi <= hist->GetNbinsY(); ++yi)
            {
                for(int yj(1); yj <= hist->GetNbinsY(); ++yj)
                    cov_ij += (hist->GetBinContent(xi, yi) - central->GetBinContent(xi)) * (hist->GetBinContent(xj, yj) - central->GetBinContent(xj));
            } // End loop over the universes.
            cov_ij /= hist->GetNbinsY();
            cov->SetBinContent(xi, xj, cov_ij);
        } // End (inner) loop over the bins.
    } // End (outer) loop over the bins.
}

#endif
