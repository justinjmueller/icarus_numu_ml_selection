/**
 * @file variation.h
 * @brief Header file defining the calculation of systematics implemented as variation
 * samples.
 * @author justin.mueller@colostate.edu
*/

#ifndef VARIATION_H
#define VARIATION_H

#include <random>
#include <string>
#include <map>
#include <vector>
#include "TH2D.h"
#include "TH1D.h"
#include "types.h"
#include "vars.h"
#include "utilities.h"

/**
 * Calculates the histograms associated with the variation systematics. A
 * bootstrapping method is used to better estimate the bin-to-bin correlations
 * in the effect of the systematic on the reconstructed quantities. The
 * 
*/
void calc_variation_systematics(std::string systname, std::string nominal, std::string variation, weights_t & weights)
{
    /**
     * Read the event metadata for both the nominal and variation samples.
     * We need this information to ensure that we use only the intersection
     * of events from each sample.
    */
    // Nominal sample.
    std::vector<meta_t> events_nominal;
    read_event_metadata(events_nominal, nominal);
    // Variation sample.
    std::vector<meta_t> events_variation;
    read_event_metadata(events_variation, variation);

    /**
     * Create the intersection of events from the nominal and variation samples.
    */
    std::vector<meta_t> intersection;
    for(const meta_t & meta : events_nominal)
    {
        if(std::find(events_variation.begin(), events_variation.end(), meta) != events_variation.end())
            intersection.push_back(meta);
    }

    /**
     * The bootstrapping process will be significantly easier if we do some
     * work beforehand. Specifically, we will create a std::map that maps
     * event metadata (run, subrun, event) to a vector of reconstructed
     * interactions. This will allow us to quickly access the reconstructed
     * interactions for a given event.
    */
    // Nominal sample.
    std::map<index_t, std::vector<double>> reco_nominal;
    read_selected(reco_nominal, nominal);
    std::map<meta_t, std::vector<index_t>> nominal_map;
    for(const std::pair<index_t, std::vector<double>> & entry : reco_nominal)
    {
        index_t index(entry.first);
        meta_t meta(std::get<0>(index), std::get<1>(index), std::get<2>(index));
        nominal_map[meta].push_back(index);
    }   
    // Variation sample.
    std::map<index_t, std::vector<double>> reco_variation;
    read_selected(reco_variation, variation);
    std::map<meta_t, std::vector<index_t>> variation_map;
    for(const std::pair<index_t, std::vector<double>> & entry : reco_variation)
    {
        index_t index(entry.first);
        meta_t meta(std::get<0>(index), std::get<1>(index), std::get<2>(index));
        variation_map[meta].push_back(index);
    }

    /**
     * Bootstrapping step. For each bootstrapped universe, we select N events
     * from the intersection of the samples (N events total) with replacement.
     * We will use the Mersenne Twister 19937 generator to seed a uniform
     * distribution which takes on integral values in the range [0, N-1] and
     * use this to select N events. We then loop over the bootstrapped events
     * and check if they have any selected interactions in either the nominal
     * or the variation sample. For all selected interactions in bootstrapped
     * events, we fill the TH2D histograms with the reconstructed quantities
     * (X = reconstructed quantity, Y = bootstrap universe).
    */
    //std::minstd_rand gen(0); // faster, but less random than mt19937
    std::mt19937 gen(0);
    std::uniform_int_distribution<size_t> dist(0, intersection.size() - 1);
    /**
     * Create the TH2D histograms for the reconstructed quantities for both
     * the nominal and variation samples. The histograms are named following
     * the pattern <syst_name>_bootstrap_<sample_type>_<reco_var_name>.
    */
    for(const RecoVar & r : reco_vars)
    {
        std::string name = systname + "_bootstrap_nominal_" + r.name;
        weights.insert(std::make_pair(name, new TH2D(name.c_str(), name.c_str(), r.nbins, r.xmin, r.xmax, 1000, 0, 1000)));
        name = systname + "_bootstrap_variation_" + r.name;
        weights.insert(std::make_pair(name, new TH2D(name.c_str(), name.c_str(), r.nbins, r.xmin, r.xmax, 1000, 0, 1000)));
    } // End loop over the reconstructed quantities.
    
    /**
     * Begin loop over bootstrapped universes. The variables rnom and rsys will
     * store the indices of the selected interactions in the nominal and
     * variation samples, respectively.
    */
    std::vector<index_t> rnom, rsys;
    for(size_t b(0); b < 1000; ++b)
    {
        /**
         * Select N events from the intersection of the samples with
         * replacement.
        */
        for(size_t n(0); n < intersection.size(); ++n)
        {
            // Sample event.
            meta_t meta(intersection[dist(gen)]);
            
            // Clear the vectors storing selected interactions.
            rnom.clear();
            rsys.clear();

            // Check if the event has any selected interactions in the nominal sample.
            if(nominal_map.find(meta) != nominal_map.end())
                rnom = nominal_map[meta];
            
            // Check if the event has any selected interactions in the variation sample.
            if(variation_map.find(meta) != variation_map.end())
                rsys = variation_map[meta];

            /**
             * Loop over the reconstructed quantities and fill the TH2D
             * histograms with the reconstructed values for the selected
             * interactions.
            */
            for(size_t ri(0); ri < reco_vars.size(); ++ri)
            {
                RecoVar & r = reco_vars[ri];
                // Nominal sample.
                std::string name = systname + "_bootstrap_nominal_" + r.name;
                for(size_t i(0); i < rnom.size(); ++i)
                    weights[name]->Fill(reco_nominal[rnom[i]][ri], b);
                // Variation sample.
                name = systname + "_bootstrap_variation_" + r.name;
                for(size_t i(0); i < rsys.size(); ++i)
                    weights[name]->Fill(reco_variation[rsys[i]][ri], b);
            } // End loop over the reconstructed quantities.
        } // End loop over the bootstrapped events.
    } // End loop over bootstrapped universes.

    /**
     * There are several quantities that are useful to calculate for each
     * systematic/reconstructed quantity combination. These include the
     * bootstrap mean difference and the bootstrap mean ratio.
    */
    for(const RecoVar & r : reco_vars)
    {
        /**
         * Calculate the bin-to-bin difference and ratio for the nominal and
         * variation samples along with the associated central values.
        */
        std::string name = systname + "_bootstrap_nominal_" + r.name;
        TH2D * hnom = static_cast<TH2D*>(weights[name]);
        name = systname + "_bootstrap_variation_" + r.name;
        TH2D * hsys = static_cast<TH2D*>(weights[name]);
        
        // Bootstrap difference and ratio histograms.
        name = systname + "_bootstrap_diff_" + r.name;
        weights.insert(std::make_pair(name, new TH2D(name.c_str(), name.c_str(), r.nbins, r.xmin, r.xmax, 1000, 0, 1000)));
        name = systname + "_bootstrap_ratio_" + r.name;
        weights.insert(std::make_pair(name, new TH2D(name.c_str(), name.c_str(), r.nbins, r.xmin, r.xmax, 1000, 0, 1000)));

        // Central values for difference and ratio
        name = systname + "_bootstrap_diff_" + r.name + "_cv";
        weights.insert(std::make_pair(name, new TH1D(name.c_str(), name.c_str(), r.nbins, r.xmin, r.xmax)));
        name = systname + "_bootstrap_ratio_" + r.name + "_cv";
        weights.insert(std::make_pair(name, new TH1D(name.c_str(), name.c_str(), r.nbins, r.xmin, r.xmax)));

        for(size_t i(0); i < r.nbins; ++i)
        {
            for(size_t j(0); j < 1000; ++j)
            {
                double vnom = hnom->GetBinContent(i, j);
                double vsys = hsys->GetBinContent(i, j);
                double diff = vsys - vnom;
                double ratio = vsys / (vnom != 0 ? vnom : 1);
                weights[systname + "_bootstrap_diff_" + r.name]->SetBinContent(i, j, diff);
                weights[systname + "_bootstrap_ratio_" + r.name]->SetBinContent(i, j, ratio);
                weights[systname + "_bootstrap_diff_" + r.name + "_cv"]->SetBinContent(i, weights[systname + "_bootstrap_diff_" + r.name + "_cv"]->GetBinContent(i) + diff / 1000.0);
                weights[systname + "_bootstrap_ratio_" + r.name + "_cv"]->SetBinContent(i, weights[systname + "_bootstrap_ratio_" + r.name + "_cv"]->GetBinContent(i) + diff / 1000.0);
            } // End loop over the bootstrap universes.
        } // End loop over the reconstructed quantity bins.

        /**
         * Calculate the covariance matrix for the difference and ratio
         * between the nominal and variation samples for the reconstructed
         * quantity.
        */
        name = systname + "_bootstrap_diff_" + r.name;
        calc_covariance(weights, name);
        name = systname + "_bootstrap_ratio_" + r.name;
        calc_covariance(weights, name);

    } // End loop over the reconstructed quantities.
}

void cholesky_decomposition(weights_t weights, std::string systname)
{
    TH2D * cov = static_cast<TH2D*>(weights[systname + "_cov"]);
    TH1D * central = static_cast<TH1D*>(weights[systname + "_cv"]);
    std::vector<bool> nonzero(cov->GetNbinsX(), true);
    size_t nonzero_count(0);
    for(int i(1); i <= central->GetNbinsX(); ++i)
    {
        nonzero[i-1] = central->GetBinContent(i) != 0;
        ++nonzero_count;
    }
    TMatrixDSym mcov(nonzero_count);
    for(int i(1); i <= cov->GetNbinsX(); ++i)
    {
        if(nonzero[i-1])
        {
            for(int j(1); j <= cov->GetNbinsY(); ++j)
            {
                if(nonzero[j-1])
                    mcov(i-1, j-1) = cov->GetBinContent(i, j);
            }
        }
    }
    TDecompChol decomp(mcov);
    decomp.Decompose();
}

#endif