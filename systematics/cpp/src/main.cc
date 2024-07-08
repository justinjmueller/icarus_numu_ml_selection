#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TROOT.h"
#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"
#include "types.h"
#include "vars.h"
#include "utilities.h"
//#include "variation.h"
#include "reweight.h"

int main(int argc, char ** argv)
{
    /**
     * Set various input file paths. Hardcoded for now.
    */
    //std::string caf = "/pnfs/icarus/scratch/users/mueller/systematics/sample_cv.caf.root";
    std::string base_path = "/pnfs/icarus/scratch/users/mueller/mc_run2/standard_cafs/";
    std::string nominal = "/exp/icarus/app/users/mueller/sbn_ml_cafmaker/icarus_numu_ml_selection/systematics/cpp/build/output_mc_rev3.root";
    //std::string variation = "/exp/icarus/app/users/mueller/sbn_ml_cafmaker/icarus_numu_ml_selection/systematics/cpp/build/output_sigshape.root";

    /**
     * Ignore ROOT warnings (like missing libraries, which are not necessarily
     * fatal).
    */
    gErrorIgnoreLevel = kError;

    /**
     * Prepare to store the systematic weights for each reconstructed quantity.
     * The "weights_t" is typedef that maps a string (systematic name) to a
     * TH2D* object (X = reconstructed quantity, Y = systematic universe). The
     * systematic parameters are specified in vars.h.
    */
    weights_t weights;

    /**
     * Load the selected interactions from the input file. Each selected
     * interaction has a unique index_t key (reflecting event metadata) and
     * several reconstructed quantities (specified in vars.h).
    */
    std::map<index_t, std::vector<double>> reco_map;
    read_selected(reco_map, nominal);
    double POT(0);

    //calc_variation_systematics("signal_shape", nominal, variation, weights);

    /**
     * Read the list of input files from a regular text file. Each line in the
     * file should contain the path to a CAF file.
    */
    std::vector<std::string> input_files;
    std::ifstream file_list("input_files.txt");
    std::string line;
    while(std::getline(file_list, line))
        input_files.push_back(line);
    file_list.close();

    //for(int file_index(1); file_index < argc; ++file_index)
    for(size_t file_index(0); file_index < input_files.size(); ++file_index)
    {
        std::cout << "Processing file " << file_index << std::endl;
        POT += calc_reweight_systematics(base_path + input_files[file_index], reco_map, weights);
    }
    std::cout << "Total POT: " << POT << std::endl;

    /**
     * Write the TH2D histograms to the output file. Each histogram has a name
     * following the pattern <syst_name>_<reco_var_name>. The X-axis represents
     * the reconstructed quantity and the Y-axis represents the systematic
     * universe number. We do this step before closing the input file since the
     * histograms were created while reading the input file and are therefore
     * "owned" by the input file.
    */
    TFile * output = new TFile("output_1mu1p_rev2.root", "RECREATE");
    for(std::pair<std::string, TH1*> syst : weights)
        syst.second->Write();
    output->Close();

    return 0;
}
