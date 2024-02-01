/**
 * @file data.C
 * @brief ROOT macro to be used with CAFAna executable to run the selection on
 * data events.
 * @author justin.mueller@colostate.edu
*/

#include "include/analysis.h"
#include "include/container.h"
#include "include/csv_maker.h"
#include "sbnana/CAFAna/Core/Binning.h"

using namespace ana;

/**
 * The main function of the selection. Creates a container for the CAFAna
 * Spectrum objects and populates it with a variety of variables that define
 * the selection.
 * @return none.
*/
void data()
{
    SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/run9435.flat.root", "spectra_data.root");

    spectra.add_spectrum1d("sCountParticles", Binning::Simple(20, 0, 20), kCountParticles);
    spectra.add_spectrum1d("sCountPrimaries", Binning::Simple(20, 0, 20), kCountPrimaries);

    spectra.add_spectrum1d("sDataCSV", Binning::Simple(1, 0, 2), kDataLogger);

    spectra.run();
}