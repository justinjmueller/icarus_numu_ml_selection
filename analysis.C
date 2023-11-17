/**
 * @file analysis.C
 * @brief ROOT macro to be used with CAFAna executable to run the selection.
 * @author justin.mueller@colostate.edu
*/

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"

#include "TFile.h"
#include "TH1D.h"

#include "include/analysis.h"

using namespace ana;

void analysis()
{
    // Some useful variables for later.
    const std::string input_file = "/exp/icarus/app/users/mueller/sbn_ml_cafmaker/sbn_ml_cafmaker/build/bnb_nucosmics.flat.root";

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader loader(input_file);

    // Create the binning schemes for the Vars we wish to plot.
    const Binning bCount = Binning::Simple(1, 0, 2);
    const Binning bCountParticles = Binning::Simple(20, 0, 20);
    const Binning bIs1mu1p = Binning::Simple(2, 0, 2);

    Spectrum sCountParticles("Particle Count", bCountParticles, loader, kCountParticles, kNoSpillCut);
    Spectrum sCountPrimaries("Primary Count", bCountParticles, loader, kCountPrimaries, kNoSpillCut);

    Spectrum sCountParticlesTruth("Particle Count (Truth)", bCountParticles, loader, kCountParticlesTruth, kNoSpillCut);
    Spectrum sCountPrimariesTruth("Primary Count (Truth)", bCountParticles, loader, kCountPrimariesTruth, kNoSpillCut);

    Spectrum sCount_1mu1p_NoCut("1mu1p", bCount, loader, kCount_1mu1p_NoCut, kNoSpillCut);
    Spectrum sCount_OtherNu_NoCut("OtherNu", bCount, loader, kCount_OtherNu_NoCut, kNoSpillCut);
    Spectrum sCount_Cosmic_NoCut("Cosmic", bCount, loader, kCount_Cosmic_NoCut, kNoSpillCut);

    Spectrum sCount_1mu1p_FVCut("1mu1p", bCount, loader, kCount_1mu1p_FVCut, kNoSpillCut);
    Spectrum sCount_OtherNu_FVCut("OtherNu", bCount, loader, kCount_OtherNu_FVCut, kNoSpillCut);
    Spectrum sCount_Cosmic_FVCut("Cosmic", bCount, loader, kCount_Cosmic_FVCut, kNoSpillCut);

    Spectrum sCount_1mu1p_FVConCut("1mu1p", bCount, loader, kCount_1mu1p_FVConCut, kNoSpillCut);
    Spectrum sCount_OtherNu_FVConCut("OtherNu", bCount, loader, kCount_OtherNu_FVConCut, kNoSpillCut);
    Spectrum sCount_Cosmic_FVConCut("Cosmic", bCount, loader, kCount_Cosmic_FVConCut, kNoSpillCut);

    Spectrum sCount_1mu1p_FVConTopCut("1mu1p", bCount, loader, kCount_1mu1p_FVConTopCut, kNoSpillCut);
    Spectrum sCount_OtherNu_FVConTopCut("OtherNu", bCount, loader, kCount_OtherNu_FVConTopCut, kNoSpillCut);
    Spectrum sCount_Cosmic_FVConTopCut("Cosmic", bCount, loader, kCount_Cosmic_FVConTopCut, kNoSpillCut);

    Spectrum sCount_1mu1p_AllCut("1mu1p", bCount, loader, kCount_1mu1p_AllCut, kNoSpillCut);
    Spectrum sCount_OtherNu_AllCut("OtherNu", bCount, loader, kCount_OtherNu_AllCut, kNoSpillCut);
    Spectrum sCount_Cosmic_AllCut("Cosmic", bCount, loader, kCount_Cosmic_AllCut, kNoSpillCut);

    loader.Go();

    TFile output_file("spectra.root", "recreate");
    
    TH1D* hCountParticles;
    hCountParticles = sCountParticles.ToTH1(1);
    hCountParticles->Write("sCountParticles");
    
    TH1D* hCountPrimaries;
    hCountPrimaries = sCountPrimaries.ToTH1(1);
    hCountPrimaries->Write("sCountPrimaries");
    TH1D* hCountParticlesTruth;
    hCountParticlesTruth = sCountParticlesTruth.ToTH1(1);
    hCountParticlesTruth->Write("sCountParticlesTruth");
    TH1D* hCountPrimariesTruth;
    hCountPrimariesTruth = sCountPrimariesTruth.ToTH1(1);
    hCountPrimariesTruth->Write("sCountPrimariesTruth");

    TH1D* hCount_1mu1p_NoCut;
    hCount_1mu1p_NoCut = sCount_1mu1p_NoCut.ToTH1(sCount_1mu1p_NoCut.POT());
    hCount_1mu1p_NoCut->Write("sCount_1mu1p_NoCut");
    TH1D* hCount_OtherNu_NoCut;
    hCount_OtherNu_NoCut = sCount_OtherNu_NoCut.ToTH1(sCount_OtherNu_NoCut.POT());
    hCount_OtherNu_NoCut->Write("sCount_OtherNu_NoCut");
    TH1D* hCount_Cosmic_NoCut;
    hCount_Cosmic_NoCut = sCount_Cosmic_NoCut.ToTH1(sCount_Cosmic_NoCut.POT());
    hCount_Cosmic_NoCut->Write("sCount_Cosmic_NoCut");

    TH1D* hCount_1mu1p_FVCut;
    hCount_1mu1p_FVCut = sCount_1mu1p_FVCut.ToTH1(sCount_1mu1p_FVCut.POT());
    hCount_1mu1p_FVCut->Write("sCount_1mu1p_FVCut");
    TH1D* hCount_OtherNu_FVCut;
    hCount_OtherNu_FVCut = sCount_OtherNu_FVCut.ToTH1(sCount_OtherNu_FVCut.POT());
    hCount_OtherNu_FVCut->Write("sCount_OtherNu_FVCut");
    TH1D* hCount_Cosmic_FVCut;
    hCount_Cosmic_FVCut = sCount_Cosmic_FVCut.ToTH1(sCount_Cosmic_FVCut.POT());
    hCount_Cosmic_FVCut->Write("sCount_Cosmic_FVCut");

    TH1D* hCount_1mu1p_FVConCut;
    hCount_1mu1p_FVConCut = sCount_1mu1p_FVConCut.ToTH1(sCount_1mu1p_FVConCut.POT());
    hCount_1mu1p_FVConCut->Write("sCount_1mu1p_FVConCut");
    TH1D* hCount_OtherNu_FVConCut;
    hCount_OtherNu_FVConCut = sCount_OtherNu_FVConCut.ToTH1(sCount_OtherNu_FVConCut.POT());
    hCount_OtherNu_FVConCut->Write("sCount_OtherNu_FVConCut");
    TH1D* hCount_Cosmic_FVConCut;
    hCount_Cosmic_FVConCut = sCount_Cosmic_FVConCut.ToTH1(sCount_Cosmic_FVConCut.POT());
    hCount_Cosmic_FVConCut->Write("sCount_Cosmic_FVConCut");

    TH1D* hCount_1mu1p_FVConTopCut;
    hCount_1mu1p_FVConTopCut = sCount_1mu1p_FVConTopCut.ToTH1(sCount_1mu1p_FVConTopCut.POT());
    hCount_1mu1p_FVConTopCut->Write("sCount_1mu1p_FVConTopCut");
    TH1D* hCount_OtherNu_FVConTopCut;
    hCount_OtherNu_FVConTopCut = sCount_OtherNu_FVConTopCut.ToTH1(sCount_OtherNu_FVConTopCut.POT());
    hCount_OtherNu_FVConTopCut->Write("sCount_OtherNu_FVConTopCut");
    TH1D* hCount_Cosmic_FVConTopCut;
    hCount_Cosmic_FVConTopCut = sCount_Cosmic_FVConTopCut.ToTH1(sCount_Cosmic_FVConTopCut.POT());
    hCount_Cosmic_FVConTopCut->Write("sCount_Cosmic_FVConTopCut");

    TH1D* hCount_1mu1p_AllCut;
    hCount_1mu1p_AllCut = sCount_1mu1p_AllCut.ToTH1(sCount_1mu1p_AllCut.POT());
    hCount_1mu1p_AllCut->Write("sCount_1mu1p_AllCut");
    TH1D* hCount_OtherNu_AllCut;
    hCount_OtherNu_AllCut = sCount_OtherNu_AllCut.ToTH1(sCount_OtherNu_AllCut.POT());
    hCount_OtherNu_AllCut->Write("sCount_OtherNu_AllCut");
    TH1D* hCount_Cosmic_AllCut;
    hCount_Cosmic_AllCut = sCount_Cosmic_AllCut.ToTH1(sCount_Cosmic_AllCut.POT());
    hCount_Cosmic_AllCut->Write("sCount_Cosmic_AllCut");

    output_file.Close();
}