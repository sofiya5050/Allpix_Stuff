// SPDX-FileCopyrightText: 2021-2023 CERN and the Allpix Squared authors
// SPDX-License-Identifier: MIT

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSystem.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>

#include "RtypesCore.h"

#include <stdexcept>
#include <memory>
#include <map> // To access the map object (dictionary)
#include <cmath>
#include <fstream> // To read and write to files
#include <sstream> // To get a string buffer from input stream

#include "/cvmfs/clicdp.cern.ch/software/allpix-squared/3.0.3/x86_64-centos7-gcc12-opt/include/objects/MCParticle.hpp"
#include "/cvmfs/clicdp.cern.ch/software/allpix-squared/3.0.3/x86_64-centos7-gcc12-opt/include/objects/PixelCharge.hpp"
#include "/cvmfs/clicdp.cern.ch/software/allpix-squared/3.0.3/x86_64-centos7-gcc12-opt/include/objects/PixelHit.hpp"
#include "/cvmfs/clicdp.cern.ch/software/allpix-squared/3.0.3/x86_64-centos7-gcc12-opt/include/objects/PropagatedCharge.hpp"
#include "/cvmfs/clicdp.cern.ch/software/allpix-squared/3.0.3/x86_64-centos7-gcc12-opt/include/objects/DepositedCharge.hpp"
#include "/cvmfs/clicdp.cern.ch/software/allpix-squared/3.0.3/x86_64-centos7-gcc12-opt/include/objects/MCTrack.hpp"


// These structs are used to get data from the functions below
struct chargeVectorTuple{
    std::vector<int> deposited;
    std::vector<int> collected;
};


struct langauFitReturn{
	double fitParams[4];
	double fitErrors[4];
	double chiSqr;
	int NDF;
	TF1 * fitFunc;
};

double XMAX = 70e3;


chargeVectorTuple getROOTAllpixData(TFile *file, std::string detector)
{
    // This function initializes the ROOT Trees that store the Allpix objects 
	// and gets necessary data out of them

    // Initialise reading of the PixelHit TTrees
    TTree *pixel_charge_tree = static_cast<TTree *>(file->Get("PixelCharge"));
    if (!pixel_charge_tree)
    {
		throw std::runtime_error("Could not read tree PixelCharge");
    }

    TBranch *pixel_charge_branch = pixel_charge_tree->FindBranch(detector.c_str());
    if (!pixel_charge_branch)
    {
		throw std::runtime_error( "Could not find the branch on tree PixelCharge for the corresponding detector, cannot continue");
    }

    // Bind the information to a predefined vector
    std::vector<allpix::PixelCharge *> pixel_charges;
    pixel_charge_branch->SetObject(&pixel_charges);

    // Initialise reading of the MCParticle TTrees
    TTree *mc_particle_tree = static_cast<TTree *>(file->Get("MCParticle"));
    if (!mc_particle_tree)
    {
		throw std::runtime_error("Could not read tree MCParticle");
    }
    TBranch *mc_particle_branch = mc_particle_tree->FindBranch(detector.c_str());
    if (!mc_particle_branch)
    {
		throw std::runtime_error("Could not find the branch on tree MCParticle for the corresponding detector, cannot continue");
    }
    // Bind the information to a predefined vector
    std::vector<allpix::MCParticle *> input_particles;
    mc_particle_branch->SetObject(&input_particles);

    std::vector<int> deposited_charges;
    std::vector<int> collected_charges;

    // Getting the charges by iterating over all events
    for (int i = 0; i < pixel_charge_tree->GetEntries(); ++i)
    {
        if (i % 1000 == 0)
        {
            std::cout << "Processing event " << i << std::endl;
        }

        // Access next event. Pushes information into pixel_charges, input_particles
        pixel_charge_tree->GetEntry(i);
        mc_particle_tree->GetEntry(i);

        int collected_charge = 0;
        for (auto &inp_hit : pixel_charges)
    	{
		collected_charge += inp_hit->getCharge();
	    }
        collected_charges.push_back(std::abs(collected_charge));

        int deposited_charge = 0;
        for(auto &mc_part: input_particles)
        {
            deposited_charge += mc_part->getTotalDepositedCharge();
        }
        deposited_charges.push_back(deposited_charge/2);
        
    }

    chargeVectorTuple v;
    v.collected = collected_charges;
    v.deposited = deposited_charges;

	std::cout << "\nFinished collecting data from Allpix's ROOT file." << std::endl;

    return v;
}


TH1F *fillHist(std::vector<int> collected_charges){
	// Small helper function to fill the histogram with values
	int nbins = 300, low = 0, high = XMAX;
	TH1F *hist = new TH1F("langauhist", "Collected Charge Histogram", nbins, low, high);
	
	for(size_t i = 0; i < collected_charges.size(); i++){
		hist->Fill(collected_charges[i]);
	}

	std::cout << "Filled the histogram." << std::endl;

	return hist;

}


double langaufun(double *x, double *par)
{
	// This is the actual fit function that will be used to find the MPV. It is a convolution of a Landau distribution and a Gaussian distribution.
	
	/*
	Fit parameters:
		0. Width (scale) parameter of Landau density
		1. MPV of Landau density
		2. Total area (integral -inf to inf, normalization constant)
		3. Width (sigma) of convoluted Gaussian function
	*/
	
	// Numeric constants
	double invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
	double mpshift = -0.22278298;      // Landau maximum location

	// Control constants
	double np = 100.0;
	double sc = 5.0;

	// Variables
	double xx;
	double mpc;
	double fland;
	double sum = 0.0;
	double xlow, xupp;
	double step;
	double i;

	// MP shift correction
	mpc = par[1] - mpshift * par[0];

	// Range of convolution integral
	xlow = x[0] - sc * par[3];
	xupp = x[0] + sc * par[3];

	step = (xupp - xlow) / np;

	for (i = 1.0; i <= np / 2; i++)
	{
		xx = xlow + (i - .5) * step;
		fland = TMath::Landau(xx, mpc, par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0], xx, par[3]);

		xx = xupp - (i - .5) * step;
		fland = TMath::Landau(xx, mpc, par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0], xx, par[3]);
	}

	return (par[2] * step * sum * invsq2pi / par[3]);
}


langauFitReturn langaufit(TH1F *his, double fitrange[2], double* startvalues)
{
	// does a fit on given histogram over a given range with initial guesses 

	/*
	Variables for langaufit call:
		his             histogram to fit
		fitrange[2]     lo and hi boundaries of fit range
		startvalues[4]  reasonable start values for the fit
	*/

	// Create a new fit function ROOT object
	std::string FunName = "LangauFitFunction";
	auto ffit = new TF1(FunName.c_str(), langaufun, fitrange[0], fitrange[1], 4);
	ffit->SetParNames("Width", "MP", "Area", "GSigma");
	ffit->SetParameters(startvalues);
	std::cout << "Initial parameter values set!\n" << std::endl;
	
	std::cout << "Fitting the Langau function to the histogram..." << std::endl;
	his->Fit(FunName.c_str(), "RQ SAME");   // fit within specified range, use ParLimits, do not plot

	// Retrieve best fit parameters and errors
	std::cout << "Retrieving best fit parameters..." << std::endl;
	langauFitReturn l;
	
	for(int p = 0; p < 4; p++){
		l.fitParams[p] = ffit->GetParameter(p);
		l.fitErrors[p] = ffit->GetParError(p);
	}

	double ChiSqr = ffit->GetChisquare();  // obtain chi^2
	int ndf = ffit->GetNDF();           // obtain ndf

	l.chiSqr = ChiSqr;
	l.NDF = ndf;
	l.fitFunc = ffit;

	return l;
 
}


void writeToFile(std::string WriteFile, int biasV, int VFD, double fluence, double MPV, double MPV_err, double chi2){
	// writes processed data to file

	fstream csvFile(WriteFile.c_str());

    if (csvFile.is_open())
    {
        ostringstream wholeFileContents;
        wholeFileContents << csvFile.rdbuf(); // reading data
        std::string str = wholeFileContents.str();

        std::string to_find = std::to_string(biasV) + "," + std::to_string(VFD) + "," + std::to_string(fluence);
		size_t pos = 6 ? biasV < 100: 7;
        if(str.rfind(std::to_string(biasV), pos) == std::string::npos) // Check if the combination of bias and VFD is already there in the file
        {
            csvFile << std::to_string(biasV) << ","
                    << std::to_string(VFD) << ","
                    << std::to_string(fluence) << ","
					<< std::to_string(MPV) << ","
					<< std::to_string(MPV_err) << ","
					<< std::to_string(chi2);
            csvFile << "\n";
        }
        csvFile.close();
    }

    else
    {
        cout << "Unable to open file";
    }
}


void writeChargeData(TFile *file, std::string detector, std::string WriteFile){
	// use this function if you want the raw deposited/collected charge data from each event

	fstream csvFile(WriteFile.c_str());
	chargeVectorTuple v = getROOTAllpixData(file, detector);

    if (csvFile.is_open())
    {
		csvFile << "Deposited,Collected\n";
		for(size_t i = 0; i < v.collected.size(); i++){
			csvFile << std::to_string(v.deposited[i]) << ","
			<< std::to_string(v.collected[i])
			<< "\n";
		}
    }

    else
    {
        std::cout << "Unable to open file";
    }

}


void processData(TFile *file, std::string detector, std::string CSVWriteFile, int biasV, int VFD, double fluence){
	// main function that gets called in runC_irrad.sh
	// the second argument refers to a detector name, which each detector created in Allpix has

    chargeVectorTuple v = getROOTAllpixData(file, detector);
    std::vector<int> deposited_charges = v.deposited;
    std::vector<int> collected_charges = v.collected;

	// Plot histogram of data
	TCanvas *c1 = new TCanvas("c1", "c1");
	std::string name = "bias " + std::to_string(biasV) + "V.pdf";
	TH1F *h1 = new TH1F("langauhist_image", name.c_str(), 200, 0, XMAX);
	
	for(size_t i = 0; i < collected_charges.size(); i++){
		h1->Fill(collected_charges[i]);
	}

	h1->GetXaxis()->SetTitle("Charge collected (e)");
	h1->GetYaxis()->SetTitle("Number of Events");
	h1->Draw();

	c1->SaveAs(name.c_str());
	delete c1; 
	delete h1;
	
	TCanvas *c2 = new TCanvas("c2", "c2");
	std::string name2 = "bias " + std::to_string(biasV) + "V-Fitted.pdf";

	TH1F *collected_hist = fillHist(collected_charges);
	double fitrange[2] = {0.0, XMAX};
	
	// Change this according to your requirements
	double CCE_guess;
	if(biasV < 200){CCE_guess = 10e3;}
	else if(biasV < 700){CCE_guess = 15e3;}
	else if(biasV < 1000){CCE_guess = 20e3;}

	double startvalues[4] = {1000.0, CCE_guess, 10000.0, 1000.0};
	langauFitReturn l = langaufit(collected_hist, fitrange, &startvalues[0]);

	collected_hist->GetXaxis()->SetTitle("Charge collected (e)");
	collected_hist->GetYaxis()->SetTitle("Number of Events");
	collected_hist->Draw();

	c2->SaveAs(name2.c_str());

	double chi2 = l.chiSqr/((double) l.NDF);
	
	std::cout << "Bias: " << biasV 
			<< " V, VFD: " << VFD 
			<< " V, fluence: " << fluence 
			<< ", MPV: (" <<  l.fitParams[1] << " ± " << l.fitErrors[1] << ")e, "
			<< "Reduced chi^2: " << chi2 << ".\n" << std::endl;

	writeToFile(CSVWriteFile, biasV, VFD, fluence, l.fitParams[1],  l.fitErrors[1], chi2);
}