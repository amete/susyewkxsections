//\brief : Simple ROOT macro to get the EWK cross-sections
//         First run the fit via "root -l fit_c1n2_wino.C"
//         if the necessary ROOT file doesn't exist
//\author: A. S. Mete <amete@cern.ch>
//\date  : 24/3/2015
//         15/7/2016 updated to include slepton-pair production

#include <stdlib.h> 

// Simple function to read the cross-section from the ROOT file
void get_gaugino( TString grid        = "C1N2", 
                  TString comp        = "wino", 
                  double mass         = 1200., 
                  bool useActual      = true, // mass of C1/N2 in GeV - use actual value if exists 
                  bool printPb        = true,
                  bool relUncertainty = true ) 
{
  std::vector<double>* dummy = new std::vector<double>();

  // pb vs fb
  double multiplier = 1.;
  if(printPb) multiplier = 1.e-3;

  // Open the input file
  TFile* input = new TFile(Form("%s_%s_13TeV.root",grid.Data(),comp.Data()),"READ");
  if(!input->IsOpen()) {
    std::cout << "Cannot find input " << grid << "_" << comp << "_13TeV.root, quitting..." << std::endl;
    return;
  }

  // See if actual points exist
  std::vector<double>* masses; std::vector<double>* xsecs; std::vector<double>* xsecUncs;
  TTree* tree = (TTree*) input->Get("parameters");
  if(!tree) {
    std::cout << "Cannot find input tree, quitting..." << std::endl;
    return;
  }

  tree->SetBranchAddress("mass",&masses);
  tree->SetBranchAddress("xsec",&xsecs);
  tree->SetBranchAddress("xsecUnc",&xsecUncs);
  tree->GetEntry(0);
 
  for(unsigned int i=0; i<masses->size() && useActual; ++i) {
    if(TMath::Abs(masses->at(i)-mass)<1.e-3) {
      //std::cout << "\nCross-section for mass " << mass << " [GeV] is " << xsecs->at(i) << " +/- " << xsecUncs->at(i) << " [fb]\n" << std::endl;
      std::cout << "\nCross-section for mass " << mass << " [GeV] is " << xsecs->at(i)*multiplier;      
      if(printPb) { std::cout << " [pb]"; } else { std::cout << " [fb]"; }
      if(relUncertainty) {
          std::cout << " +/- " << xsecUncs->at(i)/xsecs->at(i) << " [rel.]\n" << std::endl;
      } else {
        std::cout << " +/- " << xsecUncs->at(i)*multiplier;
        if(printPb) { std::cout << " [pb]\n" << std::endl; } else { std::cout << " [fb]\n" << std::endl; }
      }
      return;
    }
  }

  // Fit cannot find, move on to fits
  // Read the parameter
  TParameter<int>* param = (TParameter<int>*) input->Get("nFits");
  if(!param) {
    std::cout << "Cannot find parameter nFits, quitting..." << std::endl;
    return;
  }
  int nFits = param->GetVal();

  // Find the cross-section
  TString fitName="";
  TString fitTitle="";
  TObjArray *tokens;
  double crossSectionNom = 0., crossSectionUp = 0., crossSectionDn = 0., crossSectionUnc = 0.;
  for(unsigned int i=0; i<nFits; ++i) {
    // Read the functions
    fitName.Form("fit_nom_%i",i);
    TF1* funcNom = (TF1*) input->Get(fitName);
    fitName.Form("fit_up_%i",i);
    TF1* funcUp = (TF1*) input->Get(fitName);
    fitName.Form("fit_dn_%i",i);
    TF1* funcDn = (TF1*) input->Get(fitName);
    // Find the title that holds the min_max values
    fitTitle =  funcNom->GetTitle();
    tokens   = fitTitle.Tokenize("_");
    double min = atoi(tokens->At(2)->GetName()); 
    double max = atoi(tokens->At(3)->GetName()); 
    // Check to see if the mass matches
    if(mass >= min && mass < max) {
      crossSectionNom = funcNom->Eval(mass);
      crossSectionUp  = funcUp->Eval(mass);
      crossSectionDn  = funcDn->Eval(mass);
      crossSectionUnc = crossSectionUp-crossSectionNom > crossSectionNom-crossSectionDn ? crossSectionUp-crossSectionNom : crossSectionNom-crossSectionDn;
      break;
    } 
  } 

  // Print
  if(crossSectionNom>0) {
    //std::cout << "\nCross-section for mass " << mass << " [GeV] is " << crossSectionNom << " +/- " << crossSectionUnc << " [fb]\n" << std::endl;
    std::cout << "\nCross-section for mass " << mass << " [GeV] is " << crossSectionNom*multiplier;      
    if(printPb) { std::cout << " [pb]"; } else { std::cout << " [fb]"; }
    if(relUncertainty) {
        std::cout << " +/- " << crossSectionUnc/crossSectionNom << " [rel.]\n" << std::endl;
    } else {
      std::cout << " +/- " << crossSectionUnc*multiplier;
      if(printPb) { std::cout << " [pb]\n" << std::endl; } else { std::cout << " [fb]\n" << std::endl; }
    }
  }
  else
    std::cout << "\nCouldn't find cross-section for mass " << mass << " [GeV] \n" << std::endl;

  // Close and delete
  input->Close(); delete input;
}
