//\brief : Simple ROOT macro to fit the EWK cross-sections
//\author: A. S. Mete <amete@cern.ch>
//\date  : 24/3/2015
//         15/7/2016 updated to include slepton-pair production
//\usage : Simply do "root -l fit_c1n2_wino.C"

#include <vector>
#include "Riostream.h"

// Simple Exponential Fit Function
Double_t expoFunc(Double_t *v, Double_t *par)
{
  return TMath::Exp(par[0]+par[1]*v[0]+par[2]*TMath::Log(v[0]));
}

// Load cross-section values as appears in the twiki
// https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections13TeVn2x1wino
// https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections13TeVn2x1hino
// https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections13TeVx1x1wino
// https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections13TeVx1x1hino
bool load_cross_sections(TString grid,
                         TString comp,
                         double* x, double* xe, 
                         double* y, double* ye)
{
  bool found=false;
  TString fileName="";
  if((grid.EqualTo("C1N2") || grid.EqualTo("C1pN2") || grid.EqualTo("C1mN2") || grid.EqualTo("C1C1")) &&
     (comp.EqualTo("wino") || comp.EqualTo("hino"))) {
    found=true;
  }
  if((grid.EqualTo("N1N2") || grid.EqualTo("CN")) && 
     comp.EqualTo("hino")) {
    found=true;
  }
  if((grid.EqualTo("SlepSlep")) &&
     (comp.EqualTo("left") || comp.EqualTo("right"))) {
    found=true;
  }

  if(found) {
    ifstream in; 
    in.open(Form("Inputs/xsec_%s_%s.txt",grid.Data(),comp.Data()));
    int index = 0;
    double a = 0., b = 0., c = 0.;
    while(true){
      in >> a >> b >> c;
      if(!in.good()) break;
      x[index] = a; y[index] = b; ye[index] = c;
      index++; 
    }
  }
  return found;
}

// Main application 
void fit_gaugino(TString grid = "C1N2", TString comp = "wino")
{
  bool saveROOTFile = true;
  bool doPrint      = true;

  double x[77]    = {0.}; // Mass
  double xe[77]   = {0.}; // Uncertainty of mass - set to zero since input
  double y[77]    = {0.}; // Cross-section
  double ye[77]   = {0.}; // Uncertainty of the cross-section
  double yUp[77]  = {0.}; // Cross-section + 1sigma
  double yeUp[77] = {0.}; // Pseudo-uncertainty on the cross-section + 1sigma
  double yDn[77]  = {0.}; // Cross-section - 1sigma
  double yeDn[77] = {0.}; // Pseudo-uncertainty on the cross-section - 1sigma

  if(!load_cross_sections(grid,comp,x,xe,y,ye)) {  // Read the tabulated cross-sections
    std::cout << "Couldn't find cross-sections for grid " << grid << " and composition " << comp << std::endl;
    std::cout << "Possible options for the grid are C1C1 and C1N2, while for composition are wino and hino..." << std::endl;
    return;
  } 

  int nPoints = 77;
  if(grid.EqualTo("SlepSlep")) nPoints =  10;
  for(unsigned int i=0; i<nPoints; ++i) { 
    if(grid.EqualTo("SlepSlep")) ye[i] *= y[i]; // SlepSlep has fractional uncertainties...

    yUp[i]  = y[i] + ye[i]; 
    yeUp[i] = yUp[i]*(ye[i]/y[i]); // Assume same fractional uncertainty as nominal to avoid discontiunities
                                   // in the borders of fit ranges
    yDn[i] = y[i] - ye[i]; 
    yeDn[i] = yDn[i]*(ye[i]/y[i]); // Same as above 
  }

  // Build the TGraphErrors that'll be used in the fit/drawing                                         
  TGraphErrors* crossSectionNom = new TGraphErrors(nPoints,x,y  ,xe,ye  ); 
  TGraphErrors* crossSectionUp  = new TGraphErrors(nPoints,x,yUp,xe,yeUp);  
  TGraphErrors* crossSectionDn  = new TGraphErrors(nPoints,x,yDn,xe,yeDn);  

  // Draw the canvas, use dummy histograms to better control the axes etc.
  TCanvas* canvas = new TCanvas("canvas","canvas",800,800);
  canvas->SetBorderSize(0);
  canvas->SetFillColor(0);
  TH1F* histoDummyTop = new TH1F("histoDummyTop","histoDummyTop",200,100,2000);
  histoDummyTop->GetYaxis()->SetTitle("#sigma [fb]");
  histoDummyTop->GetYaxis()->SetLabelSize(0.04);
  histoDummyTop->GetXaxis()->SetLabelOffset(1.5);
  histoDummyTop->GetYaxis()->SetRangeUser(1.e-3,1.e5);
  TPad*    topPad = new TPad("pTop","pTop",0,0.2,1,1);
  topPad->Draw();
  TH1F* histoDummyBot = new TH1F("histoDummyBot","histoDummyBot",200,100,2000);
  histoDummyBot->GetXaxis()->SetTitle("m_{#tilde{#chi}_{1}^{#pm},#tilde{#chi}_{2}^{0}} [GeV]");
  histoDummyBot->GetXaxis()->SetTitleSize(0.15);
  histoDummyBot->GetXaxis()->SetTitleOffset(1.);
  histoDummyBot->GetXaxis()->SetLabelSize(0.1);
  histoDummyBot->GetYaxis()->SetTitle("Actual/Fitted");
  histoDummyBot->GetYaxis()->SetTitleSize(0.13);
  histoDummyBot->GetYaxis()->SetTitleOffset(0.49);
  histoDummyBot->GetYaxis()->SetLabelSize(0.1);
  histoDummyBot->GetYaxis()->SetNdivisions(4);
  histoDummyBot->GetYaxis()->SetRangeUser(0.8,1.2);
  TPad*    botPad = new TPad("pBot","pBot",0,0.0,1,0.3);
  botPad->Draw();

  // Top pad
  topPad->cd();
  topPad->SetBottomMargin(0.15);
  histoDummyTop->Draw();
  crossSectionNom->SetMarkerSize(1.1);
  crossSectionNom->Draw("P&&same");      
  gPad->SetLogy(1);

  // Draw the legend
  TLegend* legend = new TLegend(0.4,0.8,0.9,0.9);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextSize(0.03);
  legend->AddEntry(crossSectionNom,Form("13 TeV %s %s cross-sections",grid.Data(),comp.Data()),"p");
  legend->Draw();

  // Preform the fits and draw each
  int    nFits =  10;
  TF1* funcsNom[10]; TF1* funcsUp[10]; TF1* funcsDn[10];
  int    fitColors[10] = {kBlue,kGreen,kOrange,kRed,kAzure-9,kYellow,kViolet,kTeal,kPink,kMagenta};
  int  fitBorders[11] = {100,150,200,300,400,600,800,1000,1200,1500,2000};
  if(grid.EqualTo("SlepSlep")) for(unsigned int i=0; i<10; i++) { fitBorders[i] = 50.*(i+1); }
  TString funcName  = "";
  TString funcTitle = "";

  for(unsigned int i=0; i<nFits; ++i){
    // Nominal
    funcName.Form("fit_nom_%i",i);
    funcsNom[i] = new TF1(funcName,expoFunc,fitBorders[i],fitBorders[i+1],3);
    funcsNom[i]->SetLineWidth(2);
    funcsNom[i]->SetLineColor(fitColors[i]);
    funcTitle.Form("fit_nom_%i_%i",fitBorders[i],fitBorders[i+1]);
    funcsNom[i]->SetTitle(funcTitle);
    crossSectionNom->Fit(funcName,"RQ");
    funcsNom[i]->Draw("same");
    // Up
    funcName.Form("fit_up_%i",i);
    funcsUp[i] = new TF1(funcName,expoFunc,fitBorders[i],fitBorders[i+1],3);
    funcsUp[i]->SetLineWidth(2);
    funcsUp[i]->SetLineColor(fitColors[i]);
    funcsUp[i]->SetLineStyle(7);
    funcTitle.Form("fit_up_%i_%i",fitBorders[i],fitBorders[i+1]);
    funcsUp[i]->SetTitle(funcTitle);
    crossSectionUp->Fit(funcName,"RQ");
    funcsUp[i]->Draw("same");
    // Down
    funcName.Form("fit_dn_%i",i);
    funcsDn[i] = new TF1(funcName,expoFunc,fitBorders[i],fitBorders[i+1],3);
    funcsDn[i]->SetLineWidth(2);
    funcsDn[i]->SetLineColor(fitColors[i]);
    funcsDn[i]->SetLineStyle(7);
    funcTitle.Form("fit_dn_%i_%i",fitBorders[i],fitBorders[i+1]);
    funcsDn[i]->SetTitle(funcTitle);
    crossSectionDn->Fit(funcName,"RQ");
    funcsDn[i]->Draw("same");
  }
  gPad->RedrawAxis();

  // Bottom Pad
  botPad->cd();
  botPad->SetBottomMargin(0.4);
  histoDummyBot->Draw();
  // Calculate the yellow band in the ratio plot
  double rx[190]  = {0.};
  double rxe[190] = {0.};
  double ry[190]  = {0.};
  double rye[190] = {0.};
  for(unsigned int i=0; i<190; ++i) {
    double mass = 100+i*10.;
    for(unsigned int j=0; j<nFits; ++j) {
      if(mass>=fitBorders[j]&&mass<fitBorders[j+1]) { 
        double var1 = funcsUp[j]->Eval(mass)-funcsNom[j]->Eval(mass);
        double var2 = funcsNom[j]->Eval(mass)-funcsDn[j]->Eval(mass);
        double uncertainty = var1>var2 ? var1:var2;
        rye[i] = uncertainty/funcsNom[j]->Eval(mass); 
      }
    }
    rx[i]  = mass;
    rxe[i] = 10.; // 100 to 2000, 190 steps == 10 GeV
    ry[i]  = 1.;
    //std::cout << rx[i] << " " << rxe[i] << " " << ry[i] << " " << rye[i] << std::endl;
  }
  TGraphErrors* fitBand = new TGraphErrors(190,rx,ry,rxe,rye); 
  fitBand->SetMarkerSize(0);
  fitBand->SetFillColor(kYellow);
  fitBand->Draw("same&&E2");
  gPad->SetGridy(1);
  // Calculate the markers in the ratio plot
  // only take into account the uncertainty on the actual value
  // since the uncertainty of the fit is shown explicitly
  TGraphErrors* ratio = new TGraphErrors(*crossSectionNom); 
  for(unsigned int i=0; i<ratio->GetN(); ++i) {
    double original_x = 0., original_y = 0.;
    crossSectionNom->GetPoint(i,original_x,original_y);
    double fitXsec = -1.;
    for(unsigned int j=0; j<nFits; ++j) {
      if(original_x>=fitBorders[j]&&original_x<fitBorders[j+1]) { 
        fitXsec = funcsNom[j]->Eval(original_x);
      } 
    }
    ratio->SetPoint(i,original_x,original_y/fitXsec); 
    ratio->SetPointError(i,0.,crossSectionNom->GetErrorY(i)/fitXsec); 
    //std::cout << original_x << " " << original_y << " " << fitXsec << " " << original_y/fitXsec << std::endl; 
  }
  TLine* line = new TLine(histoDummyBot->GetXaxis()->GetXmin(),1,histoDummyBot->GetXaxis()->GetXmax(),1);
  line->SetLineColor(kRed);
  line->SetLineStyle(7);
  line->SetLineWidth(2);
  line->Draw("same");
  ratio->Draw("same&&P");
  gPad->RedrawAxis();

  // Print the values if user asks for it
  if(doPrint) { 
    std::cout << "==================================================================================" << std::endl;
    std::cout << std::setw(27) << "" << grid << " " << comp << " cross-sections [fb] " << std::endl;
    std::cout << "==================================================================================" << std::endl;
    std::cout << std::setw(13) << ""
              << " ::    Actual" 
              << " -   Fitted" 
              << " - " << std::setw(8) << "" 
              << " ::   Actual" 
              << " -   Fitted" 
              << " - "  
              << std::endl;
    std::cout << "  Mass [GeV] " 
              << " ::     xsec " 
              << " -    xsec " 
              << " - Diff [%]" 
              << " ::     unc " 
              << " -     unc " 
              << " - Diff [%]" 
              << std::endl;
    std::cout << "==================================================================================" << std::endl;
    for(unsigned int i=0; i<nPoints; ++i){
      double mass=x[i];
      if(mass<1.e-3) continue;
      for(unsigned int j=0; j<nFits; ++j){
        if(mass>=fitBorders[j]&&mass<=fitBorders[j+1]) { 
          double var1 = funcsUp[j]->Eval(mass)-funcsNom[j]->Eval(mass);
          double var2 = funcsNom[j]->Eval(mass)-funcsDn[j]->Eval(mass);
          double uncertainty = var1>var2 ? var1:var2;
          std::cout << " "    << std::setw(8) << std::setprecision(5) << x[i] << std::setw(8);
          std::cout << " :: " << std::setw(9) << std::setprecision(5) << y[i];
          std::cout << " - "  << std::setw(8) << std::setprecision(5) << funcsNom[j]->Eval(mass);
          std::cout << " - "  << std::setw(8) << std::setprecision(2) << ( funcsNom[j]->Eval(mass) - y[i] )/y[i]*100.;
          std::cout << " :: " << std::setw(8) << std::setprecision(5) << ye[i];
          std::cout << " - "  << std::setw(8) << std::setprecision(5) << uncertainty;
          std::cout << " - "  << std::setw(8) << std::setprecision(2) << ( uncertainty - ye[i] )/ye[i]*100.;
          break;
        }
      }
      std::cout << std::endl;
    }
    std::cout << "==================================================================================" << std::endl;
  }

  // Save in a ROOT file if the users asks for it
  if(!saveROOTFile) return;

  std::cout << "Writing output ..." << std::endl;
  TFile* output = new TFile(Form("%s_%s_13TeV.root",grid.Data(),comp.Data()),"RECREATE");
  output->cd();
  TParameter<int> param("nFits",nFits);
  param.Write();
  for(unsigned int i=0; i<nFits; ++i) {
    funcsNom[i]->Write();
    funcsUp[i]->Write();
    funcsDn[i]->Write();
  }
  TTree* tree = new TTree("parameters","TTree with masses [GeV] and cross-sections [fb]");
  std::vector<double> xx , xxe , yy , yye ;
  tree->Branch("mass",&xx);
  tree->Branch("xsec",&yy);
  tree->Branch("xsecUnc",&yye);
  for(unsigned int i=0; i<nPoints; ++i) {
    xx.push_back(x[i]); xxe.push_back(xe[i]);
    yy.push_back(y[i]); yye.push_back(ye[i]);
  } 
  tree->Fill();
  output->Write(); 
  output->Close();
}    
