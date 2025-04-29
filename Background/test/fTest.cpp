#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <sstream>

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

#include "RooChi2Var.h"
#include "RooNLLVar.h"
#include "RooMinimizer.h"

#include "TFile.h"
#include "TMath.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooHist.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "RooMsgService.h"
#include "RooDataHist.h"
#include "RooExtendPdf.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TMacro.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TArrow.h"
#include "TKey.h"
#include "TGaxis.h"
#include "RooAddPdf.h"

#include "RooCategory.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"

#include "../interface/PdfModelBuilder.h"
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFunc.h>
#include <iomanip>
#include <limits> // Added for std::numeric_limits
#include "boost/program_options.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"

#include "../../tdrStyle/tdrstyle.C"
#include "../../tdrStyle/CMS_lumi.C"

using namespace std;
using namespace RooFit;
using namespace boost;

namespace po = program_options;

bool BLIND = true;
bool BLIND_FIT = false; // Flag to decide whether to blind the fit process
bool runFtestCheckWithToys=false;
int mgg_low = 95;
int mgg_high = 170;
int nBinsForMass = 4*(mgg_high-mgg_low);
double blind_low = 120; // Lower boundary of the blind region
double blind_high = 130; // Upper boundary of the blind region

RooRealVar *intLumi_ = new RooRealVar("IntLumi","hacked int lumi", 1000.);

TRandom3 *RandomGen = new TRandom3();

// Check if JSON parameter file exists
bool jsonParamsFileExists(const string& cat, const string& type, int order) {
  std::string fileName = Form("params_%s_%s%d.json", cat.c_str(), type.c_str(), order);
  std::string filePath = "params/" + fileName;
  
  // Check if the file exists
  std::ifstream f(filePath);
  return f.good();
}

// Load parameters from JSON file
bool loadParamsFromJSON(const string& cat, const string& type, int order, map<string, double>& params) {
  std::string fileName = Form("params_%s_%s%d.json", cat.c_str(), type.c_str(), order);
  std::string filePath = "params/" + fileName;
  
  try {
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(filePath, pt);
    
    // Iterate over all parameters in JSON
    for (const auto& item : pt) {
      params[item.first] = pt.get<double>(item.first);
    }
    std::cout << "[INFO] Loaded parameters from " << filePath << std::endl;
    return true;
  }
  catch (const std::exception& e) {
    std::cerr << "[WARNING] Failed to load parameters from " << filePath << ": " << e.what() << std::endl;
    return false;
  }
}

// Save parameters to JSON file
bool saveParamsToJSON(const string& cat, const string& type, int order, const RooAbsPdf* pdf) {
  // Create params directory (if it doesn't exist)
  system("mkdir -p params");
  
  std::string fileName = Form("params_%s_%s%d.json", cat.c_str(), type.c_str(), order);
  std::string filePath = "params/" + fileName;
  
  try {
    // Get PDF parameters
    RooArgSet* params = pdf->getParameters(RooArgSet());
    boost::property_tree::ptree pt;
    
    // Iterate through parameters and add to JSON tree
    TIterator* iter = params->createIterator();
    RooRealVar* param;
    while ((param = dynamic_cast<RooRealVar*>(iter->Next()))) {
      if (!param->isConstant()) {  // Only save floating parameters
        pt.put(param->GetName(), param->getVal());
      }
    }
    delete iter;
    delete params;
    
    // Write to JSON file
    boost::property_tree::write_json(filePath, pt);
    std::cout << "[INFO] Saved parameters to " << filePath << std::endl;
    return true;
  }
  catch (const std::exception& e) {
    std::cerr << "[ERROR] Failed to save parameters to " << filePath << ": " << e.what() << std::endl;
    return false;
  }
}

// Extract the suffix part of a parameter name (e.g. "cp0" from "env_pdf_0_13TeV_lau2_cp0")
std::string getParamSuffix(const std::string& paramName) {
  size_t lastUnderscorePos = paramName.rfind("_");
  if (lastUnderscorePos != std::string::npos) {
    return paramName.substr(lastUnderscorePos + 1);
  }
  return paramName; // Return the full name if no underscore found
}

// Set PDF parameter values
void setPdfParams(RooAbsPdf* pdf, const map<string, double>& params) {
  RooArgSet* pdfParams = pdf->getParameters(RooArgSet());
  TIterator* iter = pdfParams->createIterator();
  RooRealVar* param;
  
  // First pass: try exact parameter name match
  while ((param = dynamic_cast<RooRealVar*>(iter->Next()))) {
    if (!param->isConstant()) {  // Only set floating parameters
      auto it = params.find(param->GetName());
      if (it != params.end()) {
        std::cout << "[INFO] Setting parameter " << param->GetName() << " to " << it->second << std::endl;
        param->setVal(it->second);
      }
    }
  }
  
  // Second pass: try matching the suffix part (e.g., "cp0", "cp1", etc.)
  // For parameters like "env_pdf_0_13TeV_lau2_cp0", the order number (2) may vary
  // but parameters with same suffix (cp0) should use the same value
  delete iter;
  iter = pdfParams->createIterator();
  while ((param = dynamic_cast<RooRealVar*>(iter->Next()))) {
    if (!param->isConstant()) {  // Only set floating parameters
      auto it = params.find(param->GetName());
      if (it == params.end()) {  // If not already set in first pass
        std::string suffix = getParamSuffix(param->GetName());
        
        // Try to find a parameter with matching suffix
        for (const auto& p : params) {
          std::string paramSuffix = getParamSuffix(p.first);
          if (paramSuffix == suffix) {
            std::cout << "[INFO] Setting parameter " << param->GetName() 
                      << " to " << p.second << " (matched by suffix '" << suffix << "')" << std::endl;
            param->setVal(p.second);
            break;
          }
        }
      }
    }
  }
  
  delete iter;
  delete pdfParams;
}

RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, string type, int order, string *typePrefix, const char* ext="", const string& cat=""){
  RooAbsPdf* pdf = nullptr;
  
  // First create the PDF instance
  if (type=="Bernstein") {
    pdf = pdfsModel.getBernstein(Form("%s_bern%d",ext,order),order);
    *typePrefix = "bern"; 
  }
  else if (type=="Chebychev") {
    pdf = pdfsModel.getChebychev(Form("%s_cheb%d",ext,order),order);
    *typePrefix = "cheb";
  }
  else if (type=="Exponential") {
    pdf = pdfsModel.getExponentialSingle(Form("%s_exp%d",ext,order),order);
    *typePrefix = "exp";
  }
  else if (type=="PowerLaw") {
    pdf = pdfsModel.getPowerLawSingle(Form("%s_pow%d",ext,order),order);
    *typePrefix = "pow";
  }
  else if (type=="Laurent") {
    pdf = pdfsModel.getLaurentSeries(Form("eriesd",ext,order),order);
    *typePrefix = "lau";
  }
  else if (type=="BernsteinStepxGau") {
    pdf = pdfsModel.getBernsteinStepxGau(Form("%s_bern%d",ext,order),order);
    *typePrefix = "bern";
  }
  else if (type=="ExponentialStepxGau") {
    pdf = pdfsModel.getExponentialStepxGau(Form("%s_exp%d",ext,order),order);
    *typePrefix = "exp";
  }
  else if (type=="PowerLawStepxGau") {
    pdf = pdfsModel.getPowerLawStepxGau(Form("%s_pow%d",ext,order),order);
    *typePrefix = "pow";
  }
  else if (type=="LaurentStepxGau") {
    pdf = pdfsModel.getLaurentStepxGau(Form("%s_lau%d",ext,order),order);
    *typePrefix = "lau";
  }
  else if (type=="ExpModGauss"){
    if (order == 1) {
      pdf = pdfsModel.getExpModGaussian(Form("%s_expmodgauss",ext));
      *typePrefix = "expmodgauss";
    }
    else return NULL;
  }
  else if (type=="AsymGenGauss"){
    if (order == 1) {
      pdf = pdfsModel.getAsymGenGaussian(Form("%s_asymgauss",ext));
      *typePrefix = "asymgauss";
    }
    else return NULL;
  }
  else {
    cerr << "[ERROR] -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
  
  // If the PDF is successfully created and a category name is provided, try to load existing parameters
  if (pdf && cat != "") {
    map<string, double> params;
    bool paramsLoaded = false;
    
    // Try to load parameters for the given order
    if (jsonParamsFileExists(cat, *typePrefix, order)) {
      paramsLoaded = loadParamsFromJSON(cat, *typePrefix, order, params);
    } 
    // If not found, try to load parameters for order-1
    else if (order > 1 && jsonParamsFileExists(cat, *typePrefix, order-1)) {
      std::cout << "[INFO] Parameters for order " << order << " not found, trying order " << order-1 << std::endl;
      paramsLoaded = loadParamsFromJSON(cat, *typePrefix, order-1, params);
    }
    // Finally, try to load parameters for order-2
    else if (order > 2 && jsonParamsFileExists(cat, *typePrefix, order-2)) {
      std::cout << "[INFO] Parameters for orders " << order << " and " << order-1 
                << " not found, trying order " << order-2 << std::endl;
      paramsLoaded = loadParamsFromJSON(cat, *typePrefix, order-2, params);
    }
    
    // If parameters are successfully loaded, set them to the PDF
    if (paramsLoaded) {
      setPdfParams(pdf, params);
    }
  }
  
  return pdf;
}

RooFitResult* Chi2Fit(RooChi2Var chi2_){
  // Your minimization settings (assuming these variables are defined)
  int printLevel_ = -1; // Example value
  double eps_ = 100;    // Example value
  bool offSet_ = false;  // Example value
  int strategy_ = 0;   // Example value
  // Create the minimizer instance using the chi2 variable calculated in sidebands
  RooMinimizer mini(chi2_);
  mini.setPrintLevel(printLevel_);
  mini.setEps(eps_);
  mini.setOffsetting(offSet_);
  mini.setStrategy(strategy_);
  // Perform the minimization
  mini.minimize("Minuit2", "migrad"); // Or "Minuit" if preferred/needed
  mini.hesse(); // Calculate uncertainties using HESSE
  // Optional: Get the fit result object
  RooFitResult *fitTest = mini.save("fitResult_Sidebands", "Fit Result from Sideband Chi2");
  if (fitTest) {
      fitTest->Print("v"); // Print detailed fit result
  }
  return fitTest;
}

/**
 * Perform fit and calculate chi2, save parameters if necessary
 * Supports both RooDataSet and RooDataHist types of data
 */
void runFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries, bool doBlind=false, const std::string& cat = "", const std::string& typePrefix = "", int order = -1, int nBins = -1, const RooAbsPdf* pdfToSave = nullptr) {
  // Initialize variables
  int ntries = 0;
  int stat = 1;
  double minnll = 10e8;
  RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
  
  // chi2 calculation variables
  double chi2_before = std::numeric_limits<double>::max();
  RooRealVar* mass = nullptr;
  std::string data_hist_name;
  
  // Calculate chi2 before fit
  if (nBins > 0 && !cat.empty()) {
    mass = (RooRealVar*)data->get()->find("CMS_hgg_mass");
    if (mass) {
      RooPlot* frame_before = mass->frame();
      data->plotOn(frame_before, Binning(nBins));
      data_hist_name = frame_before->getHist()->GetName();
      pdf->plotOn(frame_before, Name("pdf_before"));
      
      int npars_before = pdf->getParameters(*data)->getSize();
      chi2_before = frame_before->chiSquare("pdf_before", data_hist_name.c_str(), npars_before);
      std::cout << "[INFO] Pre-fit chi2 = " << chi2_before << " for " << pdf->GetName() << std::endl;
      
      delete frame_before;
    } else {
      nBins = -1; // Disable chi2 calculation if mass variable is not found
    }
  }

  // Perform fit
  while (stat != 0 && ntries < MaxTries) {
    RooFitResult *fitTest = nullptr;
    
    if (doBlind && mass) {
      // Save original range
      double origMin = mass->getMin();
      double origMax = mass->getMax();
      
      // Set blind region range
      mass->setRange("lowSideband", origMin, blind_low);
      mass->setRange("highSideband", blind_high, origMax);
      
      // Fit only to sideband regions
      fitTest = pdf->fitTo(*data,
                          RooFit::Save(1),
                          RooFit::Range("lowSideband,highSideband"),
                          RooFit::Minimizer("Minuit2","minimize"),
                          RooFit::SumW2Error(kTRUE),
                          RooFit::PrintLevel(-1));

      // RooDataHist* data_hist = data->binnedClone();
      // RooChi2Var chi2_("chi2Var", "chi2Var", *pdf, *data_hist, RooFit::Range("lowSideband,highSideband"), RooFit::Save(true), RooFit::SumW2Error(kTRUE), RooFit::Extended(true));
      // fitTest = Chi2Fit(chi2_);
      
      // Restore original range
      mass->setRange(origMin, origMax);

    } else {
      // Regular fit
      fitTest = pdf->fitTo(*data,
                          RooFit::Save(1),
                          RooFit::Minimizer("Minuit2","minimize"),
                          RooFit::SumW2Error(kTRUE),
                          RooFit::PrintLevel(-1));
      
      // RooDataHist* data_hist = data->binnedClone();
      // RooChi2Var chi2_("chi2Var", "chi2Var", *pdf, *data_hist, RooFit::Save(true), RooFit::SumW2Error(kTRUE), RooFit::Extended(true));
      // fitTest = Chi2Fit(chi2_);

    }
    
    if (fitTest) {
      stat = fitTest->status();
      minnll = fitTest->minNll();

      // Calculate chi2 after fit and save parameters if improved
      if (stat == 0 && nBins > 0 && !cat.empty() && mass) {
        RooPlot* frame_after = mass->frame();
        data->plotOn(frame_after, Binning(nBins));
        pdf->plotOn(frame_after, Name("pdf_after"));
        
        int npars_after = pdf->getParameters(*data)->getSize();
        double chi2_after = frame_after->chiSquare("pdf_after", data_hist_name.c_str(), npars_after);
        
        // Save parameters if fit improved chi2
        if (chi2_after < chi2_before) {
          std::cout << "[INFO] Fit improved chi2: " << chi2_before << " -> " << chi2_after 
                    << " for " << typePrefix << order << std::endl;
          const RooAbsPdf* targetPdf = pdfToSave ? pdfToSave : pdf;
          saveParamsToJSON(cat, typePrefix, order, targetPdf);
        }
        
        delete frame_after;
      }

      // If fit failed, randomize parameters and try again
      if (stat != 0) {
        params_test->assignValueOnly(fitTest->randomizePars());
      }
      
      delete fitTest;
    } else {
      stat = 1; // Mark as failed
    }
    
    ntries++;
  }
  
  // Set output values
  *stat_t = stat;
  *NLL = minnll;
}

void runFit(RooAbsPdf *pdf, RooDataHist *data, double *NLL, int *stat_t, int MaxTries, bool doBlind=false, const std::string& cat = "", const std::string& typePrefix = "", int order = -1, int nBins = -1, const RooAbsPdf* pdfToSave = nullptr) {
  // Initialize variables
  int ntries = 0;
  int stat = 1;
  double minnll = 10e8;
  RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
  
  // chi2 calculation variables
  double chi2_before = std::numeric_limits<double>::max();
  RooRealVar* mass = nullptr;
  std::string data_hist_name;
  
  // Calculate chi2 before fit
  if (nBins > 0 && !cat.empty()) {
  mass = (RooRealVar*)data->get()->find("CMS_hgg_mass");
  if (mass) {
    RooPlot* frame_before = mass->frame();
    data->plotOn(frame_before, Binning(nBins));
    data_hist_name = frame_before->getHist()->GetName();
    pdf->plotOn(frame_before, Name("pdf_before"));
    
    int npars_before = pdf->getParameters(*data)->getSize();
    chi2_before = frame_before->chiSquare("pdf_before", data_hist_name.c_str(), npars_before);
    std::cout << "[INFO] Pre-fit chi2 = " << chi2_before << " for " << pdf->GetName() << std::endl;
    
    delete frame_before;
  } else {
    nBins = -1; // Disable chi2 calculation if mass variable is not found
  }
  }

  // Perform fit
  while (stat != 0 && ntries < MaxTries) {
  RooFitResult *fitTest = nullptr;
  
  if (doBlind && mass) {
    // Save original range
    double origMin = mass->getMin();
    double origMax = mass->getMax();
    
    // Set blind region range
    mass->setRange("lowSideband", origMin, blind_low);
    mass->setRange("highSideband", blind_high, origMax);
    
    // For RooDataHist, use Chi2 method to fit sideband regions
    RooChi2Var chi2_("chi2Var", "chi2Var", *pdf, *data, 
            RooFit::Range("lowSideband,highSideband"), 
            RooFit::Save(true), RooFit::SumW2Error(kTRUE));
    fitTest = Chi2Fit(chi2_);
    
    // Restore original range
    mass->setRange(origMin, origMax);
  } else {
    // Regular fit - use Chi2 method for RooDataHist
    RooChi2Var chi2_("chi2Var", "chi2Var", *pdf, *data, 
            RooFit::Save(true), RooFit::SumW2Error(kTRUE));
    fitTest = Chi2Fit(chi2_);
  }
  
  if (fitTest) {
    stat = fitTest->status();
    minnll = fitTest->minNll();

    // Calculate chi2 after fit and save parameters if improved
    if (stat == 0 && nBins > 0 && !cat.empty() && mass) {
    RooPlot* frame_after = mass->frame();
    data->plotOn(frame_after, Binning(nBins));
    pdf->plotOn(frame_after, Name("pdf_after"));
    
    int npars_after = pdf->getParameters(*data)->getSize();
    double chi2_after = frame_after->chiSquare("pdf_after", data_hist_name.c_str(), npars_after);
    
    // Save parameters if fit improved chi2
    if (chi2_after < chi2_before) {
      std::cout << "[INFO] Fit improved chi2: " << chi2_before << " -> " << chi2_after 
          << " for " << typePrefix << order << std::endl;
      saveParamsToJSON(cat, typePrefix, order, pdf);
    }
    
    delete frame_after;
    }

    // If fit failed, randomize parameters and try again
    if (stat != 0) {
    params_test->assignValueOnly(fitTest->randomizePars());
    }
    
    delete fitTest;
  } else {
    stat = 1; // Mark as failed
  }
  
  ntries++;
  }
  
  // Set output values
  *stat_t = stat;
  *NLL = minnll;
}

double getProbabilityFtest(double chi2, int ndof,RooAbsPdf *pdfNull, RooAbsPdf *pdfTest, RooRealVar *mass, RooDataSet *data, std::string name){
  mgg_low = mass->getMin();
  mgg_high = mass->getMax();
  nBinsForMass = 4*(mgg_high-mgg_low);
 
  double prob_asym = TMath::Prob(chi2,ndof);
  if (!runFtestCheckWithToys) return prob_asym;

  int ndata = data->sumEntries();
  
  // fit the pdfs to the data and keep this fit Result (for randomizing)
  RooFitResult *fitNullData = pdfNull->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
				,RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE),RooFit::PrintLevel(-1)); //FIXME
  RooFitResult *fitTestData = pdfTest->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
				,RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE),RooFit::PrintLevel(-1)); //FIXME

  // Ok we want to check the distribution in toys then 
  // Step 1, cache the parameters of each pdf so as not to upset anything 
  RooArgSet *params_null = pdfNull->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_null;
  params_null->snapshot(preParams_null);
  RooArgSet *params_test = pdfTest->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_test;
  params_test->snapshot(preParams_test);
 
  int ntoys =500;
  TCanvas *can = new TCanvas();
  can->SetLogy();
  TH1F toyhist(Form("toys_fTest_%s.pdf",pdfNull->GetName()),";Chi2;",60,-2,10);
  TH1I toyhistStatN(Form("Status_%s.pdf",pdfNull->GetName()),";FitStatus;",8,-4,4);
  TH1I toyhistStatT(Form("Status_%s.pdf",pdfTest->GetName()),";FitStatus;",8,-4,4);

  TGraph *gChi2 = new TGraph();
  gChi2->SetLineColor(kGreen+2);
  double w = toyhist.GetBinWidth(1);

  int ipoint=0;

  for (int b=0;b<toyhist.GetNbinsX();b++){
		double x = toyhist.GetBinCenter(b+1);
		if (x>0){
		  gChi2->SetPoint(ipoint,x,(ROOT::Math::chisquared_pdf(x,ndof)));
		  ipoint++;
		}
  }
  int npass =0; int nsuccesst =0;
  mass->setBins(nBinsForMass);
  for (int itoy = 0 ; itoy < ntoys ; itoy++){
    params_null->assignValueOnly(preParams_null);
    params_test->assignValueOnly(preParams_test);
    RooDataHist *binnedtoy = pdfNull->generateBinned(RooArgSet(*mass),ndata,0,1);

		int stat_n=1;
        int stat_t=1;
		int ntries = 0;
		double nllNull,nllTest;
		// Iterate on the fit 
		int MaxTries = 3;
		while (stat_n!=0){
		  if (ntries>=MaxTries) break;
		  RooFitResult *fitNull = pdfNull->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1),RooFit::SumW2Error(kTRUE) //FIXME
				,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1));
				//,RooFit::Optimize(0));

		  nllNull = fitNull->minNll();
      stat_n = fitNull->status();
		  if (stat_n!=0) params_null->assignValueOnly(fitNullData->randomizePars());
		  ntries++; 
		}
		
		ntries = 0;
		while (stat_t!=0){
		  if (ntries>=MaxTries) break;
		  RooFitResult *fitTest = pdfTest->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1),RooFit::SumW2Error(kTRUE) //FIXME
				,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1));
		  nllTest = fitTest->minNll();
      stat_t = fitTest->status();
		  if (stat_t!=0) params_test->assignValueOnly(fitTestData->randomizePars()); 
		  ntries++; 
		}
       
		toyhistStatN.Fill(stat_n);
		toyhistStatT.Fill(stat_t);

    if (stat_t !=0 || stat_n !=0) continue;
		nsuccesst++;
		double chi2_t = 2*(nllNull-nllTest);
		if (chi2_t >= chi2) npass++;
    toyhist.Fill(chi2_t);
  }

  double prob=0;
  if (nsuccesst!=0)  prob = (double)npass / nsuccesst;
  toyhist.Scale(1./(w*toyhist.Integral()));
  toyhist.Draw();
  TArrow lData(chi2,toyhist.GetMaximum(),chi2,0);
  lData.SetLineWidth(2);
  lData.Draw();
  gChi2->Draw("L");
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->DrawLatex(0.1,0.91,Form("Prob (asymptotic) = %.4f (%.4f)",prob,prob_asym));
  can->SaveAs(name.c_str());

  TCanvas *stas =new TCanvas();
  toyhistStatN.SetLineColor(2);
  toyhistStatT.SetLineColor(1); 
  TLegend *leg = new TLegend(0.2,0.6,0.4,0.87); leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(&toyhistStatN,"Null Hyp","L");
  leg->AddEntry(&toyhistStatT,"Test Hyp","L");
  toyhistStatN.Draw();
  toyhistStatT.Draw("same");
  leg->Draw();
  stas->SaveAs(Form("%s_fitstatus.pdf",name.c_str()));
  //reassign params
  params_null->assignValueOnly(preParams_null);
  params_test->assignValueOnly(preParams_test);

  delete can; delete stas;
  delete gChi2;
  delete leg;
  delete lat;

  // Still return the asymptotic prob (usually its close to the toys one)
  return prob_asym;

}

double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, std::string name){
  mgg_low = mass->getMin();
  mgg_high = mass->getMax();
  nBinsForMass = 4*(mgg_high-mgg_low);

  double prob;
  int ntoys = 100;
  // Routine to calculate the goodness of fit. 
  name+="_gofTest.pdf";
  RooRealVar norm("norm","norm",data->sumEntries(),0,10E6);
  //norm.removeRange();

  RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*mpdf,norm);

  // get The Chi2 value from the data
  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,Binning(nBinsForMass),Name("data"));

  pdf->plotOn(plot_chi2,Name("pdf"));
  int np = pdf->getParameters(*data)->getSize();

  double chi2 = plot_chi2->chiSquare("pdf","data",np);
  std::cout << "[INFO] Calculating GOF for pdf " << pdf->GetName() << ", using " <<np << " fitted parameters" <<std::endl;

  // The first thing is to check if the number of entries in any bin is < 5 
  // if so, we don't rely on asymptotic approximations
 
  if ((double)data->sumEntries()/nBinsForMass < 0 ){ // 5

    std::cout << "[INFO] Running toys for GOF test " << std::endl;
    // store pre-fit params 
    RooArgSet *params = pdf->getParameters(*data);
    RooArgSet preParams;
    params->snapshot(preParams);
    int ndata = data->sumEntries();
 
    int npass =0;
    std::vector<double> toy_chi2;
    for (int itoy = 0 ; itoy < ntoys ; itoy++){
    //  std::cout << "[INFO] " <<Form("\t.. %.1f %% complete\r",100*float(itoy)/ntoys) << std::flush;
      params->assignValueOnly(preParams);
      int nToyEvents = RandomGen->Poisson(ndata);
      RooDataHist *binnedtoy = pdf->generateBinned(RooArgSet(*mass),nToyEvents,0,1);
      pdf->fitTo(*binnedtoy,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Strategy(0),RooFit::SumW2Error(kTRUE)); //FIXME

      RooPlot *plot_t = mass->frame();
      binnedtoy->plotOn(plot_t);
      pdf->plotOn(plot_t);//,RooFit::NormRange("fitdata_1,fitdata_2"));

      double chi2_t = plot_t->chiSquare(np);
      if( chi2_t>=chi2) npass++;
      toy_chi2.push_back(chi2_t*(nBinsForMass-np));
      delete plot_t;
    }
    std::cout << "[INFO] complete" << std::endl;
    prob = (double)npass / ntoys;

    TCanvas *can = new TCanvas();
    double medianChi2 = toy_chi2[(int)(((float)ntoys)/2)];
    double rms = TMath::Sqrt(medianChi2);

    TH1F toyhist(Form("gofTest_%s.pdf",pdf->GetName()),";Chi2;",50,medianChi2-5*rms,medianChi2+5*rms);
    for (std::vector<double>::iterator itx = toy_chi2.begin();itx!=toy_chi2.end();itx++){
      toyhist.Fill((*itx));
    }
    toyhist.Draw();

    TArrow lData(chi2*(nBinsForMass-np),toyhist.GetMaximum(),chi2*(nBinsForMass-np),0);
    lData.SetLineWidth(2);
    lData.Draw();
    can->SaveAs(name.c_str());

    // back to best fit 	
    params->assignValueOnly(preParams);
  } else {
    prob = TMath::Prob(chi2*(nBinsForMass-np),nBinsForMass-np);
  }
  std::cout << "[INFO] Chi2 in Observed =  " << chi2*(nBinsForMass-np) << std::endl;
  std::cout << "[INFO] p-value  =  " << prob << std::endl;
  delete pdf;
  return prob;

}

void plot(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, string name,vector<string> flashggCats_, int status, double *prob){
  std::cout << "[INFO] Plotting mass range: " << mass->getMin() << " - " << mass->getMax() << std::endl;
  // Chi2 taken from full range fit
  RooPlot *plot_chi2 = mass->frame();
  mgg_low = mass->getMin();
  mgg_high = mass->getMax();
  nBinsForMass = 4*(mgg_high-mgg_low);
  data->plotOn(plot_chi2,Binning(nBinsForMass));
  pdf->plotOn(plot_chi2);

  int np = pdf->getParameters(*data)->getSize()+1; //Because this pdf has no extend
  double chi2 = plot_chi2->chiSquare(np);
 
  *prob = getGoodnessOfFit(mass,pdf,data,name);
  RooPlot *plot = mass->frame();
  mass->setRange("unblindReg_1",mgg_low,120);
  mass->setRange("unblindReg_2",130,mgg_high);
  if (BLIND) {
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(mgg_high-mgg_low),Invisible());
  }
  else data->plotOn(plot,Binning(nBinsForMass));

 // data->plotOn(plot,Binning(mgg_high-mgg_low));
  TCanvas *canv = new TCanvas();
  pdf->plotOn(plot);//,RooFit::NormRange("fitdata_1,fitdata_2"));
  pdf->paramOn(plot, RooFit::Layout(0.34, 0.96, 0.89), RooFit::Format("NEA", AutoPrecision(1)));
  plot->getAttText()->SetTextSize(0.025);
  if (BLIND) plot->SetMinimum(0.0001);
  plot->SetTitle("");
  plot->Draw();
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->DrawLatex(0.1,0.92,Form("#chi^{2} = %.3f, Prob = %.2f, Fit Status = %d ",chi2*(nBinsForMass-np),*prob,status));
  canv->SaveAs(name.c_str());
 	
	//plot_chi2->Draw();
  //canv->SaveAs((name+"debug").c_str());

  delete canv;
  delete lat;
}

void plot(RooRealVar *mass, RooMultiPdf *pdfs, RooCategory *catIndex, RooDataSet *data, string name, vector<string> flashggCats_, int cat, int bestFitPdf=-1){
  
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
  TLegend *leg = new TLegend(0.5,0.55,0.92,0.88);
  leg->SetFillColor(0);
  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1",mgg_low,120);
  mass->setRange("unblindReg_2",130,mgg_high);
  if (BLIND) {
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(mgg_high-mgg_low),Invisible());
  }
  else data->plotOn(plot,Binning(4*(mgg_high-mgg_low))); 
  TCanvas *canv = new TCanvas();
  ///start extra bit for ratio plot///
  RooHist *plotdata = (RooHist*)plot->getObject(plot->numItems()-1);
  bool doRatioPlot_=1;
  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.35);
  pad1->SetBottomMargin(0.18);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.25);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  // enf extra bit for ratio plot///

  int currentIndex = catIndex->getIndex();
  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,Form("Data - %s",flashggCats_[cat].c_str()),"LEP");
  int style=1;
  RooAbsPdf *pdf;
  RooCurve *nomBkgCurve;
  int bestcol= -1;
  for (int icat=0;icat<catIndex->numTypes();icat++){
    int col;
    if (icat<=6) col=color[icat];
    else {col=kBlack; style++;}
    catIndex->setIndex(icat);
    pdfs->getCurrentPdf()->fitTo(*data,RooFit::Minos(0),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE));	 //FIXME
    pdfs->getCurrentPdf()->plotOn(plot,LineColor(col),LineStyle(style));//,RooFit::NormRange("fitdata_1,fitdata_2"));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==icat) {
    ext=" (Best Fit Pdf) ";
    pdf= pdfs->getCurrentPdf();
    nomBkgCurve = (RooCurve*)plot->getObject(plot->numItems()-1);
    bestcol = col;
    }
    leg->AddEntry(pdfLeg,Form("%s%s",pdfs->getCurrentPdf()->GetName(),ext.c_str()),"L");
  }
  plot->SetTitle(Form("Category %s",flashggCats_[cat].c_str()));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->Draw();
  leg->Draw("same");
  CMS_lumi( canv, 0, 0);
  ///start extra bit for ratio plot///
  TH1D *hbplottmp = (TH1D*) pdf->createHistogram("hbplottmp",*mass,Binning(mgg_high-mgg_low,mgg_low,mgg_high));
  hbplottmp->Scale(plotdata->Integral());
  hbplottmp->Draw("same");
  int npoints = plotdata->GetN();
  double xtmp,ytmp;//
  int point =0;
  TGraphAsymmErrors *hdatasub = new TGraphAsymmErrors(npoints);
  //hdatasub->SetMarkerSize(defmarkersize);
  for (int ipoint=0; ipoint<npoints; ++ipoint) {
  //double bkgval = hbplottmp->GetBinContent(ipoint+1);
  plotdata->GetPoint(ipoint, xtmp,ytmp);
  double bkgval = nomBkgCurve->interpolate(xtmp);
  if (BLIND) {
   if ((xtmp > 120 ) && ( xtmp < 130) ) continue;
  }
  std::cout << "[INFO] plotdata->Integral() " <<  plotdata->Integral() << " ( bins " << npoints  << ") hbkgplots[i]->Integral() " << hbplottmp->Integral() << " (bins " << hbplottmp->GetNbinsX() << std::endl;
 double errhi = plotdata->GetErrorYhigh(ipoint);
 double errlow = plotdata->GetErrorYlow(ipoint);
       
 //std::cout << "[INFO]  Channel " << name  << " errhi " << errhi << " errlow " << errlow  << std::endl;
 std::cout << "[INFO] Channel  " << name << " setting point " << point <<" : xtmp "<< xtmp << "  ytmp " << ytmp << " bkgval  " << bkgval << " ytmp-bkgval " << ytmp-bkgval << std::endl;
 bool drawZeroBins_ =1;
 if (!drawZeroBins_) if(fabs(ytmp)<1e-5) continue; 
 hdatasub->SetPoint(point,xtmp,ytmp-bkgval);
 hdatasub->SetPointError(point,0.,0.,errlow,errhi );
 point++;
  } 
  pad2->cd();
  TH1 *hdummy = new TH1D("hdummyweight","",mgg_high-mgg_low,mgg_low,mgg_high);
  hdummy->SetMaximum(hdatasub->GetHistogram()->GetMaximum()+1);
  hdummy->SetMinimum(hdatasub->GetHistogram()->GetMinimum()-1);
  hdummy->GetYaxis()->SetTitle("data - best fit PDF");
  hdummy->GetYaxis()->SetTitleSize(0.12);
  hdummy->GetXaxis()->SetTitle("m_{ll#gamma} (GeV)");
  hdummy->GetXaxis()->SetTitleSize(0.12);
  hdummy->Draw("HIST");
  hdummy->GetYaxis()->SetNdivisions(808);

  TLine *line3 = new TLine(mgg_low,0.,mgg_high,0.);
  line3->SetLineColor(bestcol);
  //line3->SetLineStyle(kDashed);
  line3->SetLineWidth(5.0);
  line3->Draw();
  hdatasub->Draw("PESAME");
  // enf extra bit for ratio plot///
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  catIndex->setIndex(currentIndex);
  delete canv;
}

void plot(RooRealVar *mass, map<string,RooAbsPdf*> pdfs, RooDataSet *data, string name, vector<string> flashggCats_, int cat, int bestFitPdf=-1){
  
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
  TCanvas *canv = new TCanvas();
  TLegend *leg = new TLegend(0.6,0.65,0.88,0.88);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1",mgg_low,120);
  mass->setRange("unblindReg_2",130,mgg_high);
  if (BLIND) {
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(mgg_high-mgg_low),Invisible());
  }
  else data->plotOn(plot,Binning(mgg_high-mgg_low));

  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
	if(flashggCats_.size() >0){
  leg->AddEntry(datLeg,Form("Data - %s",flashggCats_[cat].c_str()),"LEP");
	} else {
  leg->AddEntry(datLeg,Form("Data - %d",cat),"LEP");
	}
  int i=0;
  int style=1;
  for (map<string,RooAbsPdf*>::iterator it=pdfs.begin(); it!=pdfs.end(); it++){
    int col;
    if (i<=6) col=color[i];
    else {col=kBlack; style++;}
    cout << "[INFO] Plotting " << it->first << " with color " << col << " and PDF at " << it->second << endl;
    it->second->plotOn(plot,LineColor(col),LineStyle(style));//,RooFit::NormRange("fitdata_1,fitdata_2"));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==i) ext=" (Best Fit Pdf) ";
    leg->AddEntry(pdfLeg,Form("%s%s",it->first.c_str(),ext.c_str()),"L");
    i++;
  }
  plot->SetTitle(Form(" %s",flashggCats_[cat].c_str()));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->Draw();
  leg->Draw();
  CMS_lumi( canv, 0, 0);
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  delete canv;
}

void transferMacros(TFile *inFile, TFile *outFile){
  
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())){
    if (string(key->ReadObj()->ClassName())=="TMacro") {
      //cout << key->ReadObj()->ClassName() << " : " << key->GetName() << endl;
      TMacro *macro = (TMacro*)inFile->Get(key->GetName());
      outFile->cd();
      macro->Write();
    }
  }
}

int getBestFitFunction(RooMultiPdf *bkg, RooDataSet *data, RooCategory *cat, bool silent=false){
	double global_minNll = 1E10;
	int best_index = 0;
	int number_of_indeces = cat->numTypes();
	string typePrefix = ""; // 添加typePrefix变量声明
	int order = 0; // 添加order变量声明
		
	RooArgSet snap,clean;
	RooArgSet *params = bkg->getParameters((const RooArgSet*)0);
	params->remove(*cat);
	params->snapshot(snap);
	params->snapshot(clean);
	if (!silent) {
		//params->Print("V");
	}
	
	//bkg->setDirtyInhibit(1);
	//RooAbsReal *nllm = bkg->createNLL(*data);
	//RooMinimizer minim(*nllm);
	//minim.setStrategy(1);
	
	for (int id=0;id<number_of_indeces;id++){		
		params->assignValueOnly(clean);
		cat->setIndex(id);

		//RooAbsReal *nllm = bkg->getCurrentPdf()->createNLL(*data);

		if (!silent) {
			/*
			std::cout << "BEFORE  MAKING FIT" << std::endl;
			params->Print("V");
			std::cout << "-----------------------" << std::endl;		
			*/
		}
		
		//minim.minimize("Minuit2","minimize");
		double minNll=0; //(nllm->getVal())+bkg->getCorrection();
		int fitStatus=1;		
		runFit(bkg->getCurrentPdf(),data,&minNll,&fitStatus,/*max iterations*/3,/*doBlind*/BLIND_FIT, 
      /*cat*/cat->GetName(),typePrefix,order,nBinsForMass,bkg->getCurrentPdf());
		// Add the penalty

		minNll=minNll+bkg->getCorrection();

		if (!silent) {
			std::cout << "[INFO] AFTER FITTING" << std::endl;
			std::cout << "[INFO] Function was " << bkg->getCurrentPdf()->GetName() <<std::endl;
			std::cout << "[INFO] Correction Applied is " << bkg->getCorrection() <<std::endl;
			std::cout << "[INFO] NLL + c = " <<  minNll << std::endl;
			std::cout << "-----------------------" << std::endl;
		}
			
		if (minNll < global_minNll){
        		global_minNll = minNll;
			snap.assignValueOnly(*params);
        		best_index=id;
		}
	}
    	cat->setIndex(best_index);
	params->assignValueOnly(snap);
	
	if (!silent) {
		std::cout << "[INFO] Best fit Function -- " << bkg->getCurrentPdf()->GetName() << " " << cat->getIndex() <<std::endl;
		//bkg->getCurrentPdf()->getParameters(*data)->Print("v");
	}
	return best_index;
}

/**
 * Calculates the Kolmogorov-Smirnov test probability for a PDF fitted to data
 * @param mass - The mass variable
 * @param pdf - The PDF to test
 * @param data - The dataset to test against
 * @param name - Base name for output files
 * @return KS test probability
 */
double getKSProb(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, string name) {
  // Calculate bin density based on data entries per GeV
  int mgg_low = mass->getMin(), mgg_high = mass->getMax();
  std::cout << "[INFO] KS test mass range: " << mgg_low << " - " << mgg_high << std::endl;

  // Determine bin density based on event statistics
  double binForGeV = 4;  
  int nBin = (mgg_high - mgg_low) * binForGeV;

  // Create histograms for data and PDF
  TH1F *dataHist = NULL;
  
  // Handle blinded or unblinded data
  if (BLIND) {
    // Create separate datasets for regions outside blind region
    RooDataSet *lowSideband = (RooDataSet*)data->reduce(Form("CMS_hgg_mass < %f", blind_low));
    RooDataSet *highSideband = (RooDataSet*)data->reduce(Form("CMS_hgg_mass > %f", blind_high));
    
    // Convert to histograms
    TH1F *dataHistLow = (TH1F*)lowSideband->createHistogram("dataHistLow", *mass, 
                         RooFit::Binning((blind_low-mgg_low) * binForGeV, mgg_low, blind_low));
    TH1F *dataHistHigh = (TH1F*)highSideband->createHistogram("dataHistHigh", *mass, 
                          RooFit::Binning((mgg_high-blind_high) * binForGeV, blind_high, mgg_high));
    
    // Combine into one histogram
    dataHist = new TH1F("dataHist", "Data Histogram", nBin, mgg_low, mgg_high);
    
    // Fill with data from sidebands
    for (int i = 1; i <= dataHistLow->GetNbinsX(); i++) {
      dataHist->SetBinContent(i, dataHistLow->GetBinContent(i));
      dataHist->SetBinError(i, dataHistLow->GetBinError(i));
    }
    
    int offsetBin = int((blind_high - mgg_low) * binForGeV);
    for (int i = 1; i <= dataHistHigh->GetNbinsX(); i++) {
      dataHist->SetBinContent(i + offsetBin, dataHistHigh->GetBinContent(i));
      dataHist->SetBinError(i + offsetBin, dataHistHigh->GetBinError(i));
    }
    
    delete dataHistLow;
    delete dataHistHigh;
    delete lowSideband;
    delete highSideband;
  } else {
    // Use all data if not blinded
    dataHist = (TH1F*)data->createHistogram("dataHist", *mass, RooFit::Binning(nBin, mgg_low, mgg_high));
  }
  
  // Generate histogram from PDF
  TH1F *pdfHist = (TH1F*)pdf->createHistogram("pdfHist", *mass, RooFit::Binning(nBin, mgg_low, mgg_high));
  // Set bin errors to zero
  if (pdfHist) {
    for (int i = 0; i <= pdfHist->GetNbinsX() + 1; ++i) { // Include underflow/overflow bins just in case
      pdfHist->SetBinError(i, 0.0);
    }
  }
  
  // Set signal region to zero if blinded
  if (BLIND) {
    int lowBin = pdfHist->FindBin(blind_low);
    int highBin = pdfHist->FindBin(blind_high);
    for (int i = lowBin; i < highBin; i++) {
      pdfHist->SetBinContent(i, 0);
      dataHist->SetBinContent(i, 0);
    }
  }

  // Normalize histograms
  if (pdfHist->Integral() > 0) pdfHist->Scale(1.0/pdfHist->Integral());
  if (dataHist->Integral() > 0) dataHist->Scale(1.0/dataHist->Integral());
  
  // Get cumulative distributions
  TH1F *dataCum = (TH1F*)dataHist->GetCumulative();
  TH1F *pdfCum = (TH1F*)pdfHist->GetCumulative();

  // Normalize cumulative distributions
  if (dataCum->Integral() > 0) dataCum->Scale(1.0/dataCum->GetBinContent(dataCum->GetNbinsX()));
  if (pdfCum->Integral() > 0) pdfCum->Scale(1.0/pdfCum->GetBinContent(pdfCum->GetNbinsX()));
  
  // Calculate KS test probability using the histograms directly
  double ksProb = dataHist->KolmogorovTest(pdfHist); // Use TH1::KolmogorovTest

  // Create plot
  TCanvas *canKS = new TCanvas("canKS", "canKS", 800, 800);
  
  // Define uniform text size
  Float_t axisLabelSize = 0.06;   // Axis label size
  Float_t axisTitleSize = 0.07;   // Axis title size
  Float_t axisTitleOffset = 0.9;  // Axis title offset
  Float_t legendTextSize = 0.06;  // Legend text size
  Float_t latexTextSize = 0.06;   // Other text annotation size
  
  // Create pad layout
  double padHeight = 0.4;
  TPad *pad1 = new TPad("pad1","pad1",0,padHeight,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,padHeight);
  double ratio = (1-padHeight)/padHeight;
  pad1->SetBottomMargin(0.001);
  pad1->SetLeftMargin(0.12);
  pad2->SetTopMargin(0.001);
  pad2->SetLeftMargin(0.12);
  pad2->SetBottomMargin(0.25);
  pad1->Draw();
  pad2->Draw();
  
  // Top pad: draw distributions
  pad1->cd();
  dataHist->SetLineColor(kBlack);
  dataHist->SetMarkerStyle(20);
  dataHist->SetMarkerColor(kBlack);
  dataHist->SetTitle("");
  dataHist->GetXaxis()->SetTitle(""); // Remove X axis title from top plot
  dataHist->GetYaxis()->SetTitle("Normalized Events");
  
  // Set axis text size for top pad
  dataHist->GetYaxis()->SetTitleSize(axisTitleSize);
  dataHist->GetYaxis()->SetLabelSize(axisLabelSize);
  dataHist->GetYaxis()->SetTitleOffset(axisTitleOffset);
  dataHist->GetXaxis()->SetLabelSize(0); // Hide X axis labels on top pad
  
  gPad->SetLeftMargin(0.12);
  gStyle->SetOptTitle(0);

  // Ensure Y axis uses appropriate format
  double maxY = dataHist->GetMaximum();
  if (maxY > 0 && maxY < 0.1) {
    dataHist->GetYaxis()->SetMaxDigits(2);
    dataHist->GetYaxis()->CenterTitle(true);
  }
  
  dataHist->Draw("E");
  
  pdfHist->SetLineColor(kBlue);
  pdfHist->SetLineWidth(2);
  pdfHist->Draw("HIST SAME");
  
  // Add legend with uniform text size
  TLegend *leg = new TLegend(0.6, 0.7, 0.89, 0.89);
  leg->SetTextSize(legendTextSize);
  leg->AddEntry(dataHist, "Data", "EP");
  leg->AddEntry(pdfHist, "PDF", "L");
  leg->AddEntry((TObject*)0, Form("KS Prob = %.4f", ksProb), "");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  
  // Add blind region (if needed)
  if (BLIND) {
    TBox *blindBox = new TBox(blind_low, 0, blind_high, dataHist->GetMaximum()*1.1);
    blindBox->SetFillColorAlpha(kGray, 1);
    blindBox->SetFillStyle(3004);
    blindBox->Draw("same");
    
    TLegend *blindLeg = new TLegend(0.6, 0.6, 0.89, 0.7);
    blindLeg->SetTextSize(legendTextSize);
    blindLeg->SetBorderSize(0);
    blindLeg->SetFillStyle(0);
    blindLeg->AddEntry(blindBox, "Blinded Region", "f");
    blindLeg->Draw();
  }
  
  // Bottom pad: draw cumulative distributions
  pad2->cd();
  dataCum->SetLineColor(kBlack);
  dataCum->SetMarkerStyle(20);
  dataCum->GetXaxis()->SetTitle("m_{ll#gamma} (GeV)");
  dataCum->GetYaxis()->SetTitle("Cum. Prob.");
  dataCum->GetYaxis()->SetRangeUser(0, 1.1);
  
  // Set axis text size for bottom pad, ensure consistency with top pad
  dataCum->GetXaxis()->SetTitleSize(axisTitleSize*ratio);
  dataCum->GetYaxis()->SetTitleSize(axisTitleSize*ratio);
  dataCum->GetXaxis()->SetLabelSize(axisLabelSize*ratio);
  dataCum->GetYaxis()->SetLabelSize(axisLabelSize*ratio);
  dataCum->GetYaxis()->SetTitleOffset(axisTitleOffset*ratio);
  dataCum->GetXaxis()->SetTitleOffset(1.0);
  
  dataCum->Draw("E");
  
  pdfCum->SetLineColor(kBlue);
  pdfCum->SetLineWidth(2);
  pdfCum->Draw("HIST SAME");
  
  // Draw maximum KS distance
  double maxDiff = 0;
  int maxBin = 0;
  for (int i = 1; i <= dataCum->GetNbinsX(); i++) {
    double diff = fabs(dataCum->GetBinContent(i) - pdfCum->GetBinContent(i));
    if (diff > maxDiff) {
      maxDiff = diff;
      maxBin = i;
    }
  }
  
  // Add arrow to show KS distance
  double xAtMax = dataCum->GetBinCenter(maxBin);
  double y1 = dataCum->GetBinContent(maxBin);
  double y2 = pdfCum->GetBinContent(maxBin);
  
  TArrow *arrow = new TArrow(xAtMax, y1, xAtMax, y2, 0.02, "|>");
  arrow->SetLineColor(kRed);
  arrow->SetLineWidth(2);
  arrow->Draw();
  
  // Add label with uniform text size
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->SetTextSize(latexTextSize*ratio);
  lat->DrawLatex(0.15, 0.85, Form("KS statistic = %.4f", maxDiff));
  
  // Add blind region to cumulative plot (if needed)
  if (BLIND) {
    TBox *blindBoxCum = new TBox(blind_low, 0, blind_high, 1.1);
    blindBoxCum->SetFillColorAlpha(kGray, 1);
    blindBoxCum->SetFillStyle(3004);
    blindBoxCum->Draw("same");
  }
  
  // Save image
  canKS->SaveAs(Form("%s_KSTest.pdf", name.c_str()));
  
  // Clean up
  delete canKS;
  delete dataHist;
  delete pdfHist;
  delete dataCum;
  delete pdfCum;
  delete arrow;
  delete lat;
  delete leg;
  
  return ksProb;
}

/**
 * Simplified version of the spurious signal test function
 * @param cat        - Category name (e.g. "VBF0")
 * @param funcType   - Type of function to use (e.g. "BernsteinStepxGau")
 * @param order      - Order of the function
 * @param sig        - Signal strength multiplier
 * @param runPeriod  - Run period ("run2" or "run3")
 * @param outDir     - Directory for output files
 * @return true if test passes, false otherwise
 */
bool SpurialSignalTest(RooRealVar *input_mass, TString cat, TString funcType, int order, int sig, TString runPeriod, TString outDir) {
  std::cout << "[INFO] Running spurious signal test for category " << cat << " with function " << funcType << " order " << order << std::endl;
  
  // Set mass range and binning
  // double mgg_low = (cat == "VBF1" || cat == "VBF2" || cat == "VBF3") ? 105 : 100;
  double mgg_low = input_mass->getMin(), mgg_high = input_mass->getMax();
  int bin_times = 4;
  double bin_size = (mgg_high - mgg_low) * bin_times;
  
  // Output file path
  TString outFilePath = Form("%s/%s_%s_%s%d_%dxsig", outDir.Data(),
                           runPeriod.Data(), cat.Data(), funcType.Data(), order, sig);
  ofstream output(outFilePath + ".txt", ofstream::app);
  
  // Helper function: load histogram from file
  auto loadHistogram = [](TFile* file, const TString& name, const TString& cloneName) -> TH1F* {
    if (!file || !file->GetListOfKeys()->Contains(name)) return nullptr;
    TH1F* h = (TH1F*)file->Get(name);
    return h ? (TH1F*)h->Clone(cloneName) : nullptr;
  };
  
  // Load background and signal from template file
  TFile* fbkg = TFile::Open(Form("/eos/user/j/jiehan/root/input_finalfit/templates/template_%s.root", runPeriod.Data()));
  if (!fbkg || fbkg->IsZombie()) {
    std::cerr << "[ERROR] Failed to open template file" << std::endl;
    output.close();
    return false;
  }
  
  // Load background histogram
  TString bkgHistName = Form("bkg_%s_%s", runPeriod.Data(), cat.Data());
  TH1F* hbkg = loadHistogram(fbkg, bkgHistName, Form("%s_clone", bkgHistName.Data()));
  if (!hbkg) {
    std::cerr << "[ERROR] Background template not found: " << bkgHistName << std::endl;
    fbkg->Close();
    output.close();
    return false;
  }
  
  // Load data histograms (full range and sideband)
  TString frHistName = Form("data_full_%s_%s", runPeriod.Data(), cat.Data());
  TString sbHistName = Form("data_%s_%s", runPeriod.Data(), cat.Data());
  TH1F* hfr = loadHistogram(fbkg, frHistName, Form("%s_clone", frHistName.Data()));
  TH1F* hsb = loadHistogram(fbkg, sbHistName, Form("%s_clone", sbHistName.Data()));
  
  if (!hfr || !hsb) {
    std::cerr << "[ERROR] Data templates not found" << std::endl;
    fbkg->Close();
    output.close();
    return false;
  }
  
  // Get integrals
  double dataevents = hbkg->Integral();
  double mcsbevents = hbkg->Integral(0, (120 - mgg_low) * bin_times) + 
                      hbkg->Integral(bin_times * (mgg_high - 130), bin_size);
  double sbevents = hsb->Integral();
  
  // Scale background to match sideband data
  if (mcsbevents > 0) {
    hbkg->Scale(sbevents / mcsbevents);
  } else {
    std::cerr << "[ERROR] mcsbevents is zero, cannot scale background" << std::endl;
    fbkg->Close();
    output.close();
    return false;
  }
  
  // Open signal template file and get signal PDF
  TString signalFileName = Form("/eos/user/j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Signal/outdir_combinedPDFs/CMS-HGG_combinedPDFs_%s.root", cat.Data());
  TFile* signalFile = TFile::Open(signalFileName);
  if (!signalFile) {
    std::cerr << "[ERROR] Failed to open signal template file" << std::endl;
    fbkg->Close();
    output.close();
    return false;
  }
  
  // Get workspace and variables
  RooWorkspace* ws = (RooWorkspace*)signalFile->Get("wsig_13TeV");
  if (!ws) {
    std::cerr << "[ERROR] Workspace not found" << std::endl;
    signalFile->Close();
    fbkg->Close();
    output.close();
    return false;
  }
  
  RooRealVar *mass = ws->var("CMS_hgg_mass");
  RooRealVar *mh = ws->var("MH");
  if (!mass || !mh) {
    std::cerr << "[ERROR] Mass variables not found" << std::endl;
    signalFile->Close();
    fbkg->Close();
    output.close();
    return false;
  }
  
  // Set Higgs mass
  mh->setVal(125);
  mass->setRange(mgg_low, mgg_high);

  std::cout << "[INFO] Spurious Signal Tsting: mass range " << mass->getMin() << " - " << mass->getMax() << std::endl;
  
  // Get signal PDF
  RooAbsPdf* signalPdf = ws->pdf(Form("combinedSigPdf_%s", cat.Data()));
  if (!signalPdf) {
    std::cerr << "[ERROR] Signal PDF not found" << std::endl;
    signalFile->Close();
    fbkg->Close();
    output.close();
    return false;
  }
  
  // Fix signal PDF parameters
  RooArgSet* signalParams = signalPdf->getParameters(RooArgSet(*mass));
  TIterator* iter = signalParams->createIterator();
  RooRealVar* param;
  while ((param = (RooRealVar*)iter->Next())) {
    if (!param->isConstant()) param->setConstant(true);
  }
  delete iter;
  
  // Get signal yield
  TString signalHistName = Form("sig_%s_%s", runPeriod.Data(), cat.Data());
  TH1F* hsig = loadHistogram(fbkg, signalHistName, Form("%s_clone", signalHistName.Data()));
  if (!hsig) {
    std::cerr << "[ERROR] Signal template not found: " << signalHistName << std::endl;
    signalFile->Close();
    fbkg->Close();
    output.close();
    return false;
  }
  
  double sigevents = hsig->Integral();
  
  // Create RooDataHist objects
  RooDataHist* dbkg = new RooDataHist("bkg_mc", "dataset with x", *mass, hbkg);
  RooDataHist* dsb = new RooDataHist("data_sb", "dataset with x", *mass, hsb);
  RooDataHist* dfr = new RooDataHist("data_fr", "dataset with x", *mass, hfr);
  
  // Set signal and background normalization variables
  RooRealVar nsig("nsig", "nsig", sigevents, -100 * sigevents, 100 * sigevents);
  RooRealVar nbkg("nbkg", "nbkg", dataevents, 0.01 * dataevents, 2 * dataevents);
  
  // Create PDF model builder
  PdfModelBuilder pdfsModel;
  pdfsModel.setObsVar(mass);
  
  // Set variables to track fit results
  TString status = "Pass";
  double chi2 = 0, prob = 0, nll = 0;
  double dmc = 0, dss = 0, ss = 0, ss_mc = 0;
  double tot_err = 0, ss_cor = 0, delta = 0;
  int fit_status = 0;
  
  string typePrefix;
  // Create background model
  RooAbsPdf* bkg_model = getPdf(pdfsModel, funcType.Data(), order, &typePrefix, 
                              Form("sstest_pdf_%s_%s%d", cat.Data(), funcType.Data(), order), cat.Data());
  if (!bkg_model) {
    std::cerr << "[ERROR] Failed to create background PDF" << std::endl;
    signalFile->Close();
    fbkg->Close();
    output.close();
    return false;
  }
  
  // Fit background model to background template
  RooFitResult* bkg_model_fit = nullptr;
  bkg_model_fit = bkg_model->fitTo(*dbkg, Save(1), Minimizer("Minuit2", "minimize"), 
                                   SumW2Error(kTRUE), PrintLevel(-1));
  
  // Plot background fit
  RooPlot* frame_bkg = mass->frame(Title(Form("Background %s Fit", funcType.Data())));
  dbkg->plotOn(frame_bkg);
  bkg_model->plotOn(frame_bkg);
  bkg_model->paramOn(frame_bkg, Layout(0.34, 0.96, 0.89), Format("NEA", AutoPrecision(1)));
  bkg_model->SetName(funcType);
  
  nll = bkg_model_fit->minNll();
  int bkg_npars = bkg_model_fit->floatParsFinal().getSize();
  int bkg_ndof = bin_size - bkg_npars;
  
  output << "\t" << funcType.Data() << "\tbkg:\tnpars = " << bkg_npars 
         << " \tchi^2 = " << frame_bkg->chiSquare(bkg_npars) 
         << "\tprob = " << TMath::Prob(frame_bkg->chiSquare(bkg_npars) * bkg_ndof, bkg_ndof) 
         << "\tnll: " << nll << "\tstatus = " << bkg_model_fit->status() << std::endl;
  
  // Create pseudo data: inject signal
  TH1F* hdata = (TH1F*)hbkg->Clone("hdata");
  hsig->Reset();
  
  // Normalize signal PDF and fill histogram
  RooAbsReal* normIntegral = signalPdf->createIntegral(*mass);
  if (!normIntegral || normIntegral->getVal() <= 0) {
    std::cerr << "[ERROR] Signal PDF normalization invalid" << std::endl;
    delete bkg_model;
    signalFile->Close();
    fbkg->Close();
    output.close();
    return false;
  }
  
  double normVal = normIntegral->getVal();
  
  // Fill signal histogram
  for (int i = 1; i <= hsig->GetNbinsX(); i++) {
    double x = hsig->GetBinCenter(i);
    mass->setVal(x);
    double val = signalPdf->getVal() * sigevents * hsig->GetBinWidth(i) / normVal;
    if (val < 0 || val > 100) val = 0;
    hsig->SetBinContent(i, val);
  }
  
  // Add signal to background
  hdata->Add(hsig, sig);
  for (int i=0; i<=hdata->GetNbinsX(); i++) {
    double x = hdata->GetBinCenter(i);
    double content = hdata->GetBinContent(i);
    // std::cout << "[INFO] Bin " << i << ": x = " << x << ", content = " << content << std::endl;
  }
  RooDataHist* ddata = new RooDataHist("data_bin", "dataset with x", *mass, hdata);
  
  // Create combined model
  RooAddPdf* model = new RooAddPdf("model", "model", RooArgList(*signalPdf, *bkg_model), 
                                 RooArgList(nsig, nbkg));
  
  // Set canvas and plot frames
  TCanvas* canv = new TCanvas();
  RooPlot* frame_data = mass->frame();
  RooPlot* frame_data_trash = mass->frame();
  
  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  TPad* pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.35);
  pad1->SetBottomMargin(0.18);
  pad2->SetBottomMargin(0.25);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  
  // Fit combined model to data
  RooFitResult* model_fit = nullptr;
  int MaxTries = 3;

  // runFit(model, ddata, &nll, &fit_status, /*max iterations*/3, /*doBlind*/false, /*cat*/cat.Data(), typePrefix, order, bin_size, model);

  for (int i = 0; i < MaxTries; i++) {
    model_fit = model->fitTo(*ddata, Save(1), Minimizer("Minuit2", "minimize"), 
                            SumW2Error(kTRUE), PrintLevel(-1));
    // if (model_fit && model_fit->status() == 0) break;
    if (i == MaxTries - 1) {
      std::cerr << "[ERROR] Fit failed after " << MaxTries << " attempts" << std::endl;
    }
  }
  fit_status = model_fit->status();
  
  // Extract signal strength and error
  ss_mc = nsig.getVal();
  dmc = nsig.getError();
  
  // Plot data and fit
  ddata->plotOn(frame_data, Name("data"), DataError(RooAbsData::SumW2));
  RooHist* plotdata = (RooHist*)frame_data->getObject(frame_data->numItems() - 1);
  
  ddata->plotOn(frame_data_trash, Name("data"), DataError(RooAbsData::SumW2));
  model->plotOn(frame_data_trash, Name("fit"));
  
  int data_npars = model_fit ? model_fit->floatParsFinal().getSize() : 1;
  int data_ndof = bin_size - data_npars;
  chi2 = frame_data_trash->chiSquare(data_npars);
  prob = TMath::Prob(chi2 * data_ndof, data_ndof);

  output << "\t" << funcType.Data() << "\tdata(MC):\tnpars = " << data_npars 
         << " \tchi^2 = " << chi2 << "\tprob = " << prob 
         << "\tnll: " << nll << "\tstatus = " << fit_status
         << std::endl;
  
  // runFit(model, ddata, &nll, &fit_status, /*max iterations*/3, /*doBlind*/false, /*cat*/cat.Data(), typePrefix, order, bin_size, model);

  // Refit using unweighted errors
  for (int i = 0; i < MaxTries; i++) {
    model_fit = model->fitTo(*ddata, Save(1), Minimizer("Minuit2", "minimize"), 
                           SumW2Error(kFALSE), PrintLevel(-1));
    // if (model_fit && model_fit->status() == 0) break;
    if (i == MaxTries - 1) {
      std::cerr << "[ERROR] Fit failed after " << MaxTries << " attempts" << std::endl;
    }
  }
  fit_status = model_fit->status();

  // Get signal strength and error
  ss = nsig.getVal();
  dss = nsig.getError();
  tot_err = sqrt(dss * dss + ss * ss);
  delta = fabs(ss) - 2 * dmc;
  ss_cor = (delta < 0) ? 0 : ((ss > 0) ? delta : -delta);
  
  // Test failure criteria
  if (delta > 0.2 * dss) status = "Fail";

  output << "\t" << funcType.Data() << "\tdata(Po):\tnpars = " << data_npars 
         << " \tchi^2 = " << chi2 << "\tprob = " << prob 
         << "\tnll: " << nll << "\tstatus = " << fit_status
         << "\tsig: " << ss_mc << " +/- " << dmc 
         << "\tsig (unweighted): " << ss << " +/- " << dss 
         << "\tdelta: " << delta 
         << "\ttot_err: " << tot_err 
         << "\tss_cor: " << ss_cor 
         << std::endl;
  
  // Plot fit results
  model->plotOn(frame_data, Name("fit"));
  model->plotOn(frame_data, Name("background"), Components(bkg_model->GetName()), 
              LineStyle(ELineStyle::kDashed), LineColor(kGreen));
  RooCurve* nomBkgCurve = (RooCurve*)frame_data->getObject(frame_data->numItems() - 1);
  
  // Create legend
  TLegend* leg = new TLegend(0.6, 0.65, 0.88, 0.88);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->AddEntry(frame_data->findObject("data"), "Data", "ep");
  leg->AddEntry(frame_data->findObject("fit"), "Bkg+Sig", "l");
  leg->AddEntry(frame_data->findObject("background"), "Bkg", "l");
  
  // Set title and style
  model->SetName(Form("%s_model", funcType.Data()));
  frame_data->SetTitle(Form("Pseudo Data (x%d Sig), Prob: %.3f", sig, prob));
  frame_data->SetXTitle("");
  frame_data->SetLabelSize(0.042, "XY");
  frame_data->SetTitleSize(0.056, "Y");
  frame_data->SetTitleOffset(0.75, "Y");
  frame_data->Draw();
  leg->Draw("same");
  
  // Add signal curve
  signalPdf->plotOn(frame_data, Name("signal"),Normalization(ss, RooAbsReal::NumEvent), 
                  LineColor(kRed), LineWidth(4));
  
  // Create residuals panel
  pad2->cd();
  TH1D* hdummy = new TH1D("hdummyweight", "", mgg_high - mgg_low, mgg_low, mgg_high);
  hdummy->SetStats(0);
  hdummy->SetMaximum(5);
  hdummy->SetMinimum(-5);
  hdummy->GetYaxis()->SetTitle("Data - Bkg PDF");
  hdummy->GetYaxis()->SetTitleOffset(0.35);
  hdummy->GetYaxis()->SetTitleSize(0.12);
  hdummy->GetYaxis()->SetLabelSize(0.09);
  hdummy->GetXaxis()->SetTitle("m_{ll#gamma} (GeV)");
  hdummy->GetXaxis()->SetTitleSize(0.12);
  hdummy->GetXaxis()->SetLabelSize(0.09);
  hdummy->Draw("HIST");
  
  // Plot residuals and signal curve
  TGraphAsymmErrors* hdatasub = new TGraphAsymmErrors();
  int point = 0;
  for (int i = 0; i < plotdata->GetN(); i++) {
    double x, y;
    plotdata->GetPoint(i, x, y);
    double bkg = nomBkgCurve->interpolate(x);
    hdatasub->SetPoint(point, x, y - bkg);
    hdatasub->SetPointError(point, 0, 0, plotdata->GetErrorYlow(i), plotdata->GetErrorYhigh(i));
    point++;
  }
  
  hdatasub->SetMarkerStyle(8);
  hdatasub->Draw("PESAME");
  
  // Save image
  canv->SaveAs(Form("%s/pesudo_data_shape_%s_%s_%s%d.png", outDir.Data(),
                      runPeriod.Data(), cat.Data(), funcType.Data(), order));
  
  // Clean up memory
  delete bkg_model;
  delete canv;
  delete hdummy;
  delete hdatasub;
  delete leg;
  signalFile->Close();
  fbkg->Close();
  output.close();
  
  // Return test result
  return (status == "Pass");
}

int main(int argc, char* argv[]){
 
  setTDRStyle();
  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  string year_ = "2016";
  //int year_ = 2017;

  string fileName;
  int ncats;
  int singleCategory;
  int catOffset;
  string datfile;
  string outDir;
  string outfilename;
  bool is2011=false;
  bool verbose=false;
  bool saveMultiPdf=false;
  int isFlashgg_ =1;
  string flashggCatsStr_;
  vector<string> flashggCats_;
  bool isData_ =0;
  double mgg_low, mgg_high;
  int nBinsForMass;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilename,i", po::value<string>(&fileName),                                              "In file name")
    ("ncats,c", po::value<int>(&ncats)->default_value(5),                                       "Number of categories")
    ("singleCat", po::value<int>(&singleCategory)->default_value(-1),                           "Run A single Category")
    ("datfile,d", po::value<string>(&datfile)->default_value("dat/fTest.dat"),                  "Right results to datfile for BiasStudy")
    ("outDir,D", po::value<string>(&outDir)->default_value("plots/fTest"),                      "Out directory for plots")
    ("saveMultiPdf", po::value<string>(&outfilename),         					"Save a MultiPdf model with the appropriate pdfs")
    ("runFtestCheckWithToys", 									"When running the F-test, use toys to calculate pvals (and make plots) ")
    ("is2011",                                                                                  "Run 2011 config")
    ("is2012",                                                                                  "Run 2012 config")
    ("unblind",  									        "Dont blind plots")
    ("blindFit",                                                                               "Blind fits in signal region (120-130 GeV)")
    ("isFlashgg",  po::value<int>(&isFlashgg_)->default_value(1),  								    	        "Use Flashgg output ")
    ("isData",  po::value<bool>(&isData_)->default_value(0),  								    	        "Use Data not MC ")
		("flashggCats,f", po::value<string>(&flashggCatsStr_)->default_value("UntaggedTag_0,UntaggedTag_1,UntaggedTag_2,UntaggedTag_3,UntaggedTag_4,VBFTag_0,VBFTag_1,VBFTag_2,TTHHadronicTag,TTHLeptonicTag,VHHadronicTag,VHTightTag,VHLooseTag,VHEtTag"),       "Flashgg category names to consider")
    ("year", po::value<string>(&year_)->default_value("2016"),       "Dataset year")
    ("catOffset", po::value<int>(&catOffset)->default_value(0),       "Category numbering scheme offset")
    ("mgg_low", po::value<double>(&mgg_low)->default_value(95),                            "Lower bound for mgg")
    ("mgg_high", po::value<double>(&mgg_high)->default_value(170),                                                "Upper bound for mgg")
    ("verbose,v",                                                                               "Run with more output")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  nBinsForMass = 4*(mgg_high-mgg_low);
  if (vm.count("help")) { cout << desc << endl; exit(1); }
  if (vm.count("is2011")) is2011=true;
	if (vm.count("unblind")) BLIND=false;
  if (vm.count("blindFit")) BLIND_FIT=true;
  saveMultiPdf = vm.count("saveMultiPdf");

  if (vm.count("verbose")) verbose=true;
  if (vm.count("runFtestCheckWithToys")) runFtestCheckWithToys=true;

  if (!verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
    gErrorIgnoreLevel=kWarning;
  }
	split(flashggCats_,flashggCatsStr_,boost::is_any_of(","));
  
	int startingCategory=0;
  if (singleCategory >-1){
    ncats=singleCategory+1;	
    startingCategory=singleCategory;
  }
	if (isFlashgg_==1){
	
	ncats= flashggCats_.size();

	}

  if(verbose) std::cout << "[INFO] SaveMultiPdf? " << saveMultiPdf << std::endl;
  TFile *outputfile;
  RooWorkspace *outputws;

  if (saveMultiPdf){
    outputfile = new TFile(outfilename.c_str(),"RECREATE");
    outputws = new RooWorkspace(); outputws->SetName("multipdf");
  }

  system(Form("mkdir -p %s",outDir.c_str()));
  TFile *inFile = TFile::Open(fileName.c_str());
  RooWorkspace *inWS;
	if(isFlashgg_){
		if (isData_){
			inWS = (RooWorkspace*)inFile->Get("tagsDumper/cms_hgg_13TeV");
		} else {
			inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
		}
	} else {
		inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");//FIXME
	}
	if (verbose) std::cout << "[INFO]  inWS open " << inWS << std::endl;
	if (saveMultiPdf){
		transferMacros(inFile,outputfile);

		RooRealVar *intL; 
		RooRealVar *sqrts;

		if (isFlashgg_){
			//intL  = (RooRealVar*)inWS->var("IntLumi");
			intL  = intLumi_;
			sqrts = (RooRealVar*)inWS->var("SqrtS");
			if (!sqrts){ sqrts = new RooRealVar("SqrtS","SqrtS",13); }
		std::cout << "[INFO] got intL and sqrts " << intL << ", " << sqrts << std::endl;


		} else {
			//intL  = (RooRealVar*)inWS->var("IntLumi");
			intL  = intLumi_;
			sqrts = (RooRealVar*)inWS->var("Sqrts");
		}
		outputws->import(*intL);
		outputws->import(*sqrts);
		std::cout << "[INFO] got intL and sqrts " << intL << ", " << sqrts << std::endl;
	}

	vector<string> functionClasses;
	// functionClasses.push_back("Bernstein");
	// functionClasses.push_back("Exponential");
	// functionClasses.push_back("PowerLaw");
	// functionClasses.push_back("Laurent");
  functionClasses.push_back("PowerLawStepxGau");
  functionClasses.push_back("LaurentStepxGau");
  functionClasses.push_back("ExponentialStepxGau");
  functionClasses.push_back("BernsteinStepxGau");
  functionClasses.push_back("ExpModGauss");
  // functionClasses.push_back("AsymGenGauss");
	map<string,string> namingMap;
	namingMap.insert(pair<string,string>("Bernstein","pol"));
	namingMap.insert(pair<string,string>("Exponential","exp"));
	namingMap.insert(pair<string,string>("PowerLaw","pow"));
	namingMap.insert(pair<string,string>("Laurent","lau"));
  namingMap.insert(pair<string,string>("BernsteinStepxGau","bern"));
  namingMap.insert(pair<string,string>("ExponentialStepxGau","exp"));
  namingMap.insert(pair<string,string>("PowerLawStepxGau","pow"));
  namingMap.insert(pair<string,string>("LaurentStepxGau","lau"));
  namingMap.insert(pair<string,string>("ExpModGauss","modgau"));
  namingMap.insert(pair<string,string>("AsymGenGauss","asymgau"));

	// store results here

	FILE *resFile = NULL;
	if  (singleCategory >-1) resFile = fopen(Form("%s/fTestResults_%s.txt",outDir.c_str(),flashggCats_[singleCategory].c_str()),"w");
	else resFile = fopen(Form("%s/fTestResults.txt",outDir.c_str()),"w");
	
	if (!resFile) {
		std::cerr << "[ERROR] Failed to open results file for writing. Check permissions and path: " 
			  << outDir.c_str() << std::endl;
		return 1;
	}
	
	vector<map<string,int> > choices_vec;
	vector<map<string,std::vector<int> > > choices_envelope_vec;
	vector<map<string,RooAbsPdf*> > pdfs_vec;

	PdfModelBuilder pdfsModel;
	RooRealVar *mass = (RooRealVar*)inWS->var("CMS_hgg_mass");
  mass->setRange(mgg_low,mgg_high);
  mass->setBins(nBinsForMass);
	std:: cout << "[INFO] Got mass from ws " << mass << std::endl;
	pdfsModel.setObsVar(mass);
	double upperEnvThreshold = 0.05; // upper threshold on delta(chi2) to include function in envelope (looser than truth function)

	fprintf(resFile,"Truth Model & d.o.f & $\\Delta NLL_{N+1}$ & $p(\\chi^{2}>\\chi^{2}_{(N\\rightarrow N+1)})$ \\\\\n");
	fprintf(resFile,"\\hline\n");

	std::string ext = is2011 ? "7TeV" : "8TeV";
  if( isFlashgg_ ){
    if( year_.find("run2")) { ext = "13TeV"; }
    else if ( year_.find("run3")) { ext = "13.6TeV"; }
    else if ( year_.find("all")) { ext = "13TeV & 13.6TeV"; }
    //else{ ext = "13TeV"; } //FIXME 
    else{ ext = Form("%s_13TeV",year_.c_str()); }
  }
	//if (isFlashgg_) ext = "13TeV";
        //FIXME trying to remove duplicated names for 2016+2017 combination
	//if (isFlashgg_) ext = Form("13TeV_%d",year_);
	for (int cat=startingCategory; cat<ncats; cat++){

		map<string,int> choices;
		map<string,std::vector<int> > choices_envelope;
		map<string,RooAbsPdf*> pdfs;
    map<string,RooAbsPdf*> allPdfs;
		string catname;
		if (isFlashgg_){
			catname = Form("%s",flashggCats_[cat].c_str());
		} else {
			catname = Form("cat%d",cat);
		}
		RooDataSet *dataFull;
		RooDataSet *dataFull0;
		if (isData_) {
    dataFull0 = (RooDataSet*)inWS->data(Form("Data_13TeV_%s",catname.c_str()));
    /*dataFull= (RooDataSet*) dataFull0->emptyClone();
    for (int i =0 ; i < dataFull0->numEntries() ; i++){
    double m = dataFull0->get(i)->getRealValue("CMS_hgg_mass");
    //if (m <(mgg_low+0.01) or m > (mgg_high-0.01)) 

    if (m==mgg_low){
    std::cout << "dataset mass m="<< m << std::endl;
    continue;
    }
    dataFull->add(*dataFull0->get(),1.0);
    }*/
		if (verbose) std::cout << "[INFO] opened data for  "  << Form("Data_%s",catname.c_str()) <<" - " << dataFull <<std::endl;
    }
		else 
    {dataFull0 = (RooDataSet*)inWS->data(Form("data_mass_%s",catname.c_str()));
		if (verbose) std::cout << "[INFO] opened data for  "  << Form("data_mass_%s",catname.c_str()) <<" - " << dataFull <<std::endl;
    }
    dataFull = (RooDataSet*)dataFull0->emptyClone();
    for (int i =0 ; i < dataFull0->numEntries() ; i++){
      double m = dataFull0->get(i)->getRealValue("CMS_hgg_mass");
      if (m <(mgg_low+0.01) or m > (mgg_high-0.01)) continue;
      dataFull->add(*dataFull0->get(),1.0);
    }

		RooDataSet *data;
		//	RooDataHist thisdataBinned(Form("roohist_data_mass_cat%d",cat),"data",*mass,*dataFull);
		//	RooDataSet *data = (RooDataSet*)&thisdataBinned;
    string thisdataBinned_name;

		if ( isFlashgg_){
			thisdataBinned_name =Form("roohist_data_mass_%s",flashggCats_[cat].c_str());
			//	RooDataHist thisdataBinned(Form("roohist_data_mass_cat%d",cat),"data",*mass,*dataFull);
			//	data = (RooDataSet*)&thisdataBinned;
			//		std::cout << "debug " << thisdataBinned.GetName() << std::endl;

			//RooDataSet *data = (RooDataSet*)dataFull;
		} 
    else {
			thisdataBinned_name= Form("roohist_data_mass_cat%d",cat);
			//RooDataSet *data = (RooDataSet*)dataFull;
		}
		RooDataHist thisdataBinned(thisdataBinned_name.c_str(),"data",*mass,*dataFull);
		data = (RooDataSet*)&thisdataBinned;

		RooArgList storedPdfs("store");

		fprintf(resFile,"\\multicolumn{4}{|c|}{\\textbf{Category %d}} \\\\\n",cat);
		fprintf(resFile,"\\hline\n");

		double MinimimNLLSoFar=1e10;
		int simplebestFitPdfIndex = 0;

    // Determine runPeriod based on year parameter
    TString runPeriod = "all";
    if (year_.find("run3") != std::string::npos) {
      runPeriod = "run3";
    }
    else if (year_.find("run2") != std::string::npos) {
      runPeriod = "run2";
    }

		// Standard F-Test to find the truth functions
		for (vector<string>::iterator funcType=functionClasses.begin(); funcType!=functionClasses.end(); funcType++){

			double thisNll=0.; double prevNll=0.; double chi2=0.; double prob=0.; 
			int order=1; int prev_order=0; int cache_order=0;

			RooAbsPdf *prev_pdf=NULL;
			RooAbsPdf *cache_pdf=NULL;
			std::vector<int> pdforders;

      // Run SpurialSignalTest
      TString catName;
      if (isFlashgg_) {
        catName = flashggCats_[cat].c_str();
      } else {
        catName = Form("cat%d", cat);
      }
      // Convert function type and order to format needed for SpurialSignalTest
      TString funcTypeStr = funcType->c_str();
      // Run SpurialSignalTest
      bool passedSSTest = true;

			int counter =0;
      string typePrefix;
      // while (prob<upperEnvThreshold){
			while (prob<upperEnvThreshold && order < 7){
        RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order,&typePrefix,Form("ftest_pdf_%d_%s",(cat+catOffset),ext.c_str()), catname);
        if (!bkgPdf){
          // assume this order is not allowed
          order++;
        }
        else {
          passedSSTest = SpurialSignalTest(mass, catName, funcTypeStr, order, 0, runPeriod, outDir);
          
          if (!passedSSTest) {
            std::cout << "[INFO] " << funcTypeStr << " with order " << order 
                      << " failed SpurialSignalTest, skipping to next order" << std::endl;
            order++;
            continue;
          }
          
          std::cout << "[INFO] " << funcTypeStr << " with order " << order 
                  << " passed SpurialSignalTest, proceeding with F-test" << std::endl;

          int fitStatus = 0;
          bkgPdf->Print();
          runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/3, BLIND_FIT, flashggCats_[cat].c_str(), typePrefix, order, nBinsForMass, bkgPdf);
          if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;       
          chi2 = 2.*(prevNll-thisNll);
          if (chi2<0. && order>1) chi2=0.;
          if (prev_pdf!=NULL){
            prob = getProbabilityFtest(chi2,order-prev_order,prev_pdf,bkgPdf,mass,data
                ,Form("%s/Ftest_from_%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,(cat+catOffset)));
            std::cout << "[INFO]  F-test Prob(chi2>chi2(data)) == " << prob << std::endl;
          } else {
            prob = 0;
          }
          double gofProb=0;
          // otherwise we get it later ...
          if (!saveMultiPdf) plot(mass,bkgPdf,data,Form("%s/%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,(cat+catOffset)),flashggCats_,fitStatus,&gofProb);
          cout << "[INFO]\t " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
          //fprintf(resFile,"%15s && %d && %10.2f && %10.2f && %10.2f \\\\\n",funcType->c_str(),order,thisNll,chi2,prob);
          prevNll=thisNll;
          if (prev_pdf==NULL) {
            cache_pdf=bkgPdf;
            cache_order=order;
          }
          else {
            cache_pdf=prev_pdf;
            cache_order=prev_order;
          }
          prev_order=order;
          prev_pdf=bkgPdf;
          cout << "[INFO] Ftest\t " << *funcType << " " << cache_order << " " << prev_order << endl;
          order++;
        }
        counter++;
      }

			fprintf(resFile,"%15s & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
      cout << "[INFO] Ftest upper limit " << cache_order << " " << prob << endl;
      if (cache_order > 0) {
			  choices.insert(pair<string,int>(*funcType,cache_order));
        pdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),cache_order),cache_pdf));
      }
			int truthOrder = cache_order;

			// Now run loop to determine functions inside envelope
			if (saveMultiPdf){
				chi2=0.;
				thisNll=0.;
				prevNll=0.;
				prob=0.;
				order=1;
				prev_order=0;
				cache_order=0;
				std::cout << "[INFO] Determining Envelope Functions for Family " << *funcType << ", cat " << cat << std::endl;
				std::cout << "[INFO] Upper end Threshold for highest order function " << upperEnvThreshold <<std::endl;

				while (prob<upperEnvThreshold && order < 7){
					RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order,&typePrefix,Form("env_pdf_%d_%s",(cat+catOffset),ext.c_str()), catname);
          cout << "[INFO] get pdf called " << *funcType << " " << order << " " << bkgPdf << endl;
					if (!bkgPdf ){
						// assume this order is not allowed
						order++;
					}
					else {
            passedSSTest = SpurialSignalTest(mass, catName, funcTypeStr, order, 0, runPeriod, outDir);
          
            if (!passedSSTest) {
              std::cout << "[INFO] " << funcTypeStr << " with order " << order 
                        << " failed SpurialSignalTest, skipping to next order" << std::endl;
              order++;
              continue;
            }
            std::cout << "[INFO] SSTest passed for " << funcTypeStr << " with order " << order 
                      << ", proceeding with F-test" << std::endl;

						//RooFitResult *fitRes;
						int fitStatus=0;
						runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/3, BLIND_FIT, flashggCats_[cat].c_str(), typePrefix, order, nBinsForMass, bkgPdf);
						//thisNll = fitRes->minNll();
						if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;
						double myNll = 2.*thisNll;
						chi2 = 2.*(prevNll-thisNll);
						// if (chi2<0. && order>1) chi2=0.;
						prob = TMath::Prob(chi2,order-prev_order);

						cout << "[INFO] \t " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
						prevNll=thisNll;
						cache_order=prev_order;
						cache_pdf=prev_pdf;

						// Calculate goodness of fit for the thing to be included (will use toys for lowstats)!
						double gofProb =0; 
						plot(mass,bkgPdf,data,Form("%s/%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,(cat+catOffset)),flashggCats_,fitStatus,&gofProb);
                        
            // Calculate KS test probability
            double ksProb = getKSProb(mass, bkgPdf, dataFull, Form("%s/%s%d_cat%d",outDir.c_str(),funcType->c_str(),order,(cat+catOffset)));
            cout << "[INFO] \t KS test probability = " << ksProb << endl;

						if ((prob < upperEnvThreshold) ) { // Looser requirements for the envelope

							// if (gofProb > 0.01 || order == truthOrder ) {  // Good looking fit or one of our regular truth functions
              if (gofProb > 0.1 && ksProb > 0.2) {  // Only good looking fit with acceptable KS test

								std::cout << "[INFO] Adding to Envelope " << bkgPdf->GetName() << " "<< gofProb 
									<< " 2xNLL + c is " << myNll + bkgPdf->getVariables()->getSize() << " truth order " << truthOrder << std::endl;
								allPdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),order),bkgPdf));
								storedPdfs.add(*bkgPdf);
								pdforders.push_back(order);

								// Keep track but we shall redo this later
								if ((myNll + bkgPdf->getVariables()->getSize()) < MinimimNLLSoFar) {
									simplebestFitPdfIndex = storedPdfs.getSize()-1;
									MinimimNLLSoFar = myNll + bkgPdf->getVariables()->getSize();
								}
							}
						}

						prev_order=order;
						prev_pdf=bkgPdf;
						order++;
            cout << "[INFO] Ftest envelope " << *funcType << " " << order << " " << prob << endl;
					}
				}

				fprintf(resFile,"%15s & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
				choices_envelope.insert(pair<string,std::vector<int> >(*funcType,pdforders));
			}
		}

		fprintf(resFile,"\\hline\n");
		choices_vec.push_back(choices);
		choices_envelope_vec.push_back(choices_envelope);
		pdfs_vec.push_back(pdfs);
    for (map<string,RooAbsPdf*>::iterator it=pdfs.begin(); it!=pdfs.end(); it++){
      std::cout << "[INFO]  " << it->first << " " << it->second << std::endl;
    }
    // Some categories don't have any functions in the envelope
		plot(mass,pdfs,data,Form("%s/truths_cat%d",outDir.c_str(),(cat+catOffset)),flashggCats_,cat);

		if (saveMultiPdf){
			// Put selectedModels into a MultiPdf
			string catindexname;
			string catname;
			if (isFlashgg_){
				catindexname = Form("pdfindex_%s_%s",flashggCats_[cat].c_str(),ext.c_str());
				catname = Form("%s",flashggCats_[cat].c_str());
			} else {
				catindexname = Form("pdfindex_%d_%s",(cat+catOffset),ext.c_str());
				catname = Form("cat%d",(cat+catOffset));
			}
			RooCategory catIndex(catindexname.c_str(),"c");
			RooMultiPdf *pdf = new RooMultiPdf(Form("CMS_hgg_%s_%s_bkgshape",catname.c_str(),ext.c_str()),"all pdfs",catIndex,storedPdfs);
			//RooRealVar nBackground(Form("CMS_hgg_%s_%s_bkgshape_norm",catname.c_str(),ext.c_str()),"nbkg",data->sumEntries(),0,10E8);
			RooRealVar nBackground(Form("CMS_hgg_%s_%s_bkgshape_norm",catname.c_str(),ext.c_str()),"nbkg",data->sumEntries(),0,3*data->sumEntries());
			//nBackground.removeRange(); // bug in roofit will break combine until dev branch brought in
			//double check the best pdf!
			int bestFitPdfIndex = getBestFitFunction(pdf,data,&catIndex,!verbose);
			catIndex.setIndex(bestFitPdfIndex);
			std::cout << "// ------------------------------------------------------------------------- //" <<std::endl; 
			std::cout << "[INFO] Created MultiPdf " << pdf->GetName() << ", in Category " << cat << " with a total of " << catIndex.numTypes() << " pdfs"<< std::endl;
			storedPdfs.Print();
			std::cout << "[INFO] Best Fit Pdf = " << bestFitPdfIndex << ", " << storedPdfs.at(bestFitPdfIndex)->GetName() << std::endl;
			std::cout << "// ------------------------------------------------------------------------- //" <<std::endl;
			std::cout << "[INFO] Simple check of index "<< simplebestFitPdfIndex <<std::endl;

			mass->setBins(nBinsForMass);
			RooDataHist dataBinned(Form("roohist_data_mass_%s",catname.c_str()),"data",*mass,*dataFull);

			// Save it (also a binned version of the dataset
			outputws->import(*pdf);
			outputws->import(nBackground);
			outputws->import(catIndex);
			outputws->import(dataBinned);
			outputws->import(*data);
			plot(mass,pdf,&catIndex,data,Form("%s/multipdf_%s",outDir.c_str(),catname.c_str()),flashggCats_,cat,bestFitPdfIndex);

		}

  }
  if (saveMultiPdf){
    outputfile->cd();
    outputws->Write();
    outputfile->Close();	
  }

  FILE *dfile = fopen(datfile.c_str(),"w");
  if (!dfile) {
    std::cerr << "[ERROR] Failed to open dat file for writing: " << datfile << std::endl;
    return 1;
  }
  
  cout << "[RESULT] Recommended options" << endl;

  for (int cat=startingCategory; cat<ncats; cat++){
    cout << "Cat " << cat << endl;
    fprintf(dfile,"cat=%d\n",(cat+catOffset)); 
    for (map<string,int>::iterator it=choices_vec[cat-startingCategory].begin(); it!=choices_vec[cat-startingCategory].end(); it++){
      cout << "\t" << it->first << " - " << it->second << endl;
      fprintf(dfile,"truth=%s:%d:%s%d\n",it->first.c_str(),it->second,namingMap[it->first].c_str(),it->second);
    }
    for (map<string,std::vector<int> >::iterator it=choices_envelope_vec[cat-startingCategory].begin(); it!=choices_envelope_vec[cat-startingCategory].end(); it++){
      std::vector<int> ords = it->second;
      for (std::vector<int>::iterator ordit=ords.begin(); ordit!=ords.end(); ordit++){
        fprintf(dfile,"paul=%s:%d:%s%d\n",it->first.c_str(),*ordit,namingMap[it->first].c_str(),*ordit);
      }
    }
    fprintf(dfile,"\n");
  }
  inFile->Close();
  if (dfile) fclose(dfile);
  fclose(resFile);
  cout << "[INFO] Program completed successfully" << endl;
  return 0;
}
