#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"

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
#include "RooAddPdf.h"

#include "RooCategory.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"

#include "../interface/PdfModelBuilder.h"
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFunc.h>
#include <iomanip>
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

// Global parameters
bool runFtestCheckWithToys = false;
int mgg_low = 100;
int mgg_high = 180;
int nBinsForMass = 4 * (mgg_high - mgg_low);

RooRealVar *intLumi_ = new RooRealVar("IntLumi", "hacked int lumi", 1000.);

TRandom3 *RandomGen = new TRandom3();

// Function to get PDF based on type and order
RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, string type, int order, const char* ext="") {
  if (type=="Bernstein") 
    return pdfsModel.getBernstein(Form("%s_bern%d", ext, order), order); 
  else if (type=="Chebychev") 
    return pdfsModel.getChebychev(Form("%s_cheb%d", ext, order), order); 
  else if (type=="Exponential") 
    return pdfsModel.getExponentialSingle(Form("%s_exp%d", ext, order), order); 
  else if (type=="PowerLaw") 
    return pdfsModel.getPowerLawSingle(Form("%s_pow%d", ext, order), order); 
  else if (type=="Laurent") 
    return pdfsModel.getLaurentSeries(Form("%s_lau%d", ext, order), order);
  else if (type=="BernsteinStepxGau") 
    return pdfsModel.getBernsteinStepxGau(Form("%s_bern%d", ext, order), order); 
  else if (type=="ExponentialStepxGau") 
    return pdfsModel.getExponentialStepxGau(Form("%s_exp%d", ext, order), order); 
  else if (type=="PowerLawStepxGau") 
    return pdfsModel.getPowerLawStepxGau(Form("%s_pow%d", ext, order), order); 
  else if (type=="LaurentStepxGau") 
    return pdfsModel.getLaurentStepxGau(Form("%s_lau%d", ext, order), order);
  else if (type=="ExpModGauss") {
    if (order == 1) 
      return pdfsModel.getExpModGaussian(Form("%s_expmodgauss", ext));
    else 
      return NULL;
  }
  else {
    cerr << "[ERROR] -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
}

// Function to run fit for both RooDataSet and RooDataHist
void runFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries) {
  int ntries = 0;
  RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
  int stat = 1;
  double minnll = 10e8;
  
  while (stat != 0 && ntries < MaxTries) {
    RooFitResult *fitTest = pdf->fitTo(*data, RooFit::Save(1),
                                      RooFit::Minimizer("Minuit2", "minimize"),
                                      RooFit::SumW2Error(kTRUE));
    stat = fitTest->status();
    minnll = fitTest->minNll();
    if (stat != 0) params_test->assignValueOnly(fitTest->randomizePars());
    delete fitTest; // Properly delete fit result
    ntries++; 
  }
  *stat_t = stat;
  *NLL = minnll;
  delete params_test; // Clean up parameters
}

// Version for RooDataHist
void runFit(RooAbsPdf *pdf, RooDataHist *data, double *NLL, int *stat_t, int MaxTries) {
  int ntries = 0;
  RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
  int stat = 1;
  double minnll = 10e8;
  
  while (stat != 0 && ntries < MaxTries) {
    RooFitResult *fitTest = pdf->fitTo(*data, RooFit::Save(1),
                                      RooFit::Minimizer("Minuit2", "minimize"),
                                      RooFit::SumW2Error(kTRUE));
    stat = fitTest->status();
    minnll = fitTest->minNll();
    if (stat != 0) params_test->assignValueOnly(fitTest->randomizePars());
    delete fitTest; // Properly delete fit result
    ntries++; 
  }
  *stat_t = stat;
  *NLL = minnll;
  delete params_test; // Clean up parameters
}

// Provides goodness of fit test
double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, std::string name) {
  double prob;
  int ntoys = 100;
  // Routine to calculate the goodness of fit. 
  name += "_gofTest.pdf";
  RooRealVar norm("norm", "norm", data->sumEntries(), 0, 10E6);
  //norm.removeRange();

  RooExtendPdf *pdf = new RooExtendPdf("ext", "ext", *mpdf, norm);

  // get The Chi2 value from the data
  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2, Binning(nBinsForMass), Name("data"));

  pdf->plotOn(plot_chi2, Name("pdf"));
  int np = pdf->getParameters(*data)->getSize();

  double chi2 = plot_chi2->chiSquare("pdf", "data", np);
  std::cout << "[INFO] Calculating GOF for pdf " << pdf->GetName() << ", using " << np << " fitted parameters" << std::endl;

  // The first thing is to check if the number of entries in any bin is < 5 
  // if so, we don't rely on asymptotic approximations
  if ((double)data->sumEntries() / nBinsForMass < 5) {
    std::cout << "[INFO] Running toys for GOF test " << std::endl;
    // store pre-fit params 
    RooArgSet *params = pdf->getParameters(*data);
    RooArgSet preParams;
    params->snapshot(preParams);
    int ndata = data->sumEntries();
 
    int npass = 0;
    std::vector<double> toy_chi2;
    for (int itoy = 0; itoy < ntoys; itoy++) {
      params->assignValueOnly(preParams);
      int nToyEvents = RandomGen->Poisson(ndata);
      RooDataHist *binnedtoy = pdf->generateBinned(RooArgSet(*mass), nToyEvents, 0, 1);
      
      // Check if binnedtoy is valid
      if (!binnedtoy) {
        std::cerr << "[ERROR] Failed to generate binned toy" << std::endl;
        continue;
      }
      
      // Fit with proper error handling
      RooFitResult *fitResult = nullptr;
      try {
        fitResult = pdf->fitTo(*binnedtoy, RooFit::Minimizer("Minuit2", "minimize"), 
                            RooFit::Minos(0), RooFit::Hesse(0), RooFit::PrintLevel(-1), 
                            RooFit::Strategy(0), RooFit::SumW2Error(kTRUE), RooFit::Save(true));
      } catch (std::exception& e) {
        std::cerr << "[ERROR] Fit failed: " << e.what() << std::endl;
        delete binnedtoy;
        continue;
      }
      
      // Check if fit was successful
      if (!fitResult || fitResult->status() != 0) {
        delete fitResult;
        delete binnedtoy;
        continue;
      }
      delete fitResult;

      RooPlot *plot_t = mass->frame();
      binnedtoy->plotOn(plot_t);
      pdf->plotOn(plot_t);

      double chi2_t = plot_t->chiSquare(np);
      if (chi2_t >= chi2) npass++;
      toy_chi2.push_back(chi2_t * (nBinsForMass - np));
      delete plot_t;
      delete binnedtoy;
    }
    std::cout << "[INFO] complete" << std::endl;
    if (ntoys > 0) {
      prob = (double)npass / ntoys;
    } else {
      prob = 0;
    }

    // Sort chi2 values for median calculation
    std::sort(toy_chi2.begin(), toy_chi2.end());
    
    TCanvas *can = new TCanvas();
    double medianChi2 = toy_chi2.size() > 0 ? 
                         toy_chi2[std::min((int)(toy_chi2.size() / 2), (int)(toy_chi2.size() - 1))] : 0;
    double rms = TMath::Sqrt(medianChi2);

    TH1F toyhist(Form("gofTest_%s.pdf", pdf->GetName()), ";Chi2;", 50, 
                 std::max(0.0, medianChi2 - 5 * rms), medianChi2 + 5 * rms);
    for (std::vector<double>::iterator itx = toy_chi2.begin(); itx != toy_chi2.end(); itx++) {
      toyhist.Fill((*itx));
    }
    toyhist.Draw();

    if (chi2 * (nBinsForMass - np) <= toyhist.GetXaxis()->GetXmax()) {
      TArrow lData(chi2 * (nBinsForMass - np), toyhist.GetMaximum(), 
                  chi2 * (nBinsForMass - np), 0);
      lData.SetLineWidth(2);
      lData.Draw();
    }
    can->SaveAs(name.c_str());
    delete can;

    // back to best fit 
    params->assignValueOnly(preParams);
  } else {
    prob = TMath::Prob(chi2 * (nBinsForMass - np), nBinsForMass - np);
  }
  std::cout << "[INFO] Chi2 in Observed =  " << chi2 * (nBinsForMass - np) << std::endl;
  std::cout << "[INFO] p-value  =  " << prob << std::endl;
  
  delete plot_chi2;
  delete pdf;
  return prob;
}

/**
 * Implementation of the spurious signal test function
 * @param cat        - Category name (e.g. "VBF0")
 * @param funcType   - Type of function to use (e.g. "BernsteinStepxGau")
 * @param order      - Order of the function
 * @param sig        - Signal strength multiplier
 * @param runPeriod  - Run period ("run2" or "run3")
 * @return true if test passes, false otherwise
 */
bool SpurialSignalTest(TString cat, TString funcType, int order, int sig, TString runPeriod, TString outDir) {
  std::cout << "[INFO] Running spurious signal test for category " << cat 
            << " with function " << funcType << " order " << order
            << " (sig = " << sig << ", " << runPeriod << ")" << std::endl;
  
  // Store function type as bkg_fun for compatibility with rest of the code
  TString bkg_fun = funcType;
  
  std::cout << "[INFO] Using " << bkg_fun << " function with order " << order << std::endl;
  
  // Mass range and binning setup
  double mgg_low = 100, mgg_high = 180, bin_times = 4;
  double bin_size = (mgg_high - mgg_low) * bin_times;
  
  // Output container for various objects we need to clean up
  std::vector<TObject*> cleanup_objects;

  // Load background MC template from new location with new naming convention
  TH1F *hbkg = nullptr, *hsig = nullptr;
  TFile* fbkg = TFile::Open("/eos/home-j/jiehan/root/templates/template_all.root");
  if (!fbkg || fbkg->IsZombie()) {
    std::cerr << "[ERROR] Failed to open template file" << std::endl;
    return false;
  }
  
  // Use new template naming - format: bkg_run2_VBF0 or bkg_run3_VBF0
  TString bkgHistName = Form("bkg_%s_%s", runPeriod.Data(), cat.Data());
  if (fbkg->GetListOfKeys()->Contains(bkgHistName)) {
    TH1F* temp = (TH1F*)fbkg->Get(bkgHistName);
    if (!temp) {
      std::cerr << "[ERROR] Failed to retrieve background histogram: " << bkgHistName << std::endl;
      fbkg->Close();
      return false;
    }
    // Clone to own the object
    hbkg = (TH1F*)temp->Clone(Form("%s_clone", bkgHistName.Data()));
    cleanup_objects.push_back(hbkg);
  } else {
    std::cerr << "[ERROR] Background template not found: " << bkgHistName << std::endl;
    fbkg->Close();
    return false;
  }
  
  double dataevents = hbkg->Integral();
  double mcsbevents = hbkg->Integral(0, (122 - mgg_low) * bin_times) + 
                      hbkg->Integral(bin_times * (mgg_high - 128), bin_size);

  // Get signal from template file - NEW METHOD
  TString signalFileName = Form("/eos/user/j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Signal/outdir_combinedPDFs/CMS-HGG_combinedPDFs_%s.root", cat.Data());
  TFile* signalFile = TFile::Open(signalFileName);
  if (!signalFile || signalFile->IsZombie()) {
    std::cerr << "[ERROR] Failed to open signal template file: " << signalFileName << std::endl;
    fbkg->Close();
    for (auto obj : cleanup_objects) delete obj;
    return false;
  }
  
  // Get the workspace and the signal PDF
  RooWorkspace* ws = (RooWorkspace*)signalFile->Get("wsig_13TeV");
  if (!ws) {
    std::cerr << "[ERROR] Workspace not found in file: " << signalFileName << std::endl;
    signalFile->Close();
    fbkg->Close();
    for (auto obj : cleanup_objects) delete obj;
    return false;
  }

  RooRealVar *mass = ws->var("CMS_hgg_mass");
  if (!mass) {
    std::cerr << "[ERROR] Mass variable not found in workspace" << std::endl;
    signalFile->Close();
    fbkg->Close();
    for (auto obj : cleanup_objects) delete obj;
    return false;
  }
  
  RooRealVar *mh = ws->var("MH");
  if (!mh) {
    std::cerr << "[ERROR] MH variable not found in workspace" << std::endl;
    signalFile->Close();
    fbkg->Close();
    for (auto obj : cleanup_objects) delete obj;
    return false;
  }
  mh->setVal(125);
  
  // Load the signal PDF from workspace - should be named something like 'sigpdf_<category>_13TeV'
  RooAbsPdf* signalPdf = ws->pdf(Form("combinedSigPdf_%s", cat.Data()));
  if (!signalPdf) {
    std::cerr << "[ERROR] Signal PDF not found in workspace" << std::endl;
    signalFile->Close();
    fbkg->Close();
    for (auto obj : cleanup_objects) delete obj;
    return false;
  }
  
  // Fix all signal PDF parameters (freeze the shape)
  RooArgSet* signalParams = signalPdf->getParameters(RooArgSet(*mass));
  TIterator* sigParamIter = signalParams->createIterator();
  RooRealVar* param;
  
  std::cout << "[INFO] Fixing all signal PDF parameters:" << std::endl;
  while ((param = (RooRealVar*)sigParamIter->Next())) {
    if (!param->isConstant()) {
      std::cout << "  - Fixing parameter: " << param->GetName() << std::endl;
      param->setConstant(true);
    }
  }
  delete sigParamIter;
  
  // Get expected signal yield from workspace
  TString signalHistName = Form("sig_%s_%s", runPeriod.Data(), cat.Data());
  if (fbkg->GetListOfKeys()->Contains(signalHistName)) {
    TH1F* temp = (TH1F*)fbkg->Get(signalHistName);
    if (!temp) {
      std::cerr << "[ERROR] Failed to retrieve signal histogram: " << signalHistName << std::endl;
      delete signalParams;
      signalFile->Close();
      fbkg->Close();
      for (auto obj : cleanup_objects) delete obj;
      return false;
    }
    // Clone to own the object
    hsig = (TH1F*)temp->Clone(Form("%s_clone", signalHistName.Data()));
    cleanup_objects.push_back(hsig);
  } else {
    std::cerr << "[ERROR] Signal template not found: " << signalHistName << std::endl;
    delete signalParams;
    signalFile->Close();
    fbkg->Close();
    for (auto obj : cleanup_objects) delete obj;
    return false;
  }
  
  double sigevents = hsig->Integral();
  std::cout << "[INFO] Retrieved signal model with expected yield: " << sigevents << std::endl;

  // Load data (full range and sidebands) from the same template file
  TH1F* hfr = nullptr;
  TString frHistName = Form("data_full_%s_%s", runPeriod.Data(), cat.Data());
  if (fbkg->GetListOfKeys()->Contains(frHistName)) {
    TH1F* temp = (TH1F*)fbkg->Get(frHistName);
    if (!temp) {
      std::cerr << "[ERROR] Failed to retrieve full range data histogram: " << frHistName << std::endl;
      delete signalParams;
      signalFile->Close();
      fbkg->Close();
      for (auto obj : cleanup_objects) delete obj;
      return false;
    }
    // Clone to own the object
    hfr = (TH1F*)temp->Clone(Form("%s_clone", frHistName.Data()));
    cleanup_objects.push_back(hfr);
  } else {
    std::cerr << "[ERROR] Full range data template not found: " << frHistName << std::endl;
    delete signalParams;
    signalFile->Close();
    fbkg->Close();
    for (auto obj : cleanup_objects) delete obj;
    return false;
  }
  double frevents = hfr->Integral();

  TH1F* hsb = nullptr;
  TString sbHistName = Form("data_%s_%s", runPeriod.Data(), cat.Data());
  if (fbkg->GetListOfKeys()->Contains(sbHistName)) {
    TH1F* temp = (TH1F*)fbkg->Get(sbHistName);
    if (!temp) {
      std::cerr << "[ERROR] Failed to retrieve sideband data histogram: " << sbHistName << std::endl;
      delete signalParams;
      signalFile->Close();
      fbkg->Close();
      for (auto obj : cleanup_objects) delete obj;
      return false;
    }
    // Clone to own the object
    hsb = (TH1F*)temp->Clone(Form("%s_clone", sbHistName.Data()));
    cleanup_objects.push_back(hsb);
  } else {
    std::cerr << "[ERROR] Sideband data template not found: " << sbHistName << std::endl;
    delete signalParams;
    signalFile->Close();
    fbkg->Close();
    for (auto obj : cleanup_objects) delete obj;
    return false;
  }
  double sbevents = hsb->Integral();

  // We can close the template file now that we've cloned all histograms
  fbkg->Close();

  // Scale background to match sideband data
  if (mcsbevents > 0) {
    hbkg->Scale(sbevents / mcsbevents);
  } else {
    std::cerr << "[ERROR] mcsbevents is zero, cannot scale background" << std::endl;
    delete signalParams;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    return false;
  }

  cout << "\n\t[INFO] Finished preparing data!!!\n" << endl;

  mass->setMin(mgg_low);
  mass->setMax(mgg_high);

  // Create RooDataHist objects
  RooDataHist* dbkg = new RooDataHist("bkg_mc", "dataset with x", *mass, hbkg);
  RooDataHist* dsb = new RooDataHist("data_sb", "dataset with x", *mass, hsb);
  RooDataHist* dfr = new RooDataHist("data_fr", "dataset with x", *mass, hfr);
  
  cleanup_objects.push_back(dbkg);
  cleanup_objects.push_back(dsb);
  cleanup_objects.push_back(dfr);

  // Set up signal and background normalization variables
  RooRealVar nsig("nsig", "nsig", sigevents, -100 * sigevents, 100 * sigevents);
  RooRealVar nbkg("nbkg", "nbkg", dataevents, 0.01 * dataevents, 2 * dataevents);

  // Create output files with updated naming
  cout << "\n\t[INFO] Creating output files at " << outDir << endl;
  TString outFilePath = Form("%s/%s_%s_%s%d_%dxsig", outDir.Data(),
                             runPeriod.Data(), cat.Data(), funcType.Data(), order, sig);
  ofstream output(outFilePath + ".txt", ofstream::app);
  
  // Initialize PDF model builder
  PdfModelBuilder pdfsModel;
  pdfsModel.setObsVar(mass);

  // Variables to track fit results
  int flag = 1;
  TString status = "Pass";
  double chi2 = 0, prob = 0, nll = 0;
  double dmc = 0, dss = 0, ss = 0, ss_mc = 0;
  double tot_err = 0, ss_cor = 0, delta = 0;
  int fit_status = 0, tries = 0;
  
  // Set up ranges for mass variable
  mass->setRange("range_low", mgg_low, 122);
  mass->setRange("signal", 122, 128);
  mass->setRange("range_high", 128, mgg_high);
  
  // Create background model
  RooAbsPdf* bkg_model = getPdf(pdfsModel, bkg_fun.Data(), order, 
                              Form("sstest_pdf_%s_%s%d", cat.Data(), funcType.Data(), order));
  if (bkg_model == NULL) {
    cout << "[ERROR] Failed to create background PDF" << endl;
    delete signalParams;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    return false;
  }
  
  // Fit background model to background template
  int bkg_npars = 0;
  int bkg_ndof = 0;
  RooFitResult* bkg_model_fit = nullptr;
  
  // Full range fitting with error handling
  try {
    bkg_model_fit = bkg_model->fitTo(*dbkg, Save(1), Minimizer("Minuit2", "minimize"), 
                                    SumW2Error(kTRUE), EvalErrorWall(false));
  } catch (std::exception& e) {
    std::cerr << "[ERROR] Background fit failed: " << e.what() << std::endl;
    delete signalParams;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    delete bkg_model;
    return false;
  }
  
  if (!bkg_model_fit || bkg_model_fit->status() != 0) {
    int maxRetries = 3;
    bool fitSucceeded = false;
    
    for (int retry = 0; retry < maxRetries; retry++) {
      std::cout << "[INFO] Retrying background fit (attempt " << retry+1 << "/" << maxRetries << ")" << std::endl;
      
      // Create new initial values for parameters
      RooArgSet* params = bkg_model->getParameters(*dbkg);
      TIterator* iter = params->createIterator();
      RooRealVar* var;
      
      while ((var = (RooRealVar*)iter->Next())) {
        if (!var->isConstant()) {
          double val = var->getVal();
          double err = var->getError();
          // Randomize starting value within ±2σ
          var->setVal(val + gRandom->Gaus(0, 2*err));
        }
      }
      delete iter;
      delete params;
      
      // Try fit again
      try {
        delete bkg_model_fit;
        bkg_model_fit = bkg_model->fitTo(*dbkg, Save(1), Minimizer("Minuit2", "minimize"), 
                                        SumW2Error(kTRUE), EvalErrorWall(false));
        
        if (bkg_model_fit && bkg_model_fit->status() == 0) {
          fitSucceeded = true;
          break;
        }
      } catch (std::exception& e) {
        std::cerr << "[ERROR] Background fit retry failed: " << e.what() << std::endl;
      }
    }
    
    if (!fitSucceeded) {
      std::cerr << "[ERROR] Background fit failed after " << maxRetries << " attempts" << std::endl;
      delete signalParams;
      signalFile->Close();
      for (auto obj : cleanup_objects) delete obj;
      delete bkg_model;
      delete bkg_model_fit;
      return false;
    }
  }
  
  bkg_npars = bkg_model_fit->floatParsFinal().getSize();
  bkg_ndof = bin_size - bkg_npars;
  
  // Plot background fit
  RooPlot* frame_bkg = mass->frame(Title(Form("Background with %s pdf", bkg_fun.Data())));
  if (!frame_bkg) {
    std::cerr << "[ERROR] Failed to create frame for background fit" << std::endl;
    delete signalParams;
    delete bkg_model_fit;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    delete bkg_model;
    return false;
  }
  
  dbkg->plotOn(frame_bkg);
  bkg_model->plotOn(frame_bkg);
  bkg_model->paramOn(frame_bkg, Layout(0.34, 0.96, 0.89), Format("NEA", AutoPrecision(1)));
  frame_bkg->getAttText()->SetTextSize(0.03);
  bkg_model->SetName(bkg_fun);
  
  nll = bkg_model_fit->minNll();
  output << "\t" << bkg_fun.Data() << "\tbkg:\tnpars = " << bkg_npars 
         << " \tchi^2 = " << frame_bkg->chiSquare(bkg_npars) 
         << "\tprob = " << TMath::Prob(frame_bkg->chiSquare(bkg_npars) * bkg_ndof, bkg_ndof) 
         << "\tnll: " << nll << endl;
  
  frame_bkg->Draw();
  gPad->Print(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/test/mc_bkg_shape_%s_%s_%s%d.png",
               runPeriod.Data(), cat.Data(), funcType.Data(), order));

  cout << "\t=================================\n";
  cout << "\n\t[INFO] Finished background function fit\n" << endl;

  // Inject signal
  TH1F* hdata = (TH1F*)hbkg->Clone("hdata");
  cleanup_objects.push_back(hdata);
  
  // Reset signal histogram for reuse
  hsig->Reset();
  
  // Ensure proper normalization of the signal PDF
  RooAbsReal* normIntegral = signalPdf->createIntegral(*mass);
  if (!normIntegral) {
    std::cerr << "[ERROR] Failed to create normalization integral" << std::endl;
    delete signalParams;
    delete bkg_model_fit;
    delete frame_bkg;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    delete bkg_model;
    return false;
  }
  
  double normVal = normIntegral->getVal();
  std::cout << "[INFO] Signal PDF normalization integral: " << normVal << std::endl;
  
  // Make sure mass variable is set to the correct range
  mass->setRange(mgg_low, mgg_high);
  
  // Check for invalid normalization value
  if (normVal <= 0) {
    std::cerr << "[ERROR] Invalid normalization value: " << normVal << std::endl;
    delete normIntegral;
    delete signalParams;
    delete bkg_model_fit;
    delete frame_bkg;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    delete bkg_model;
    return false;
  }
  
  // Fill signal histogram from PDF
  for (int i = 1; i <= hsig->GetNbinsX(); i++) {
    double x = hsig->GetBinCenter(i);
    mass->setVal(x);
    double val = signalPdf->getVal() * sigevents * hsig->GetBinWidth(i) / normVal;
    hsig->SetBinContent(i, val);
  }
  
  // Add signal to background
  hdata->Add(hsig, sig);
  hdata->SetName("asimov_data");
  
  RooDataHist* ddata = new RooDataHist("data_bin", "dataset with x", *mass, hdata);
  if (!ddata) {
    std::cerr << "[ERROR] Failed to create data histogram" << std::endl;
    delete normIntegral;
    delete signalParams;
    delete bkg_model_fit;
    delete frame_bkg;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    delete bkg_model;
    return false;
  }
  cleanup_objects.push_back(ddata);

  // Create combined model
  RooAddPdf* model = new RooAddPdf("model", "model", RooArgList(*signalPdf, *bkg_model), 
                                 RooArgList(nsig, nbkg));
  if (!model) {
    std::cerr << "[ERROR] Failed to create combined model" << std::endl;
    delete normIntegral;
    delete signalParams;
    delete bkg_model_fit;
    delete frame_bkg;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    delete bkg_model;
    return false;
  }

  // Setup canvas for plotting
  TCanvas* canv = new TCanvas();
  RooPlot* frame_data = mass->frame();
  RooPlot* frame_data_trash = mass->frame();
  
  if (!frame_data || !frame_data_trash) {
    std::cerr << "[ERROR] Failed to create data frames" << std::endl;
    delete normIntegral;
    delete signalParams;
    delete bkg_model_fit;
    delete frame_bkg;
    delete model;
    delete canv;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    delete bkg_model;
    return false;
  }
  
  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  TPad* pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.35);
  pad1->SetBottomMargin(0.18);
  pad2->SetBottomMargin(0.25);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  // Fit combined model to data with error handling
  RooFitResult* model_fit = nullptr;
  try {
    model_fit = model->fitTo(*ddata, Save(1), Minimizer("Minuit2", "minimize"), 
                            SumW2Error(kTRUE), PrintLevel(-1));
  } catch (std::exception& e) {
    std::cerr << "[ERROR] Combined model fit failed: " << e.what() << std::endl;
    delete normIntegral;
    delete signalParams;
    delete bkg_model_fit;
    delete frame_bkg;
    delete frame_data;
    delete frame_data_trash;
    delete model;
    delete canv;
    delete pad1;
    delete pad2;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    delete bkg_model;
    return false;
  }
  
  if (!model_fit) {
    std::cerr << "[ERROR] Combined model fit failed to return a valid result" << std::endl;
    delete normIntegral;
    delete signalParams;
    delete bkg_model_fit;
    delete frame_bkg;
    delete frame_data;
    delete frame_data_trash;
    delete model;
    delete canv;
    delete pad1;
    delete pad2;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    delete bkg_model;
    return false;
  }
  
  fit_status = model_fit->status();

  // Fix background parameters and only fit signal with proper handling
  RooArgSet* floatPars = nullptr;
  try {
    floatPars = dynamic_cast<RooArgSet*>(model->getParameters(*ddata)->selectByAttrib("Constant", false));
  } catch (std::exception& e) {
    std::cerr << "[ERROR] Failed to get floating parameters: " << e.what() << std::endl;
    delete normIntegral;
    delete signalParams;
    delete bkg_model_fit;
    delete model_fit;
    delete frame_bkg;
    delete frame_data;
    delete frame_data_trash;
    delete model;
    delete canv;
    delete pad1;
    delete pad2;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    delete bkg_model;
    return false;
  }
  
  if (!floatPars) {
    std::cerr << "[ERROR] Failed to get floating parameters" << std::endl;
    delete normIntegral;
    delete signalParams;
    delete bkg_model_fit;
    delete model_fit;
    delete frame_bkg;
    delete frame_data;
    delete frame_data_trash;
    delete model;
    delete canv;
    delete pad1;
    delete pad2;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    delete bkg_model;
    return false;
  }
  
  TIterator* iter = floatPars->createIterator();
  for (RooRealVar* var = (RooRealVar*)iter->Next(); var != nullptr; 
       var = (RooRealVar*)iter->Next()) {
    if (TString(var->GetName()) != "nsig") {
      var->setConstant(true);
    }
  }
  delete iter;

  // Refit with only nsig floating
  fit_status = -1;
  tries = 0;
  RooFitResult* signal_fit = nullptr;
  
  while ((fit_status != 0) && (tries < 10)) {
    try {
      if (signal_fit) delete signal_fit;
      signal_fit = model->fitTo(*ddata, RooFit::Save(1), 
                              RooFit::Minimizer("Minuit2", "minimize"), 
                              RooFit::SumW2Error(kTRUE), RooFit::PrintLevel(-1));
      
      if (signal_fit) {
        fit_status = signal_fit->status();
      } else {
        fit_status = -1;
      }
    } catch (std::exception& e) {
      std::cerr << "[ERROR] Signal-only fit failed: " << e.what() << std::endl;
      fit_status = -1;
    }
    tries++;
    
    // Add a small delay to avoid CPU spikes
    if (fit_status != 0 && tries < 10) {
      sleep(1);
    }
  }
  
  if (fit_status != 0) {
    std::cerr << "[WARNING] Signal-only fit did not converge after " << tries << " attempts" << std::endl;
  }

  // Extract signal strength and error with scaled error
  ss_mc = nsig.getVal();
  dmc = nsig.getError();
  int data_npars = (signal_fit && signal_fit->floatParsFinal().getSize() > 0) ? 
                    signal_fit->floatParsFinal().getSize() : 1;
  int data_ndof = bin_size - data_npars;

  // Plot data vs fit
  ddata->plotOn(frame_data, Name("data"), DataError(RooAbsData::SumW2));
  RooHist* plotdata = (RooHist*)frame_data->getObject(frame_data->numItems() - 1);
  if (!plotdata) {
    std::cerr << "[ERROR] Failed to get plot data" << std::endl;
    delete normIntegral;
    delete signalParams;
    delete bkg_model_fit;
    delete model_fit;
    if (signal_fit) delete signal_fit;
    delete floatPars;
    delete frame_bkg;
    delete frame_data;
    delete frame_data_trash;
    delete model;
    delete canv;
    delete pad1;
    delete pad2;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    delete bkg_model;
    return false;
  }
  
  ddata->plotOn(frame_data_trash, Name("data"), DataError(RooAbsData::SumW2));
  model->plotOn(frame_data_trash, Name("fit"));
  chi2 = frame_data_trash->chiSquare(data_npars);
  prob = TMath::Prob(chi2 * data_ndof, data_ndof);
  
  output << "\t" << bkg_fun.Data() << "\tdata(MC):\tnpars = " << data_npars 
         << "\tchi^2 = " << chi2 << "\tprob = " << prob 
         << "\tfitting status = " << fit_status << endl;

  // Unfix background parameters for final fit
  iter = floatPars->createIterator();
  for (RooRealVar* var = (RooRealVar*)iter->Next(); var != nullptr; 
       var = (RooRealVar*)iter->Next()) {
    var->setConstant(false);
  }
  delete iter;

  // Refit with unweighted error for final result
  RooFitResult* final_fit = nullptr;
  try {
    final_fit = model->fitTo(*ddata, Save(1), Minimizer("Minuit2", "minimize"), 
                           SumW2Error(kFALSE), PrintLevel(-1));
  } catch (std::exception& e) {
    std::cerr << "[ERROR] Final fit failed: " << e.what() << std::endl;
    delete normIntegral;
    delete signalParams;
    delete bkg_model_fit;
    delete model_fit;
    if (signal_fit) delete signal_fit;
    delete floatPars;
    delete frame_bkg;
    delete frame_data;
    delete frame_data_trash;
    delete model;
    delete canv;
    delete pad1;
    delete pad2;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    delete bkg_model;
    return false;
  }

  // Fix background parameters again for signal measurement
  iter = floatPars->createIterator();
  for (RooRealVar* var = (RooRealVar*)iter->Next(); var != nullptr; 
       var = (RooRealVar*)iter->Next()) {
    if (TString(var->GetName()) != "nsig") {
      var->setConstant(true);
    }
  }
  delete iter;

  // Refit with only nsig floating (final measurement)
  fit_status = -1;
  tries = 0;
  RooFitResult* final_signal_fit = nullptr;
  
  while ((fit_status != 0) && (tries < 10)) {
    try {
      if (final_signal_fit) delete final_signal_fit;
      final_signal_fit = model->fitTo(*ddata, RooFit::Save(1), 
                                    RooFit::Minimizer("Minuit2", "minimize"), 
                                    RooFit::SumW2Error(kFALSE), RooFit::PrintLevel(-1));
      
      if (final_signal_fit) {
        fit_status = final_signal_fit->status();
      } else {
        fit_status = -1;
      }
    } catch (std::exception& e) {
      std::cerr << "[ERROR] Final signal-only fit failed: " << e.what() << std::endl;
      fit_status = -1;
    }
    tries++;
    
    // Add a small delay to avoid CPU spikes
    if (fit_status != 0 && tries < 10) {
      sleep(1);
    }
  }

  // Get signal strength and error
  ss = nsig.getVal();
  dss = nsig.getError();
  tot_err = sqrt(dss * dss + ss * ss);
  delta = fabs(ss) - 2 * dmc;
  ss_cor = (delta < 0) ? 0 : ((ss > 0) ? delta : -delta);
  
  // Test failure criterion
  if (delta > 0.2 * dss) status = "Fail";

  // Get final chi2 and probability
  data_npars = (final_signal_fit && final_signal_fit->floatParsFinal().getSize() > 0) ? 
               final_signal_fit->floatParsFinal().getSize() : 1;
  data_ndof = bin_size - data_npars;
  model->plotOn(frame_data, Name("fit"));
  chi2 = frame_data->chiSquare(data_npars);
  prob = TMath::Prob(chi2 * data_ndof, data_ndof);
  
  // Plot background component
  model->plotOn(frame_data, Name("background"), Components(bkg_model->GetName()), 
                LineStyle(ELineStyle::kDashed), LineColor(kGreen));
  RooCurve* nomBkgCurve = (RooCurve*)frame_data->getObject(frame_data->numItems() - 1);
  if (!nomBkgCurve) {
    std::cerr << "[ERROR] Failed to get background curve" << std::endl;
    delete normIntegral;
    delete signalParams;
    delete bkg_model_fit;
    delete model_fit;
    if (signal_fit) delete signal_fit;
    if (final_fit) delete final_fit;
    if (final_signal_fit) delete final_signal_fit;
    delete floatPars;
    delete frame_bkg;
    delete frame_data;
    delete frame_data_trash;
    delete model;
    delete canv;
    delete pad1;
    delete pad2;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    delete bkg_model;
    return false;
  }

  // Create legend and finalize plot
  model->SetName(Form("%s_model", bkg_fun.Data()));
  frame_data->SetTitle(Form("Pesudo data with with x%d signal, prob: %.3f", sig, prob));
  frame_data->SetXTitle("");
  frame_data->SetLabelSize(0.042, "XY");
  frame_data->SetTitleSize(0.056, "Y");
  frame_data->SetTitleOffset(0.75, "Y");
  
  output << "\t" << bkg_fun.Data() << "\tdata(Psu):\tnpars = " << data_npars 
         << "\tchi^2 = " << chi2 << "\tprob = " << prob 
         << "\tstatus = " << fit_status << endl;
  output << "\t" << bkg_fun.Data() << "\tSS:\tnsig = " << ss_mc << ":" << ss 
         << "\tdmc = " << dmc << "\tss_cor = " << ss_cor << "\tdss = " << dss 
         << "\ttot_err = " << tot_err << "\tstatus = " << status.Data() << "\n" << endl;

  // Add legend to plot
  TLegend* leg = new TLegend(0.6, 0.65, 0.88, 0.88);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->AddEntry(frame_data->findObject("data"), "MC", "ep");
  leg->AddEntry(frame_data->findObject("fit"), "Bkg + Sig", "l");
  leg->AddEntry(frame_data->findObject("background"), "Bkg", "l");
  frame_data->Draw();
  leg->Draw("same");

  // Add signal curve to plot
  signalPdf->plotOn(frame_data, RooFit::Name("signal"), 
                  RooFit::Normalization(ss, RooAbsReal::NumEvent), 
                  LineColor(kRed), LineWidth(4));
  RooCurve* nomSigCurve = (RooCurve*)frame_data->getObject(frame_data->numItems() - 1);
  if (!nomSigCurve) {
    std::cerr << "[ERROR] Failed to get signal curve" << std::endl;
    delete normIntegral;
    delete signalParams;
    delete bkg_model_fit;
    delete model_fit;
    if (signal_fit) delete signal_fit;
    if (final_fit) delete final_fit;
    if (final_signal_fit) delete final_signal_fit;
    delete floatPars;
    delete frame_bkg;
    delete frame_data;
    delete frame_data_trash;
    delete model;
    delete canv;
    delete pad1;
    delete pad2;
    delete leg;
    signalFile->Close();
    for (auto obj : cleanup_objects) delete obj;
    delete bkg_model;
    return false;
  }

  // Create bottom pad with residuals
  pad2->cd();
  int npoints = plotdata->GetN();
  double xtmp, ytmp;
  int point = 0;
  TGraphAsymmErrors* hdatasub = new TGraphAsymmErrors(npoints);
  
  for (int ipoint = 0; ipoint < npoints; ++ipoint) {
    plotdata->GetPoint(ipoint, xtmp, ytmp);
    double bkgval = nomBkgCurve->interpolate(xtmp);
    double errhi = plotdata->GetErrorYhigh(ipoint);
    double errlow = plotdata->GetErrorYlow(ipoint);

    std::cout << "[INFO] Category " << cat.Data() << " setting point " << point 
              << " : xtmp " << xtmp << "  ytmp " << ytmp 
              << " bkgval  " << bkgval << " ytmp-bkgval " << ytmp - bkgval << std::endl;
    
    hdatasub->SetPoint(point, xtmp, ytmp - bkgval);
    hdatasub->SetPointError(point, 0., 0., errlow, errhi);
    point++;
  }

  // Setup dummy histogram for axis formatting
  TH1D* hdummy = new TH1D("hdummyweight", "", mgg_high - mgg_low, mgg_low, mgg_high);
  hdummy->SetStats(0);
  
  // Check if we can safely access the histogram from hdatasub
  if (hdatasub->GetHistogram()) {
    hdummy->SetMaximum(hdatasub->GetHistogram()->GetMaximum() + 1);
    hdummy->SetMinimum(hdatasub->GetHistogram()->GetMinimum() - 1);
  } else {
    // Set default range if histogram not available
    hdummy->SetMaximum(5);
    hdummy->SetMinimum(-5);
  }
  
  hdummy->GetYaxis()->SetTitle("data - bkg PDF");
  hdummy->GetYaxis()->SetTitleOffset(0.35);
  hdummy->GetYaxis()->SetTitleSize(0.12);
  hdummy->GetYaxis()->SetLabelSize(0.09);
  hdummy->GetXaxis()->SetTitle("m_{ll#gamma} (GeV)");
  hdummy->GetXaxis()->SetTitleSize(0.12);
  hdummy->GetXaxis()->SetLabelSize(0.09);
  hdummy->Draw("HIST");
  hdummy->GetYaxis()->SetNdivisions(808);

  // Plot residuals and signal curve
  hdatasub->SetMarkerStyle(8);
  hdatasub->Draw("PESAME");
  nomSigCurve->Draw("L SAME");
  
  // Save plot safely
  bool saveFailed = false;
  try {
    canv->SaveAs(Form("%s/pesudo_data_shape_%s_%s_%s%d.png", outDir.Data(),
                      runPeriod.Data(), cat.Data(), funcType.Data(), order));
  } catch (std::exception& e) {
    std::cerr << "[ERROR] Failed to save canvas: " << e.what() << std::endl;
    saveFailed = true;
  }
  
  if (saveFailed) {
    // Try to save to a default location as fallback
    try {
      canv->SaveAs(Form("/tmp/pesudo_data_shape_%s_%s_%s%d.png",
                        runPeriod.Data(), cat.Data(), funcType.Data(), order));
      std::cout << "[INFO] Fallback: Canvas saved to /tmp directory" << std::endl;
    } catch (...) {
      std::cerr << "[ERROR] Also failed to save canvas to fallback location" << std::endl;
    }
  }

  // Clean up all resources
  delete normIntegral;
  delete signalParams;
  delete bkg_model_fit;
  delete model_fit;
  if (signal_fit) delete signal_fit;
  if (final_fit) delete final_fit;
  if (final_signal_fit) delete final_signal_fit;
  delete floatPars;
  delete frame_bkg;
  delete frame_data;
  delete frame_data_trash;
  delete model;
  delete hdummy;
  delete hdatasub;
  delete leg;
  delete canv;
  delete pad1;
  delete pad2;
  signalFile->Close();
  
  for (auto obj : cleanup_objects) {
    delete obj;
  }
  
  delete bkg_model;
  
  output.close();
  
  // Return test result
  return (status == "Pass");
}
