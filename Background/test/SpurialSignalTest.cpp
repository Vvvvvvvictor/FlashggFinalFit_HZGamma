#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

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

bool BLIND = true;
bool runFtestCheckWithToys=false;
int mgg_low = 100;
int mgg_high = 180;
int nBinsForMass = 4*(mgg_high-mgg_low);

RooRealVar *intLumi_ = new RooRealVar("IntLumi","hacked int lumi", 1000.);

TRandom3 *RandomGen = new TRandom3();

RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, string type, int order, const char* ext=""){
  
  if (type=="Bernstein") return pdfsModel.getBernstein(Form("%s_bern%d",ext,order),order); 
  else if (type=="Chebychev") return pdfsModel.getChebychev(Form("%s_cheb%d",ext,order),order); 
  else if (type=="Exponential") return pdfsModel.getExponentialSingle(Form("%s_exp%d",ext,order),order); 
  else if (type=="PowerLaw") return pdfsModel.getPowerLawSingle(Form("%s_pow%d",ext,order),order); 
  else if (type=="Laurent") return pdfsModel.getLaurentSeries(Form("%s_lau%d",ext,order),order);
  else if (type=="BernsteinStepxGau") return pdfsModel.getBernsteinStepxGau(Form("%s_bern%d",ext,order),order); 
  else if (type=="ExponentialStepxGau") return pdfsModel.getExponentialStepxGau(Form("%s_exp%d",ext,order),order); 
  else if (type=="PowerLawStepxGau") return pdfsModel.getPowerLawStepxGau(Form("%s_pow%d",ext,order),order); 
  else if (type=="LaurentStepxGau") return pdfsModel.getLaurentStepxGau(Form("%s_lau%d",ext,order),order);
  else {
    cerr << "[ERROR] -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
}

void runFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries){

	int ntries=0;
  	RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
	//params_test->Print("v");
	int stat=1;
	double minnll=10e8;
	while (stat!=0){
	  if (ntries>=MaxTries) break;
	  RooFitResult *fitTest = pdf->fitTo(*data,RooFit::Save(1)
    ,RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
          stat = fitTest->status();
	  minnll = fitTest->minNll();
	  if (stat!=0) params_test->assignValueOnly(fitTest->randomizePars());
	  ntries++; 
	}
	*stat_t = stat;
	*NLL = minnll;
}

double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, std::string name){

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
 
  if ((double)data->sumEntries()/nBinsForMass < 5 ){

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

void SpurialSignalTest(int cat = 0, int sig = 0, TString channel = "two_jet"){
  double mgg_low=100, mgg_high=180, bin_times=4;
  double bin_size = (mgg_high-mgg_low)*bin_times;
  // if (channel.EqualTo("two_jet") & (cat == 0))
  // {
  //     mgg_low = 100;
  //     mgg_high = 180;
  //     bin_size = 320;
  // }
  //background MC template
  TH1F* hbkg;
  TFile* fbkg = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_sig_template.root"));
  if (fbkg->GetListOfKeys()->Contains((Form("bkg_%s_cat%d",channel.Data(), cat))))
      hbkg = (TH1F*)fbkg->Get(Form("bkg_%s_cat%d",channel.Data(), cat));
  else
      abort();
  // TFile* fbkg = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_template_v3_cut2/bkg/bkg_0sig_cat%d.root", cat));
  // hbkg = (TH1F*)fbkg->Get(Form("mass_cat%d", cat));
  double dataevents = hbkg->Integral();
  double mcsbevents = hbkg->Integral(0,(122-mgg_low)*bin_times)+hbkg->Integral(bin_times*(mgg_high-128),bin_size);

  //signal MC
  TFile* fsig = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_sig_template.root"));
  TH1F* hsig = (TH1F*)( fsig->Get(Form("sig_%s_cat%d", channel.Data(), cat)) );
  double sigevents = hsig->Integral();

  //data full range
  TFile* ffr = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_sig_template.root"));
  TH1F* hfr = (TH1F*)( ffr->Get(Form("data_full_%s_cat%d", channel.Data(), cat)) );
  double frevents = hfr->Integral();

  //data side band
  TFile* fsb = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_sig_template.root"));
  TH1F* hsb = (TH1F*)( fsb->Get(Form("data_%s_cat%d", channel.Data(), cat)) );
  double sbevents = hsb->Integral();

  hbkg->Scale(sbevents/mcsbevents);

  cout << "\n\tFinish preparing data!!!\n" << endl;

  //initializing
  RooRealVar mH("mH", "mH", 125.38, mgg_low, mgg_high);
  mH.setMin( mgg_low );
  mH.setMax( mgg_high );
  // mH.setMax( hbkg->GetBinCenter(hbkg->GetNbinsX()) + hbkg->GetBinWidth(hbkg->GetNbinsX()) );

  //background function fit
  RooDataHist* dbkg = new RooDataHist("bkg_mc","dataset with x", mH, hbkg);
  RooDataHist* dsb = new RooDataHist("data_sb","dataset with x", mH, hsb);
  RooDataHist* dfr = new RooDataHist("data_fr","dataset with x", mH, hfr);
  // cout<<"bkg hist integral "<<hbkg->Integral()<<" "<<dbkg->sumEntries()<<endl;

  RooRealVar nsig("nsig","nsig",sigevents,-100*sigevents,100*sigevents);
  RooRealVar nbkg("nbkg","nbkg",dataevents, 0.01*dataevents, 2*dataevents);

  ofstream output(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/outputs/%s_%d_%dxsig.txt", channel.Data(), cat, sig), ofstream::app);

  TFile *f = new TFile(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/outputs/%s_%d_%dxsig.root", channel.Data(), cat, sig),"UPDATE"); 
  hbkg->Write(hbkg->GetName(), TObject::kOverwrite);
  hsb->Write(hsb->GetName(), TObject::kOverwrite);
  hfr->Write(hfr->GetName(), TObject::kOverwrite);
  hsig->Write(hsig->GetName(), TObject::kOverwrite);

  PdfModelBuilder pdfsModel;
  std:: cout << "[INFO] Got mass from ws " << mH.getVal() << endl;
  pdfsModel.setObsVar(&mH);

  int flag; TString status;
  double dmc, dss, ss, ss_mc, tot_err, ss_cor, delta, chi2, prob, nll;
  RooAbsCollection *floatPars; TIterator *iter;

  int bkg_npars, bkg_ndof;
  RooPlot *frame_bkg;

  RooAbsPdf* bkg_model;
  RooFitResult *model_fit, *bkg_model_fit;
  int data_npars, data_ndof, fit_status, tries;

  //signal function fit
  RooAbsPdf* signal = pdfsModel.getDoubleCB();
  RooDataHist* dsig = new RooDataHist("data_bin","dataset with x", mH, hsig);
  RooFitResult *signal_fit;
  signal_fit = signal->fitTo(*dsig,Save(1),Minimizer("Minuit2","minimize"),SumW2Error(kTRUE),EvalErrorWall(false)); //FIXME kTRUE or kFALSE
  int sig_npars = signal_fit->floatParsFinal().getSize();
  int sig_ndof = bin_size-sig_npars;
  RooPlot *frame_sig = mH.frame(Title("Signal with distorted Gaussian pdf"));
  dsig->plotOn(frame_sig, DataError(RooAbsData::SumW2));
  signal->plotOn(frame_sig);
  signal->paramOn(frame_sig, Layout(0.4,0.96,0.89),Format("NEA",AutoPrecision(1)));
  signal->Write(signal->GetName(), TObject::kOverwrite);
  output << "\tsig:\tnpars = " << sig_npars << "\tchi^2 = " << frame_sig->chiSquare(sig_npars) << "\tprob = " << TMath::Prob(frame_sig->chiSquare(sig_npars)*sig_ndof, sig_ndof) << endl;
  frame_sig->Draw();
  gPad->Print(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/test/signal_shape_%s_cat%d.png",channel.Data(),cat));

  // fix the float parameters in model
  floatPars = signal->getParameters(*dsig)->selectByAttrib("Constant", false);
  iter = floatPars->createIterator();
  for (RooRealVar* var = (RooRealVar*)iter->Next(); var != nullptr; var = (RooRealVar*)iter->Next()) {
      var->setConstant(true);
  }

  cout << "\t=================================" << endl;
  cout << "\n\t Finish signal function fit\n" << endl;

  vector<string> functionClasses;
  functionClasses.push_back("LaurentStepxGau");
  functionClasses.push_back("BernsteinStepxGau");
  functionClasses.push_back("ExponentialStepxGau");
  functionClasses.push_back("PowerLawStepxGau");
  map<string,string> namingMap;
  namingMap.insert(pair<string,string>("BernsteinStepxGau","bern"));
  namingMap.insert(pair<string,string>("ExponentialStepxGau","exp"));
  namingMap.insert(pair<string,string>("PowerLawStepxGau","pow"));
  namingMap.insert(pair<string,string>("LaurentStepxGau","lau"));

  mH.setRange("range_low",mgg_low,122);
  mH.setRange("signal",122,128);
  mH.setRange("range_high",128,mgg_high);

  int order, catOffset=0; TString bkg_fun; string ext;

  for (vector<string>::iterator func=functionClasses.begin(); func!=functionClasses.end(); func++){
    flag = 1; status = "Pass"; chi2=0;
    order = 1;
    bkg_fun = *func; ext = namingMap[*func];
    while (order<7) {
      bkg_model = getPdf(pdfsModel,*func,order,Form("sstest_pdf_%d_%s",(cat+catOffset),ext.c_str()));
      if (bkg_model == NULL) {
        cout << "[ERROR] -- getPdf() -- pdf is NULL. Exiting." << endl;
        order++;
        continue;
      }

      //background function fit
      bkg_model_fit = bkg_model->fitTo(*dfr,Range("range"),SplitRange(true),Save(1),Minimizer("Minuit2","minimize"),SumW2Error(kTRUE),EvalErrorWall(false)); 
      bkg_npars = bkg_model_fit->floatParsFinal().getSize();
      frame_bkg = mH.frame(Title(Form("Data side band with %s pdf", bkg_fun.Data())));
      bkg_ndof = bin_size-bin_times*7-bkg_npars;
      dfr->plotOn(frame_bkg, Cut("mH>128 || mH<122"));
      bkg_model->plotOn(frame_bkg,NormRange("range_low,signal,range_high"));
      chi2 = frame_bkg->chiSquare(bkg_npars);
      bkg_model->SetName(bkg_fun);
      bkg_model->Write(bkg_model->GetName(), TObject::kOverwrite);
      nll = bkg_model_fit->minNll();
      prob = TMath::Prob(chi2*bkg_ndof, bkg_ndof);
      // if(prob<0.01) status = "Fail";
      output << "\t" << bkg_fun.Data() << "\tsb:\tnpars = " << bkg_npars << " \tchi^2 = " << chi2 << "\tprob = " << prob << "\tnll: " << nll << endl;
      bkg_model->paramOn(frame_bkg, Layout(0.34,0.96,0.89),Format("NEA",AutoPrecision(1)));
      frame_bkg->getAttText()->SetTextSize(0.03);
      frame_bkg->Draw();
      gPad->Print(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/test/data_sb_shape_%s_cat%d_%s.png",channel.Data(),cat,bkg_fun.Data()));

      cout << "\t=================================" << endl;
      cout << "\n\t Finish data side band fit\n" << endl;

      //background function fit
      bkg_model_fit = bkg_model->fitTo(*dbkg,Save(1),Minimizer("Minuit2","minimize"),SumW2Error(kTRUE),EvalErrorWall(false));
      bkg_npars = bkg_model_fit->floatParsFinal().getSize();
      bkg_ndof = bin_size-bkg_npars;
      frame_bkg = mH.frame(Title(Form("Background with %s pdf", bkg_fun.Data())));
      dbkg->plotOn(frame_bkg);
      bkg_model->plotOn(frame_bkg);
      bkg_model->paramOn(frame_bkg, Layout(0.34,0.96,0.89),Format("NEA",AutoPrecision(1)));
      frame_bkg->getAttText()->SetTextSize(0.03);
      bkg_model->SetName(bkg_fun);
      bkg_model->Write(bkg_model->GetName(), TObject::kOverwrite);
      nll = bkg_model_fit->minNll();
      output << "\t" << bkg_fun.Data() << "\tbkg:\tnpars = " << bkg_npars << " \tchi^2 = " << frame_bkg->chiSquare(bkg_npars) << "\tprob = " << TMath::Prob(frame_bkg->chiSquare(bkg_npars)*bkg_ndof, bkg_ndof) << "\tnll: " << nll << endl;
      frame_bkg->Draw();
      // RooHist *hpull = frame_bkg->pullHist();
      // RooPlot *frame3 = mH.frame(Title("Pull Distribution"));
      // frame3->addPlotable(hpull, "P");
      // frame3->Draw();
      gPad->Print(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/test/mc_bkg_shape_%s_cat%d_%s.png",channel.Data(),cat,bkg_fun.Data()));

      // RooAbsCollection *m_bkgParameters = bkg_model->getParameters(RooArgSet())->selectByAttrib("Constant", false);
      // TIterator *it = m_bkgParameters->createIterator();
      // for (RooRealVar *p = (RooRealVar *)it->Next(); p != 0; p = (RooRealVar *)it->Next()) {
      //     p->setConstant(kTRUE);
      // }

      cout << "\t=================================\n";
      cout << "\n\tFinish background function fit\n" << endl;

      //inject signal
      TH1F* hdata = (TH1F*)hbkg->Clone();
      hdata->Add(hsig, sig);
      hdata->SetName("asimov data");
      hdata->Write(hdata->GetName(), TObject::kOverwrite);
      RooDataHist* ddata = new RooDataHist("data_bin","dataset with x", mH, hdata);

      // asimov dataset fit
      RooAddPdf* model = new RooAddPdf("model","model",RooArgList(*signal, *bkg_model),RooArgList(nsig,nbkg));

      TCanvas *canv = new TCanvas();
      RooPlot *frame_data, *frame_data_trash;
      frame_data = mH.frame();
      frame_data_trash = mH.frame();
      TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
      TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.35);
      pad1->SetBottomMargin(0.18);
      // pad2->SetTopMargin(0.00001);
      pad2->SetBottomMargin(0.25);
      pad1->Draw();
      pad2->Draw();
      pad1->cd();

      model_fit = model->fitTo(*ddata,Save(1),Minimizer("Minuit2","minimize"),SumW2Error(kTRUE)); //FIXME kTRUE or kFALSE
      fit_status = model_fit->status();

      // fix the float parameters in model except nsig
      floatPars = model->getParameters(*ddata)->selectByAttrib("Constant", false);
      iter = floatPars->createIterator();
      for (RooRealVar* var = (RooRealVar*)iter->Next(); var != nullptr; var = (RooRealVar*)iter->Next()) {
          if (var->GetName() != TString("nsig")) {
          var->setConstant(true);
          }
      }

      fit_status = -1;
      tries = 0;
      while((fit_status != 0) && (tries < 10)){
          // model_fit = model->chi2FitTo(*ddata,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),DataError(RooAbsData::SumW2)); //FIXME kTRUE or kFALSE
          model_fit = model->fitTo(*ddata,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME kTRUE or kFALSE
          fit_status = model_fit->status();
          tries++;
      }

      ss_mc = nsig.getVal();
      dmc = nsig.getError();
      data_ndof = bin_size-data_npars;

      // if(ss < 0) frame_data_trash->SetMinimum(ss);
      ddata->plotOn(frame_data, Name("data"), DataError(RooAbsData::SumW2));
      RooHist *plotdata = (RooHist*)frame_data->getObject(frame_data->numItems()-1);
      ddata->plotOn(frame_data_trash, Name("data"), DataError(RooAbsData::SumW2));
      model->plotOn(frame_data_trash, Name("fit"));
      chi2 = frame_data_trash->chiSquare(data_npars);
      prob = TMath::Prob(chi2*data_ndof, data_ndof);
      output << "\t" << bkg_fun.Data() << "\tdata(MC):\tnpars = " << data_npars << "\tchi^2 = " << chi2 << ": " << chi2 << "\tprob = " << prob << "\tfitting status = " << fit_status << endl;

      // unfix the float parameters in model
      iter = floatPars->createIterator();
      for (RooRealVar* var = (RooRealVar*)iter->Next(); var != nullptr; var = (RooRealVar*)iter->Next()) {
          var->setConstant(false);
      }

      // dmc = nsig.getError();
      // data_npars = model_fit->floatParsFinal().getSize();
      // data_ndof = bin_size-data_npars;
      // ddata->plotOn(frame_data_trash, Name("data"), DataError(RooAbsData::SumW2));
      // RooHist *plotdata = (RooHist*)frame_data_trash->getObject(frame_data_trash->numItems()-1);
      // model->plotOn(frame_data_trash, Name("fitmcweight"));
      // chi2 = frame_data_trash->chiSquare(data_npars);
      // prob = TMath::Prob(chi2*data_ndof, data_ndof);
      // frame_data_trash->remove("fitmcweight");
      // output << "\t" << bkg_fun.Data() << "\tdata(MC):\tnpars = " << data_npars << "\tchi^2 = " << chi2 << "\tprob = " << prob << "\tstatus = " << fit_status << endl;

      model_fit = model->fitTo(*ddata,Save(1),Minimizer("Minuit2","minimize"),SumW2Error(kFALSE)); //FIXME kTRUE or kFALSE

      // fix the float parameters in model except nsig
      floatPars = model->getParameters(*ddata)->selectByAttrib("Constant", false);
      iter = floatPars->createIterator();
      for (RooRealVar* var = (RooRealVar*)iter->Next(); var != nullptr; var = (RooRealVar*)iter->Next()) {
          if (var->GetName() != TString("nsig")) {
          var->setConstant(true);
          }
      }

      while((fit_status != 0) && (tries < 10)){
          model_fit = model->fitTo(*ddata,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kFALSE)); //FIXME kTRUE or kFALSE
          fit_status = model_fit->status();
          tries++;
      }

      ss = nsig.getVal();
      dss = nsig.getError();
      tot_err = sqrt(dss * dss + ss * ss);
      delta = fabs(ss) - 2 * dmc;
      ss_cor = (delta < 0) ? 0 : ((ss > 0) ? delta : -delta);
      if (delta > 0.2 * dss) status = "Fail";

      data_npars = model_fit->floatParsFinal().getSize();
      data_ndof = bin_size-data_npars;

      // if(ss < 0) frame_data->SetMinimum(ss);
      // ddata->plotOn(frame_data, Name("data"), DataError(RooAbsData::SumW2));
      model->plotOn(frame_data, Name("fit"));
      chi2 = frame_data->chiSquare(data_npars);
      prob = TMath::Prob(chi2*data_ndof, data_ndof);

      // RooHist *plotdata = (RooHist*)frame_data->getObject(frame_data->numItems()-1);
      // model->plotOn(frame_data, Name("signal"), Components(signal->GetName()), LineStyle(ELineStyle::kDashed), LineColor(kGreen));
      model->plotOn(frame_data, Name("background"), Components(bkg_model->GetName()), LineStyle(ELineStyle::kDashed), LineColor(kGreen));
      RooCurve* nomBkgCurve = (RooCurve*)frame_data->getObject(frame_data->numItems()-1);

      // if(prob<0.01) status="Fail";
      model->SetName(Form("%s_model", bkg_fun.Data()));
      frame_data->SetTitle(Form("Pesudo data with with x%d signal, prob: %.3f", sig, prob));
      frame_data->SetXTitle("");
      frame_data->SetLabelSize(0.042, "XY");
      frame_data->SetTitleSize(0.056, "Y");
      frame_data->SetTitleOffset(0.75, "Y");
      // model->Write(model->GetName(), TObject::kOverwrite);
      output << "\t" << bkg_fun.Data() << "\tdata(Psu):\tnpars = " << data_npars << "\tchi^2 = " << chi2 << "\tprob = " << prob << "\tstatus = " << fit_status << endl;
      output << "\t" << bkg_fun.Data() << "\tSS:\tnsig = " << ss_mc << ":" << ss << "\tdmc = " << dmc << "\tss_cor = " << ss_cor << "\tdss = " << dss << "\ttot_err = " << tot_err << "\tstatus = " << status.Data() << "\n" << endl;
      // output << "\tnbkg = " << nbkg.getVal() << "\tnbkg_err = " << nbkg.getError() << "\n" << endl;

      TLegend *leg = new TLegend(0.6,0.65,0.88,0.88);
      leg->SetFillColor(0);
      leg->SetLineColor(0);
      leg->AddEntry(frame_data->findObject("data"),"MC","ep");
      leg->AddEntry(frame_data->findObject("fit"),"Bkg + Sig","l");
      // leg->AddEntry(frame_data->findObject("signal"),"Sig","l");
      leg->AddEntry(frame_data->findObject("background"),"Bkg","l");
      frame_data->Draw();
      leg->Draw("same");

      signal->plotOn(frame_data, RooFit::Name("signal"), RooFit::Normalization(ss,RooAbsReal::NumEvent), LineColor(kRed),LineWidth(4));
      // model->plotOn(frame_data, Name("signal"), Components(signal->GetName()), LineStyle(ELineStyle::kDashed), LineColor(kGreen));
      RooCurve* nomSigCurve = (RooCurve*)frame_data->getObject(frame_data->numItems()-1);

      pad2->cd();
      int npoints = plotdata->GetN();
      double xtmp,ytmp;//
      int point =0;
      TGraphAsymmErrors *hdatasub = new TGraphAsymmErrors(npoints);
      for (int ipoint=0; ipoint<npoints; ++ipoint) {
          plotdata->GetPoint(ipoint, xtmp,ytmp);
          double bkgval = nomBkgCurve->interpolate(xtmp);
          // if ((xtmp > 122 ) && ( xtmp < 128) ) continue;
          double errhi = plotdata->GetErrorYhigh(ipoint);
          double errlow = plotdata->GetErrorYlow(ipoint);

          std::cout << "[INFO] Channel  " << channel.Data() << " setting point " << point <<" : xtmp "<< xtmp << "  ytmp " << ytmp << " bkgval  " << bkgval << " ytmp-bkgval " << ytmp-bkgval << std::endl;
          // if(fabs(ytmp)<1e-5) continue;
          hdatasub->SetPoint(point,xtmp,ytmp-bkgval);
          hdatasub->SetPointError(point,0.,0.,errlow,errhi );
          point++;
      }

      TH1 *hdummy = new TH1D("hdummyweight","",mgg_high-mgg_low,mgg_low,mgg_high);
      hdummy->SetStats(0);
      hdummy->SetMaximum(hdatasub->GetHistogram()->GetMaximum()+1);
      hdummy->SetMinimum(hdatasub->GetHistogram()->GetMinimum()-1);
      hdummy->GetYaxis()->SetTitle("data - bkg PDF");
      hdummy->GetYaxis()->SetTitleOffset(0.35);
      hdummy->GetYaxis()->SetTitleSize(0.12);
      hdummy->GetYaxis()->SetLabelSize(0.09);
      hdummy->GetXaxis()->SetTitle("m_{ll#gamma} (GeV)");//bing
      hdummy->GetXaxis()->SetTitleSize(0.12);
      hdummy->GetXaxis()->SetLabelSize(0.09);
      hdummy->Draw("HIST");
      hdummy->GetYaxis()->SetNdivisions(808);

      // TLine *line3 = new TLine(mgg_low,0.,mgg_high,0.);
      // line3->SetLineColor(kRed);
      // //line3->SetLineStyle(kDashed);
      // line3->SetLineWidth(5.0);
      // line3->Draw();

      hdatasub->SetMarkerStyle(8);
      hdatasub->Draw("PESAME");
      nomSigCurve->Draw("L SAME");
      canv->SaveAs(Form("./test/pesudo_data_shape_%s_cat%d_%s.png",channel.Data(),cat,bkg_fun.Data()));

      // gPad->Print(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/test/pesudo_data_shape_%s_cat%d_%s.png",channel.Data(),cat,bkg_fun.Data()));

      // RooHist *hpull = frame_data->pullHist();
      // RooPlot *frame3 = mH.frame(Title("Pull Distribution"));
      // frame3->addPlotable(hpull, "P");
      // frame3->Draw();
      // gPad->Print(Form("test/%s_cat%d_%s_%dsig.png",channel.Data(),cat,bkg_fun.Data(),sig));

      cout << "\t=================================" << endl;
      cout << "\n\t Finish asimov data fit\n" << endl;

      order++;
    }
  }
}

int main(){
  for (int i = 0; i < 4; i++){
    SpurialSignalTest(i, 0, "two_jet");
  }
  return 0;
}

// double getProbabilityFtest(double chi2, int ndof,RooAbsPdf *pdfNull, RooAbsPdf *pdfTest, RooRealVar *mass, RooDataSet *data, std::string name){
 
//   double prob_asym = TMath::Prob(chi2,ndof);
//   if (!runFtestCheckWithToys) return prob_asym;

//   int ndata = data->sumEntries();
  
//   // fit the pdfs to the data and keep this fit Result (for randomizing)
//   RooFitResult *fitNullData = pdfNull->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
// 		,RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE),RooFit::PrintLevel(-1)); //FIXME
//   RooFitResult *fitTestData = pdfTest->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
// 		,RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE),RooFit::PrintLevel(-1)); //FIXME

//   // Ok we want to check the distribution in toys then 
//   // Step 1, cache the parameters of each pdf so as not to upset anything 
//   RooArgSet *params_null = pdfNull->getParameters((const RooArgSet*)(0));
//   RooArgSet preParams_null;
//   params_null->snapshot(preParams_null);
//   RooArgSet *params_test = pdfTest->getParameters((const RooArgSet*)(0));
//   RooArgSet preParams_test;
//   params_test->snapshot(preParams_test);
 
//   int ntoys =5000;
//   TCanvas *can = new TCanvas();
//   can->SetLogy();
//   TH1F toyhist(Form("toys_fTest_%s.pdf",pdfNull->GetName()),";Chi2;",60,-2,10);
//   TH1I toyhistStatN(Form("Status_%s.pdf",pdfNull->GetName()),";FitStatus;",8,-4,4);
//   TH1I toyhistStatT(Form("Status_%s.pdf",pdfTest->GetName()),";FitStatus;",8,-4,4);

//   TGraph *gChi2 = new TGraph();
//   gChi2->SetLineColor(kGreen+2);
//   double w = toyhist.GetBinWidth(1);

//   int ipoint=0;

//   for (int b=0;b<toyhist.GetNbinsX();b++){
// 	double x = toyhist.GetBinCenter(b+1);
// 	if (x>0){
// 	  gChi2->SetPoint(ipoint,x,(ROOT::Math::chisquared_pdf(x,ndof)));
// 	  ipoint++;
// 	}
//   }
//   int npass =0; int nsuccesst =0;
//   mass->setBins(nBinsForMass);
//   for (int itoy = 0 ; itoy < ntoys ; itoy++){

//         params_null->assignValueOnly(preParams_null);
//         params_test->assignValueOnly(preParams_test);
//   	RooDataHist *binnedtoy = pdfNull->generateBinned(RooArgSet(*mass),ndata,0,1);

// 	int stat_n=1;
//         int stat_t=1;
// 	int ntries = 0;
// 	double nllNull,nllTest;
// 	// Iterate on the fit 
// 	int MaxTries = 2;
// 	while (stat_n!=0){
// 	  if (ntries>=MaxTries) break;
// 	  RooFitResult *fitNull = pdfNull->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1),RooFit::SumW2Error(kTRUE) //FIXME
// 		,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1));
// 		//,RooFit::Optimize(0));

// 	  nllNull = fitNull->minNll();
//           stat_n = fitNull->status();
// 	  if (stat_n!=0) params_null->assignValueOnly(fitNullData->randomizePars());
// 	  ntries++; 
// 	}
	
// 	ntries = 0;
// 	while (stat_t!=0){
// 	  if (ntries>=MaxTries) break;
// 	  RooFitResult *fitTest = pdfTest->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1),RooFit::SumW2Error(kTRUE) //FIXME
// 		,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1));
// 	  nllTest = fitTest->minNll();
//           stat_t = fitTest->status();
// 	  if (stat_t!=0) params_test->assignValueOnly(fitTestData->randomizePars()); 
// 	  ntries++; 
// 	}
       
// 	toyhistStatN.Fill(stat_n);
// 	toyhistStatT.Fill(stat_t);

//         if (stat_t !=0 || stat_n !=0) continue;
// 	nsuccesst++;
// 	double chi2_t = 2*(nllNull-nllTest);
// 	if (chi2_t >= chi2) npass++;
//         toyhist.Fill(chi2_t);
//   }

//   double prob=0;
//   if (nsuccesst!=0)  prob = (double)npass / nsuccesst;
//   toyhist.Scale(1./(w*toyhist.Integral()));
//   toyhist.Draw();
//   TArrow lData(chi2,toyhist.GetMaximum(),chi2,0);
//   lData.SetLineWidth(2);
//   lData.Draw();
//   gChi2->Draw("L");
//   TLatex *lat = new TLatex();
//   lat->SetNDC();
//   lat->SetTextFont(42);
//   lat->DrawLatex(0.1,0.91,Form("Prob (asymptotic) = %.4f (%.4f)",prob,prob_asym));
//   can->SaveAs(name.c_str());

//   TCanvas *stas =new TCanvas();
//   toyhistStatN.SetLineColor(2);
//   toyhistStatT.SetLineColor(1); 
//   TLegend *leg = new TLegend(0.2,0.6,0.4,0.87); leg->SetFillColor(0);
//   leg->SetTextFont(42);
//   leg->AddEntry(&toyhistStatN,"Null Hyp","L");
//   leg->AddEntry(&toyhistStatT,"Test Hyp","L");
//   toyhistStatN.Draw();
//   toyhistStatT.Draw("same");
//   leg->Draw();
//   stas->SaveAs(Form("%s_fitstatus.pdf",name.c_str()));
//   //reassign params
//   params_null->assignValueOnly(preParams_null);
//   params_test->assignValueOnly(preParams_test);

//   delete can; delete stas;
//   delete gChi2;
//   delete leg;
//   delete lat;

//   // Still return the asymptotic prob (usually its close to the toys one)
//   return prob_asym;

// }

// void plot(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, string name,vector<string> flashggCats_, int status, double *prob){
  
//   // Chi2 taken from full range fit
//   RooPlot *plot_chi2 = mass->frame();
//   data->plotOn(plot_chi2,Binning(nBinsForMass));
//   pdf->plotOn(plot_chi2);

//   int np = pdf->getParameters(*data)->getSize()+1; //Because this pdf has no extend
//   double chi2 = plot_chi2->chiSquare(np);
 
//   *prob = getGoodnessOfFit(mass,pdf,data,name);
//   RooPlot *plot = mass->frame();
//   mass->setRange("unblindReg_1",mgg_low,120);
//   mass->setRange("unblindReg_2",130,mgg_high);
//   if (BLIND) {
//     data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_1"));
//     data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_2"));
//     data->plotOn(plot,Binning(mgg_high-mgg_low),Invisible());
//   }
//   else data->plotOn(plot,Binning(mgg_high-mgg_low));

//  // data->plotOn(plot,Binning(mgg_high-mgg_low));
//   TCanvas *canv = new TCanvas();
//   pdf->plotOn(plot);//,RooFit::NormRange("fitdata_1,fitdata_2"));
//   pdf->paramOn(plot, RooFit::Layout(0.34, 0.96, 0.89), RooFit::Format("NEA", AutoPrecision(1)));
//   plot->getAttText()->SetTextSize(0.025);
//   if (BLIND) plot->SetMinimum(0.0001);
//   plot->SetTitle("");
//   plot->Draw();
//   TLatex *lat = new TLatex();
//   lat->SetNDC();
//   lat->SetTextFont(42);
//   lat->DrawLatex(0.1,0.92,Form("#chi^{2} = %.3f, Prob = %.2f, Fit Status = %d ",chi2*(nBinsForMass-np),*prob,status));
//   canv->SaveAs(name.c_str());
 	
// 	//plot_chi2->Draw();
//   //canv->SaveAs((name+"debug").c_str());

//   delete canv;
//   delete lat;
// }
// void plot(RooRealVar *mass, RooMultiPdf *pdfs, RooCategory *catIndex, RooDataSet *data, string name, vector<string> flashggCats_, int cat, int bestFitPdf=-1){
  
//   int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
//   TLegend *leg = new TLegend(0.5,0.55,0.92,0.88);
//   leg->SetFillColor(0);
//   leg->SetLineColor(1);
//   RooPlot *plot = mass->frame();

//   mass->setRange("unblindReg_1",mgg_low,120);
//   mass->setRange("unblindReg_2",130,mgg_high);
//   if (BLIND) {
//     data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_1"));
//     data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_2"));
//     data->plotOn(plot,Binning(mgg_high-mgg_low),Invisible());
//   }
//   else data->plotOn(plot,Binning(mgg_high-mgg_low)); 
//   TCanvas *canv = new TCanvas();
//   ///start extra bit for ratio plot///
//   RooHist *plotdata = (RooHist*)plot->getObject(plot->numItems()-1);
//   bool doRatioPlot_=1;
//   TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
//   TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.35);
//   pad1->SetBottomMargin(0.18);
//   pad2->SetTopMargin(0.00001);
//   pad2->SetBottomMargin(0.25);
//   pad1->Draw();
//   pad2->Draw();
//   pad1->cd();
//   // enf extra bit for ratio plot///

//   int currentIndex = catIndex->getIndex();
//   TObject *datLeg = plot->getObject(int(plot->numItems()-1));
//   leg->AddEntry(datLeg,Form("Data - %s",flashggCats_[cat].c_str()),"LEP");
//   int style=1;
//   RooAbsPdf *pdf;
//   RooCurve *nomBkgCurve;
//   int bestcol= -1;
//   for (int icat=0;icat<catIndex->numTypes();icat++){
//     int col;
//     if (icat<=6) col=color[icat];
//     else {col=kBlack; style++;}
//     catIndex->setIndex(icat);
//     pdfs->getCurrentPdf()->fitTo(*data,RooFit::Minos(0),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE));	 //FIXME
//     pdfs->getCurrentPdf()->plotOn(plot,LineColor(col),LineStyle(style));//,RooFit::NormRange("fitdata_1,fitdata_2"));
//     TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
//     std::string ext = "";
//     if (bestFitPdf==icat) {
//     ext=" (Best Fit Pdf) ";
//     pdf= pdfs->getCurrentPdf();
//     nomBkgCurve = (RooCurve*)plot->getObject(plot->numItems()-1);
//     bestcol = col;
//     }
//     leg->AddEntry(pdfLeg,Form("%s%s",pdfs->getCurrentPdf()->GetName(),ext.c_str()),"L");
//   }
//   plot->SetTitle(Form("Category %s",flashggCats_[cat].c_str()));
//   if (BLIND) plot->SetMinimum(0.0001);
//   plot->Draw();
//   leg->Draw("same");
//   CMS_lumi( canv, 0, 0);
//   ///start extra bit for ratio plot///
//   TH1D *hbplottmp = (TH1D*) pdf->createHistogram("hbplottmp",*mass,Binning(mgg_high-mgg_low,mgg_low,mgg_high));
//   hbplottmp->Scale(plotdata->Integral());
//   hbplottmp->Draw("same");
//   int npoints = plotdata->GetN();
//   double xtmp,ytmp;//
//   int point =0;
//   TGraphAsymmErrors *hdatasub = new TGraphAsymmErrors(npoints);
//   //hdatasub->SetMarkerSize(defmarkersize);
//   for (int ipoint=0; ipoint<npoints; ++ipoint) {
//   //double bkgval = hbplottmp->GetBinContent(ipoint+1);
//   plotdata->GetPoint(ipoint, xtmp,ytmp);
//   double bkgval = nomBkgCurve->interpolate(xtmp);
//   if (BLIND) {
//    if ((xtmp > 120 ) && ( xtmp < 130) ) continue;
//   }
//   std::cout << "[INFO] plotdata->Integral() " <<  plotdata->Integral() << " ( bins " << npoints  << ") hbkgplots[i]->Integral() " << hbplottmp->Integral() << " (bins " << hbplottmp->GetNbinsX() << std::endl;
//  double errhi = plotdata->GetErrorYhigh(ipoint);
//  double errlow = plotdata->GetErrorYlow(ipoint);
       
//  //std::cout << "[INFO]  Channel " << name  << " errhi " << errhi << " errlow " << errlow  << std::endl;
//  std::cout << "[INFO] Channel  " << name << " setting point " << point <<" : xtmp "<< xtmp << "  ytmp " << ytmp << " bkgval  " << bkgval << " ytmp-bkgval " << ytmp-bkgval << std::endl;
//  bool drawZeroBins_ =1;
//  if (!drawZeroBins_) if(fabs(ytmp)<1e-5) continue; 
//  hdatasub->SetPoint(point,xtmp,ytmp-bkgval);
//  hdatasub->SetPointError(point,0.,0.,errlow,errhi );
//  point++;
//   } 
//   pad2->cd();
//   TH1 *hdummy = new TH1D("hdummyweight","",mgg_high-mgg_low,mgg_low,mgg_high);
//   hdummy->SetMaximum(hdatasub->GetHistogram()->GetMaximum()+1);
//   hdummy->SetMinimum(hdatasub->GetHistogram()->GetMinimum()-1);
//   hdummy->GetYaxis()->SetTitle("data - best fit PDF");
//   hdummy->GetYaxis()->SetTitleSize(0.12);
//   hdummy->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
//   hdummy->GetXaxis()->SetTitleSize(0.12);
//   hdummy->Draw("HIST");
//   hdummy->GetYaxis()->SetNdivisions(808);

//   TLine *line3 = new TLine(mgg_low,0.,mgg_high,0.);
//   line3->SetLineColor(bestcol);
//   //line3->SetLineStyle(kDashed);
//   line3->SetLineWidth(5.0);
//   line3->Draw();
//   hdatasub->Draw("PESAME");
//   // enf extra bit for ratio plot///
//   canv->SaveAs(Form("%s.pdf",name.c_str()));
//   canv->SaveAs(Form("%s.png",name.c_str()));
//   catIndex->setIndex(currentIndex);
//   delete canv;
// }

// void plot(RooRealVar *mass, map<string,RooAbsPdf*> pdfs, RooDataSet *data, string name, vector<string> flashggCats_, int cat, int bestFitPdf=-1){
  
//   int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
//   TCanvas *canv = new TCanvas();
//   TLegend *leg = new TLegend(0.6,0.65,0.88,0.88);
//   leg->SetFillColor(0);
//   leg->SetLineColor(0);
//   RooPlot *plot = mass->frame();

//   mass->setRange("unblindReg_1",mgg_low,120);
//   mass->setRange("unblindReg_2",130,mgg_high);
//   if (BLIND) {
//     data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_1"));
//     data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_2"));
//     data->plotOn(plot,Binning(mgg_high-mgg_low),Invisible());
//   }
//   else data->plotOn(plot,Binning(mgg_high-mgg_low));

//   TObject *datLeg = plot->getObject(int(plot->numItems()-1));
// 	if(flashggCats_.size() >0){
//   leg->AddEntry(datLeg,Form("Data - %s",flashggCats_[cat].c_str()),"LEP");
// 	} else {
//   leg->AddEntry(datLeg,Form("Data - %d",cat),"LEP");
// 	}
//   int i=0;
//   int style=1;
//   for (map<string,RooAbsPdf*>::iterator it=pdfs.begin(); it!=pdfs.end(); it++){
//     int col;
//     if (i<=6) col=color[i];
//     else {col=kBlack; style++;}
//     it->second->plotOn(plot,LineColor(col),LineStyle(style));//,RooFit::NormRange("fitdata_1,fitdata_2"));
//     TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
//     std::string ext = "";
//     if (bestFitPdf==i) ext=" (Best Fit Pdf) ";
//     leg->AddEntry(pdfLeg,Form("%s%s",it->first.c_str(),ext.c_str()),"L");
//     i++;
//   }
//   plot->SetTitle(Form(" %s",flashggCats_[cat].c_str()));
//   if (BLIND) plot->SetMinimum(0.0001);
//   plot->Draw();
//   leg->Draw("same");
//   CMS_lumi( canv, 0, 0);
//   canv->SaveAs(Form("%s.pdf",name.c_str()));
//   canv->SaveAs(Form("%s.png",name.c_str()));
//   delete canv;
// }

// void transferMacros(TFile *inFile, TFile *outFile){
  
//   TIter next(inFile->GetListOfKeys());
//   TKey *key;
//   while ((key = (TKey*)next())){
//     if (string(key->ReadObj()->ClassName())=="TMacro") {
//       //cout << key->ReadObj()->ClassName() << " : " << key->GetName() << endl;
//       TMacro *macro = (TMacro*)inFile->Get(key->GetName());
//       outFile->cd();
//       macro->Write();
//     }
//   }
// }
// int getBestFitFunction(RooMultiPdf *bkg, RooDataSet *data, RooCategory *cat, bool silent=false){


// 	double global_minNll = 1E10;
// 	int best_index = 0;
// 	int number_of_indeces = cat->numTypes();
		
// 	RooArgSet snap,clean;
// 	RooArgSet *params = bkg->getParameters((const RooArgSet*)0);
// 	params->remove(*cat);
// 	params->snapshot(snap);
// 	params->snapshot(clean);
// 	if (!silent) {
// 		//params->Print("V");
// 	}
	
// 	//bkg->setDirtyInhibit(1);
// 	//RooAbsReal *nllm = bkg->createNLL(*data);
// 	//RooMinimizer minim(*nllm);
// 	//minim.setStrategy(1);
	
// 	for (int id=0;id<number_of_indeces;id++){		
// 		params->assignValueOnly(clean);
// 		cat->setIndex(id);

// 		//RooAbsReal *nllm = bkg->getCurrentPdf()->createNLL(*data);

// 		if (!silent) {
// 			/*
// 			std::cout << "BEFORE  MAKING FIT" << std::endl;
// 			params->Print("V");
// 			std::cout << "-----------------------" << std::endl;		
// 			*/
// 		}
		
// 		//minim.minimize("Minuit2","minimize");
// 		double minNll=0; //(nllm->getVal())+bkg->getCorrection();
// 		int fitStatus=1;		
// 		runFit(bkg->getCurrentPdf(),data,&minNll,&fitStatus,/*max iterations*/3);
// 		// Add the penalty

// 		minNll=minNll+bkg->getCorrection();

// 		if (!silent) {
// 			/*
// 			std::cout << "After Minimization ------------------  " <<std::endl;
// 			std::cout << bkg->getCurrentPdf()->GetName() << " " << minNll <<std::endl;
// 			bkg->Print("v");
// 			bkg->getCurrentPdf()->getParameters(*data)->Print("V");
// 			std::cout << " ------------------------------------  " << std::endl;
	
// 			params->Print("V");
// 			*/
// 			std::cout << "[INFO] AFTER FITTING" << std::endl;
// 			std::cout << "[INFO] Function was " << bkg->getCurrentPdf()->GetName() <<std::endl;
// 			std::cout << "[INFO] Correction Applied is " << bkg->getCorrection() <<std::endl;
// 			std::cout << "[INFO] NLL + c = " <<  minNll << std::endl;
// 			std::cout << "-----------------------" << std::endl;
// 		}
			
// 		if (minNll < global_minNll){
//         		global_minNll = minNll;
// 			snap.assignValueOnly(*params);
//         		best_index=id;
// 		}
// 	}
//     	cat->setIndex(best_index);
// 	params->assignValueOnly(snap);
	
// 	if (!silent) {
// 		std::cout << "[INFO] Best fit Function -- " << bkg->getCurrentPdf()->GetName() << " " << cat->getIndex() <<std::endl;
// 		//bkg->getCurrentPdf()->getParameters(*data)->Print("v");
// 	}
// 	return best_index;
// }

// int main(int argc, char* argv[]){
 
//   setTDRStyle();
//   writeExtraText = true;       // if extra text
//   extraText  = "Preliminary";  // default extra text is "Preliminary"
//   lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
//   lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
//   lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
//   string year_ = "2016";
//   //int year_ = 2017;

//   string fileName;
//   int ncats;
//   int singleCategory;
//   int catOffset;
//   string datfile;
//   string outDir;
//   string outfilename;
//   bool is2011=false;
//   bool verbose=false;
//   bool saveMultiPdf=false;
//   int isFlashgg_ =1;
//   string flashggCatsStr_;
//   vector<string> flashggCats_;
//  bool isData_ =0;

//   po::options_description desc("Allowed options");
//   desc.add_options()
//     ("help,h",                                                                                  "Show help")
//     ("infilename,i", po::value<string>(&fileName),                                              "In file name")
//     ("ncats,c", po::value<int>(&ncats)->default_value(5),                                       "Number of categories")
//     ("singleCat", po::value<int>(&singleCategory)->default_value(-1),                           "Run A single Category")
//     ("datfile,d", po::value<string>(&datfile)->default_value("dat/fTest.dat"),                  "Right results to datfile for BiasStudy")
//     ("outDir,D", po::value<string>(&outDir)->default_value("plots/fTest"),                      "Out directory for plots")
//     ("saveMultiPdf", po::value<string>(&outfilename),         					"Save a MultiPdf model with the appropriate pdfs")
//     ("runFtestCheckWithToys", 									"When running the F-test, use toys to calculate pvals (and make plots) ")
//     ("is2011",                                                                                  "Run 2011 config")
//     ("is2012",                                                                                  "Run 2012 config")
//     ("unblind",  									        "Dont blind plots")
//     ("isFlashgg",  po::value<int>(&isFlashgg_)->default_value(1),  								    	        "Use Flashgg output ")
//     ("isData",  po::value<bool>(&isData_)->default_value(0),  								    	        "Use Data not MC ")
// 		("flashggCats,f", po::value<string>(&flashggCatsStr_)->default_value("UntaggedTag_0,UntaggedTag_1,UntaggedTag_2,UntaggedTag_3,UntaggedTag_4,VBFTag_0,VBFTag_1,VBFTag_2,TTHHadronicTag,TTHLeptonicTag,VHHadronicTag,VHTightTag,VHLooseTag,VHEtTag"),       "Flashgg category names to consider")
//     ("year", po::value<string>(&year_)->default_value("2016"),       "Dataset year")
//     ("catOffset", po::value<int>(&catOffset)->default_value(0),       "Category numbering scheme offset")
//     ("verbose,v",                                                                               "Run with more output")
//   ;
//   po::variables_map vm;
//   po::store(po::parse_command_line(argc,argv,desc),vm);
//   po::notify(vm);
//   if (vm.count("help")) { cout << desc << endl; exit(1); }
//   if (vm.count("is2011")) is2011=true;
// 	if (vm.count("unblind")) BLIND=false;
//   saveMultiPdf = vm.count("saveMultiPdf");

//   if (vm.count("verbose")) verbose=true;
//   if (vm.count("runFtestCheckWithToys")) runFtestCheckWithToys=true;

//   if (!verbose) {
//     RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
//     RooMsgService::instance().setSilentMode(true);
//     gErrorIgnoreLevel=kWarning;
//   }
// 	split(flashggCats_,flashggCatsStr_,boost::is_any_of(","));
  
// 	int startingCategory=0;
//   if (singleCategory >-1){
// 	ncats=singleCategory+1;	
// 	startingCategory=singleCategory;
//   }
// 	if (isFlashgg_==1){
	
// 	ncats= flashggCats_.size();

// 	}

//   if(verbose) std::cout << "[INFO] SaveMultiPdf? " << saveMultiPdf << std::endl;
//   TFile *outputfile;
//   RooWorkspace *outputws;

//   if (saveMultiPdf){
// 	outputfile = new TFile(outfilename.c_str(),"RECREATE");
// 	outputws = new RooWorkspace(); outputws->SetName("multipdf");
//   }

//   system(Form("mkdir -p %s",outDir.c_str()));
//   TFile *inFile = TFile::Open(fileName.c_str());
//   RooWorkspace *inWS;
// 	if(isFlashgg_){
// 		if (isData_){
// 			inWS = (RooWorkspace*)inFile->Get("tagsDumper/cms_hgg_13TeV");
// 		} else {
// 			inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
// 		}
// 	} else {
// 		inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");//FIXME
// 	}
// 	if (verbose) std::cout << "[INFO]  inWS open " << inWS << std::endl;
// 	if (saveMultiPdf){
// 		transferMacros(inFile,outputfile);

// 		RooRealVar *intL; 
// 		RooRealVar *sqrts;

// 		if (isFlashgg_){
// 			//intL  = (RooRealVar*)inWS->var("IntLumi");
// 			intL  = intLumi_;
// 			sqrts = (RooRealVar*)inWS->var("SqrtS");
// 			if (!sqrts){ sqrts = new RooRealVar("SqrtS","SqrtS",13); }
// 		std::cout << "[INFO] got intL and sqrts " << intL << ", " << sqrts << std::endl;


// 		} else {
// 			//intL  = (RooRealVar*)inWS->var("IntLumi");
// 			intL  = intLumi_;
// 			sqrts = (RooRealVar*)inWS->var("Sqrts");
// 		}
// 		outputws->import(*intL);
// 		outputws->import(*sqrts);
// 		std::cout << "[INFO] got intL and sqrts " << intL << ", " << sqrts << std::endl;
// 	}

// 	vector<string> functionClasses;
// 	// functionClasses.push_back("Bernstein");
// 	// functionClasses.push_back("Exponential");
// 	// functionClasses.push_back("PowerLaw");
// 	// functionClasses.push_back("Laurent");
//   functionClasses.push_back("LaurentStepxGau");
//   functionClasses.push_back("BernsteinStepxGau");
//   functionClasses.push_back("ExponentialStepxGau");
//   functionClasses.push_back("PowerLawStepxGau");
// 	map<string,string> namingMap;
// 	namingMap.insert(pair<string,string>("Bernstein","pol"));
// 	namingMap.insert(pair<string,string>("Exponential","exp"));
// 	namingMap.insert(pair<string,string>("PowerLaw","pow"));
// 	namingMap.insert(pair<string,string>("Laurent","lau"));
//   namingMap.insert(pair<string,string>("BernsteinStepxGau","bern"));
//   namingMap.insert(pair<string,string>("ExponentialStepxGau","exp"));
//   namingMap.insert(pair<string,string>("PowerLawStepxGau","pow"));
//   namingMap.insert(pair<string,string>("LaurentStepxGau","lau"));

// 	// store results here

// 	FILE *resFile ;
// 	if  (singleCategory >-1) resFile = fopen(Form("%s/fTestResults_%s.txt",outDir.c_str(),flashggCats_[singleCategory].c_str()),"w");
// 	else resFile = fopen(Form("%s/fTestResults.txt",outDir.c_str()),"w");
// 	vector<map<string,int> > choices_vec;
// 	vector<map<string,std::vector<int> > > choices_envelope_vec;
// 	vector<map<string,RooAbsPdf*> > pdfs_vec;

// 	PdfModelBuilder pdfsModel;
// 	RooRealVar *mass = (RooRealVar*)inWS->var("CMS_hgg_mass");
// 	std:: cout << "[INFO] Got mass from ws " << mass << std::endl;
// 	pdfsModel.setObsVar(mass);
// 	double upperEnvThreshold = 0.1; // upper threshold on delta(chi2) to include function in envelope (looser than truth function)

// 	fprintf(resFile,"Truth Model & d.o.f & $\\Delta NLL_{N+1}$ & $p(\\chi^{2}>\\chi^{2}_{(N\\rightarrow N+1)})$ \\\\\n");
// 	fprintf(resFile,"\\hline\n");

// 	std::string ext = is2011 ? "7TeV" : "8TeV";
//         if( isFlashgg_ ){
//           if( year_ == "all" ){ ext = "13TeV"; }
//           //else{ ext = "13TeV"; } //FIXME 
//           else{ ext = Form("%s_13TeV",year_.c_str()); }
//         }
// 	//if (isFlashgg_) ext = "13TeV";
//         //FIXME trying to remove duplicated names for 2016+2017 combination
// 	//if (isFlashgg_) ext = Form("13TeV_%d",year_);
// 	for (int cat=startingCategory; cat<ncats; cat++){

// 		map<string,int> choices;
// 		map<string,std::vector<int> > choices_envelope;
// 		map<string,RooAbsPdf*> pdfs;
// 		map<string,RooAbsPdf*> allPdfs;
// 		string catname;
// 		if (isFlashgg_){
// 			catname = Form("%s",flashggCats_[cat].c_str());
// 		} else {
// 			catname = Form("cat%d",cat);
// 		}
// 		RooDataSet *dataFull;
// 		RooDataSet *dataFull0;
// 		if (isData_) {
//     dataFull = (RooDataSet*)inWS->data(Form("Data_13TeV_%s",catname.c_str()));
//     /*dataFull= (RooDataSet*) dataFull0->emptyClone();
//     for (int i =0 ; i < dataFull0->numEntries() ; i++){
//     double m = dataFull0->get(i)->getRealValue("CMS_hgg_mass");
//     //if (m <(mgg_low+0.01) or m > (mgg_high-0.01)) 

//     if (m==mgg_low){
//     std::cout << "dataset mass m="<< m << std::endl;
//     continue;
//     }
//     dataFull->add(*dataFull0->get(),1.0);
//     }*/
// 		if (verbose) std::cout << "[INFO] opened data for  "  << Form("Data_%s",catname.c_str()) <<" - " << dataFull <<std::endl;
//     }
// 		else 
//     {dataFull = (RooDataSet*)inWS->data(Form("data_mass_%s",catname.c_str()));
// 		if (verbose) std::cout << "[INFO] opened data for  "  << Form("data_mass_%s",catname.c_str()) <<" - " << dataFull <<std::endl;
//     }


// 		mass->setBins(nBinsForMass);
// 		RooDataSet *data;
// 		//	RooDataHist thisdataBinned(Form("roohist_data_mass_cat%d",cat),"data",*mass,*dataFull);
// 		//	RooDataSet *data = (RooDataSet*)&thisdataBinned;
// 		string thisdataBinned_name;

// 		if ( isFlashgg_){
// 			thisdataBinned_name =Form("roohist_data_mass_%s",flashggCats_[cat].c_str());
// 			//	RooDataHist thisdataBinned(Form("roohist_data_mass_cat%d",cat),"data",*mass,*dataFull);
// 			//	data = (RooDataSet*)&thisdataBinned;
// 			//		std::cout << "debug " << thisdataBinned.GetName() << std::endl;

// 			//RooDataSet *data = (RooDataSet*)dataFull;
// 		} else {
// 			thisdataBinned_name= Form("roohist_data_mass_cat%d",cat);
// 			//RooDataSet *data = (RooDataSet*)dataFull;
// 		}
// 		RooDataHist thisdataBinned(thisdataBinned_name.c_str(),"data",*mass,*dataFull);
// 		data = (RooDataSet*)&thisdataBinned;

// 		RooArgList storedPdfs("store");

// 		fprintf(resFile,"\\multicolumn{4}{|c|}{\\textbf{Category %d}} \\\\\n",cat);
// 		fprintf(resFile,"\\hline\n");

// 		double MinimimNLLSoFar=1e10;
// 		int simplebestFitPdfIndex = 0;

// 		// Standard F-Test to find the truth functions
// 		for (vector<string>::iterator funcType=functionClasses.begin(); 
// 				funcType!=functionClasses.end(); funcType++){

// 			double thisNll=0.; double prevNll=0.; double chi2=0.; double prob=0.; 
// 			int order=1; int prev_order=0; int cache_order=0;

// 			RooAbsPdf *prev_pdf=NULL;
// 			RooAbsPdf *cache_pdf=NULL;
// 			std::vector<int> pdforders;

// 			int counter =0;
// 			//	while (prob<0.05){
// 			while (prob<0.05 && order < 7){ //FIXME
// 				RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order,Form("ftest_pdf_%d_%s",(cat+catOffset),ext.c_str()));
// 				if (!bkgPdf){
// 					// assume this order is not allowed
// 					order++;
// 				}
// 				else {

// 					//RooFitResult *fitRes = bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));
// 					int fitStatus = 0;
// 					//thisNll = fitRes->minNll();
//           bkgPdf->Print();
// 					runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/3);//bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));
// 					if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;       
// 					chi2 = 2.*(prevNll-thisNll);
// 					if (chi2<0. && order>1) chi2=0.;
// 					if (prev_pdf!=NULL){
// 						prob = getProbabilityFtest(chi2,order-prev_order,prev_pdf,bkgPdf,mass,data
// 								,Form("%s/Ftest_from_%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,(cat+catOffset)));
// 						std::cout << "[INFO]  F-test Prob(chi2>chi2(data)) == " << prob << std::endl;
// 					} else {
// 						prob = 0;
// 					}
// 					double gofProb=0;
// 					// otherwise we get it later ...
// 					if (!saveMultiPdf) plot(mass,bkgPdf,data,Form("%s/%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,(cat+catOffset)),flashggCats_,fitStatus,&gofProb);
// 					cout << "[INFO]\t " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
// 					//fprintf(resFile,"%15s && %d && %10.2f && %10.2f && %10.2f \\\\\n",funcType->c_str(),order,thisNll,chi2,prob);
// 					prevNll=thisNll;
// 					cache_order=prev_order;
// 					cache_pdf=prev_pdf;
// 					prev_order=order;
// 					prev_pdf=bkgPdf;
//           cout << "[INFO] \t " << *funcType << " " << cache_order << " " << prev_order << endl;
// 					order++;
// 				}
// 				counter++;
// 			}

// 			fprintf(resFile,"%15s & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
// 			choices.insert(pair<string,int>(*funcType,cache_order));
// 			pdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),cache_order),cache_pdf));

// 			int truthOrder = cache_order;

// 			// Now run loop to determine functions inside envelope
// 			if (saveMultiPdf){
// 				chi2=0.;
// 				thisNll=0.;
// 				prevNll=0.;
// 				prob=0.;
// 				order=1;
// 				prev_order=0;
// 				cache_order=0;
// 				std::cout << "[INFO] Determining Envelope Functions for Family " << *funcType << ", cat " << cat << std::endl;
// 				std::cout << "[INFO] Upper end Threshold for highest order function " << upperEnvThreshold <<std::endl;

// 				while (prob<upperEnvThreshold){
// 					RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order,Form("env_pdf_%d_%s",(cat+catOffset),ext.c_str()));
// 					if (!bkgPdf ){
// 						// assume this order is not allowed
// 						if (order >6) { std::cout << " [WARNING] could not add ] " << std::endl; break ;}
// 						order++;
// 					}
// 					else {
// 						//RooFitResult *fitRes;
// 						int fitStatus=0;
// 						runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/3);//bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));
// 						//thisNll = fitRes->minNll();
// 						if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;
// 						double myNll = 2.*thisNll;
// 						chi2 = 2.*(prevNll-thisNll);
// 						if (chi2<0. && order>1) chi2=0.;
// 						prob = TMath::Prob(chi2,order-prev_order);

// 						cout << "[INFO] \t " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
// 						prevNll=thisNll;
// 						cache_order=prev_order;
// 						cache_pdf=prev_pdf;

// 						// Calculate goodness of fit for the thing to be included (will use toys for lowstats)!
// 						double gofProb =0; 
// 						plot(mass,bkgPdf,data,Form("%s/%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,(cat+catOffset)),flashggCats_,fitStatus,&gofProb);

// 						if ((prob < upperEnvThreshold) ) { // Looser requirements for the envelope

// 							if (gofProb > 0.01 || order == truthOrder ) {  // Good looking fit or one of our regular truth functions

// 								std::cout << "[INFO] Adding to Envelope " << bkgPdf->GetName() << " "<< gofProb 
// 									<< " 2xNLL + c is " << myNll + bkgPdf->getVariables()->getSize() << " truth order " << truthOrder << std::endl;
// 								allPdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),order),bkgPdf));
// 								storedPdfs.add(*bkgPdf);
// 								pdforders.push_back(order);

// 								// Keep track but we shall redo this later
// 								if ((myNll + bkgPdf->getVariables()->getSize()) < MinimimNLLSoFar) {
// 									simplebestFitPdfIndex = storedPdfs.getSize()-1;
// 									MinimimNLLSoFar = myNll + bkgPdf->getVariables()->getSize();
// 								}
// 							}
// 						}

// 						prev_order=order;
// 						prev_pdf=bkgPdf;
// 						order++;
// 					}
// 				}

// 				fprintf(resFile,"%15s & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
// 				choices_envelope.insert(pair<string,std::vector<int> >(*funcType,pdforders));
// 			}
// 		}

// 		fprintf(resFile,"\\hline\n");
// 		choices_vec.push_back(choices);
// 		choices_envelope_vec.push_back(choices_envelope);
// 		pdfs_vec.push_back(pdfs);

// 		plot(mass,pdfs,data,Form("%s/truths_cat%d",outDir.c_str(),(cat+catOffset)),flashggCats_,cat);

// 		if (saveMultiPdf){


// 			// Put selectedModels into a MultiPdf
// 			string catindexname;
// 			string catname;
// 			if (isFlashgg_){
// 				catindexname = Form("pdfindex_%s_%s",flashggCats_[cat].c_str(),ext.c_str());
// 				catname = Form("%s",flashggCats_[cat].c_str());
// 			} else {
// 				catindexname = Form("pdfindex_%d_%s",(cat+catOffset),ext.c_str());
// 				catname = Form("cat%d",(cat+catOffset));
// 			}
// 			RooCategory catIndex(catindexname.c_str(),"c");
// 			RooMultiPdf *pdf = new RooMultiPdf(Form("CMS_hgg_%s_%s_bkgshape",catname.c_str(),ext.c_str()),"all pdfs",catIndex,storedPdfs);
// 			//RooRealVar nBackground(Form("CMS_hgg_%s_%s_bkgshape_norm",catname.c_str(),ext.c_str()),"nbkg",data->sumEntries(),0,10E8);
// 			RooRealVar nBackground(Form("CMS_hgg_%s_%s_bkgshape_norm",catname.c_str(),ext.c_str()),"nbkg",data->sumEntries(),0,3*data->sumEntries());
// 			//nBackground.removeRange(); // bug in roofit will break combine until dev branch brought in
// 			//double check the best pdf!
// 			int bestFitPdfIndex = getBestFitFunction(pdf,data,&catIndex,!verbose);
// 			catIndex.setIndex(bestFitPdfIndex);
// 			std::cout << "// ------------------------------------------------------------------------- //" <<std::endl; 
// 			std::cout << "[INFO] Created MultiPdf " << pdf->GetName() << ", in Category " << cat << " with a total of " << catIndex.numTypes() << " pdfs"<< std::endl;
// 			storedPdfs.Print();
// 			std::cout << "[INFO] Best Fit Pdf = " << bestFitPdfIndex << ", " << storedPdfs.at(bestFitPdfIndex)->GetName() << std::endl;
// 			std::cout << "// ------------------------------------------------------------------------- //" <<std::endl;
// 			std::cout << "[INFO] Simple check of index "<< simplebestFitPdfIndex <<std::endl;

// 			mass->setBins(nBinsForMass);
// 			RooDataHist dataBinned(Form("roohist_data_mass_%s",catname.c_str()),"data",*mass,*dataFull);

// 			// Save it (also a binned version of the dataset
// 			outputws->import(*pdf);
// 			outputws->import(nBackground);
// 			outputws->import(catIndex);
// 			outputws->import(dataBinned);
// 			outputws->import(*data);
// 			plot(mass,pdf,&catIndex,data,Form("%s/multipdf_%s",outDir.c_str(),catname.c_str()),flashggCats_,cat,bestFitPdfIndex);

// 		}

// 		}
// 		if (saveMultiPdf){
// 			outputfile->cd();
// 			outputws->Write();
// 			outputfile->Close();	
// 		}

// 		FILE *dfile = fopen(datfile.c_str(),"w");
// 		cout << "[RESULT] Recommended options" << endl;

// 		for (int cat=startingCategory; cat<ncats; cat++){
// 			cout << "Cat " << cat << endl;
// 			fprintf(dfile,"cat=%d\n",(cat+catOffset)); 
// 			for (map<string,int>::iterator it=choices_vec[cat-startingCategory].begin(); it!=choices_vec[cat-startingCategory].end(); it++){
// 				cout << "\t" << it->first << " - " << it->second << endl;
// 				fprintf(dfile,"truth=%s:%d:%s%d\n",it->first.c_str(),it->second,namingMap[it->first].c_str(),it->second);
// 			}
// 			for (map<string,std::vector<int> >::iterator it=choices_envelope_vec[cat-startingCategory].begin(); it!=choices_envelope_vec[cat-startingCategory].end(); it++){
// 				std::vector<int> ords = it->second;
// 				for (std::vector<int>::iterator ordit=ords.begin(); ordit!=ords.end(); ordit++){
// 					fprintf(dfile,"paul=%s:%d:%s%d\n",it->first.c_str(),*ordit,namingMap[it->first].c_str(),*ordit);
// 				}
// 			}
// 			fprintf(dfile,"\n");
// 		}
// 		inFile->Close();

// 		return 0;
// 	}
