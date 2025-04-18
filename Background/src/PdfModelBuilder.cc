#include "TCanvas.h"

#include "RooPlot.h"
#include "RooBernstein.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
#include "RooExponential.h"
#include "../interface/RooPowerLaw.h"
#include "../interface/RooPowerLawSum.h"
#include "RooKeysPdf.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooFitResult.h"
#include "RooRandom.h"

#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"

#include "../interface/PdfModelBuilder.h"

#include "HiggsAnalysis/GBRLikelihood/interface/RooDoubleCBFast.h"
#include "RooGaussModel.h"
#include "RooFFTConvPdf.h"

#include "HiggsAnalysis/CombinedLimit/interface/HGGRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZGRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooBernsteinFast.h"

using namespace std;
using namespace RooFit;
using namespace boost;

PdfModelBuilder::PdfModelBuilder():
  obs_var_set(false),
  signal_modifier_set(false),
  signal_set(false),
  bkgHasFit(false),
  sbHasFit(false),
  keysPdfAttributesSet(false),
  verbosity(0)
{
  
  recognisedPdfTypes.push_back("Bernstein");
  recognisedPdfTypes.push_back("Exponential");
  recognisedPdfTypes.push_back("PowerLaw");
  recognisedPdfTypes.push_back("Laurent");
  recognisedPdfTypes.push_back("KeysPdf");
  recognisedPdfTypes.push_back("File");
  recognisedPdfTypes.push_back("ExpModGaussian");

  wsCache = new RooWorkspace("PdfModelBuilderCache");

};

PdfModelBuilder::~PdfModelBuilder(){};

void PdfModelBuilder::setObsVar(RooRealVar *var){
  obs_var=var;
  obs_var_set=true;
}

void PdfModelBuilder::setSignalModifier(RooRealVar *var){
  signalModifier=var;
  signal_modifier_set=true;
}

void PdfModelBuilder::setSignalModifierVal(float val){
  signalModifier->setVal(val);
}

void PdfModelBuilder::setSignalModifierConstant(bool val){
  signalModifier->setConstant(val);
}

RooAbsPdf* PdfModelBuilder::getDoubleCB(){
  RooRealVar *mean = new RooRealVar("mean", "mean", 125, 120, 130);
  RooRealVar *sigma = new RooRealVar("sigma", "sigma", 2, 0.1, 10);
  RooRealVar *alpha1 = new RooRealVar("alpha1", "alpha1", 1, 0.1, 10);
  RooRealVar *n1 = new RooRealVar("n1", "n1", 2, 0.1, 10);
  RooRealVar *alpha2 = new RooRealVar("alpha2", "alpha2", 1, 0.1, 10);
  RooRealVar *n2 = new RooRealVar("n2", "n2", 2, 0.1, 10);

  RooDoubleCBFast *doubleCB = new RooDoubleCBFast("doubleCB", "Double Sided Crystal Ball", *obs_var, *mean, *sigma, *alpha1, *n1, *alpha2, *n2);
  return doubleCB;
}

RooAbsPdf* PdfModelBuilder::getChebychev(string prefix, int order){
  
  RooArgList *coeffList = new RooArgList();
  for (int i=0; i<order; i++){
    string name = Form("%s_p%d",prefix.c_str(),i);
    //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
    //RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0*@0",RooArgList(*param));
    params.insert(pair<string,RooRealVar*>(name,param));
    //prods.insert(pair<string,RooFormulaVar*>(name,form));
    coeffList->add(*params[name]);
  }
  //RooChebychev *cheb = new RooChebychev(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  RooPolynomial *cheb = new RooPolynomial(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  return cheb;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(bern->GetName(),bern));

}

RooAbsPdf *PdfModelBuilder::getBernsteinStepxGau(string prefix, int order){
  RooRealVar *mean = new RooRealVar(Form("%s_mean", prefix.c_str()), Form("%s_mean", prefix.c_str()), 0.);
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()), 5., 1., 10.);
  RooRealVar *step = new RooRealVar(Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()), 105., 95., 120.);
  RooArgList *coeffList = new RooArgList();
  for (int i = 0; i < order; i++) {
    string name = Form("%s_b%d", prefix.c_str(), i);
    RooRealVar *param = new RooRealVar(name.c_str(), name.c_str(), 0.1 * (i + 1), -15., 15.);
    params.insert(pair<string, RooRealVar*>(name, param));
    coeffList->add(*params[name]);
  }
  if (order < 7) {
    return new RooGaussStepBernstein(prefix.c_str(), prefix.c_str(), *obs_var, *mean, *sigma, *step, *coeffList);
  }
  else {
    return NULL;
  }
}

RooAbsPdf *PdfModelBuilder::getPowerLawStepxGau(string prefix, int order){
  if (order%2==0 || order > 7) return NULL;
  RooRealVar *mean = new RooRealVar(Form("%s_mean", prefix.c_str()), Form("%s_mean", prefix.c_str()), 0.);
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()), 5., 1., 10.);
  RooRealVar *step = new RooRealVar(Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()), 105., 95., 120.);
  RooArgList *coeffList = new RooArgList();
  RooArgList formulaArgs;
  string formula = "TMath::Min(TMath::Max((@0-@1)*153.85, 0.0), 1.0)*(";

  formulaArgs.add(*obs_var);
  formulaArgs.add(*step);
  
  int xmax = obs_var->getMax();

  for (int i = 0; i < order; i += 2) {
    string pName = Form("%s_p%d", prefix.c_str(), i + 1, order);
    string cpName = Form("%s_cp%d", prefix.c_str(), i + 1, order);
    RooRealVar *p = new RooRealVar(pName.c_str(), pName.c_str(), -2.0, -10.0, 0.0);
    RooRealVar *cp = new RooRealVar(cpName.c_str(), cpName.c_str(), 0.1 * (i + 1), 0., 1.);
    formulaArgs.add(*p);
    formulaArgs.add(*cp);
    formula += Form("@%d*@0^(@%d)*(1+@%d)/(%d^(1+@%d)-@1^(1+@%d))", formulaArgs.getSize() - 1, formulaArgs.getSize() - 2, formulaArgs.getSize() - 2, xmax, formulaArgs.getSize() - 2, formulaArgs.getSize() - 2);
    if (i + 2 < order) {
      formula += "+";
    }
  }
  formula += ")";

  RooGenericPdf *stepPdf = new RooGenericPdf(Form("%s_step_pow%d", prefix.c_str(), order), Form("%s_step_pow%d", prefix.c_str(), order), formula.c_str(), formulaArgs);
  RooGaussModel *gau = new RooGaussModel(Form("%s_gau_pow%d", prefix.c_str(), order), Form("%s_gau_pow%d", prefix.c_str(), order), *obs_var, *mean, *sigma);
  RooFFTConvPdf *gauxpow = new RooFFTConvPdf(Form("%s_gauxpow%d", prefix.c_str(), order), Form("%s_gauxpow%d", prefix.c_str(), order), *obs_var, *stepPdf, *gau);
  return gauxpow;
}

RooAbsPdf *PdfModelBuilder::getExponentialStepxGau(string prefix, int order){
  if (order%2==0 || order > 7) return NULL;
  RooRealVar *mean = new RooRealVar(Form("%s_mean", prefix.c_str()), Form("%s_mean", prefix.c_str()), 0.);
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()), 5., 1., 10.);
  RooRealVar *step = new RooRealVar(Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()), 105., 95., 120.);
  RooArgList *coeffList = new RooArgList();
  RooArgList formulaArgs;
  string formula = "TMath::Min(TMath::Max((@0-@1)*153.85, 0.0), 1.0)*(";

  formulaArgs.add(*obs_var);
  formulaArgs.add(*step);

  int xmax = obs_var->getMax();

  for (int i = 0; i < order; i += 2) {
    string pName = Form("%s_p%d", prefix.c_str(), i + 1, order);
    string cpName = Form("%s_cp%d", prefix.c_str(), i + 1, order);
    RooRealVar *p = new RooRealVar(pName.c_str(), pName.c_str(), -0.05*i, -0.5, 0.0);
    RooRealVar *cp = new RooRealVar(cpName.c_str(), cpName.c_str(), 0.9, 0., 1.);
    formulaArgs.add(*p);
    formulaArgs.add(*cp);
    formula += Form("@%d*TMath::Exp(@%d*@0)*@%d/(%d^(1+@%d)-@1^(1+@%d))", formulaArgs.getSize() - 1, formulaArgs.getSize() - 2, formulaArgs.getSize() - 1, xmax, formulaArgs.getSize() - 2, formulaArgs.getSize() - 2);
    if (i + 2 < order) {
      formula += "+";
    }
  }
  formula += ")";

  RooGenericPdf *stepPdf = new RooGenericPdf(Form("%s_step_exp%d", prefix.c_str(), order), Form("%s_step_exp%d", prefix.c_str(), order), formula.c_str(), formulaArgs);
  RooGaussModel *gau = new RooGaussModel(Form("%s_gau_exp%d", prefix.c_str(), order), Form("%s_gau_exp%d", prefix.c_str(), order), *obs_var, *mean, *sigma);
  RooFFTConvPdf *gauxexp = new RooFFTConvPdf(Form("%s_gauxexp%d", prefix.c_str(), order), Form("%s_gauxexp%d", prefix.c_str(), order), *obs_var, *stepPdf, *gau);
  return gauxexp;
}

RooAbsPdf *PdfModelBuilder::getLaurentStepxGau(string prefix, int order){
  if (order > 7) return NULL;
  RooRealVar *mean = new RooRealVar(Form("%s_mean", prefix.c_str()), Form("%s_mean", prefix.c_str()), 0.);
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()), 5., 1., 10.);
  RooRealVar *step = new RooRealVar(Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()), 101., 95., 120.);
  RooArgList formulaArgs;
  string formula = "TMath::Min(TMath::Max((@0-@1)*153.85, 0.0), 1.0)*(";

  formulaArgs.add(*obs_var);
  formulaArgs.add(*step);

  int xmax = obs_var->getMax();

  // ***************************************
  // Fixed the order to be added
  // ***************************************

  int nlower = int(ceil(order / 2.));
  int nhigher = order - nlower;

  RooRealVar *cp0 = new RooRealVar(Form("%s_cp0", prefix.c_str()), Form("%s_cp0_l%d", prefix.c_str(), order), 0.01, 0., 0.5);
  formulaArgs.add(*cp0);
  formula += Form("@%d*@0^(-3)*(-2)/(%d^(-2)-@1^(-2))+", formulaArgs.getSize() - 1, xmax);

  for (int i = 1; i <= nlower; i++) {
    RooRealVar *cp = new RooRealVar(Form("%s_cp%d", prefix.c_str(), i), Form("%s_cp%d", prefix.c_str(), i), 0.001, 0., 1.);
    formulaArgs.add(*cp);
    formula += Form("@%d*@0^(%d)*(%d)/(%d^(%d)-@1^(%d))", formulaArgs.getSize() - 1, -3 - i, -2 - i, xmax, -2 - i, -2 - i);
    if (i < nlower) {
      formula += "+";
    }
  }

  for (int i = 1; i <= nhigher; i++) {
    RooRealVar *cp = new RooRealVar(Form("%s_cp%d", prefix.c_str(), i + nlower), Form("%s_cp%d", prefix.c_str(), i + nlower), 0.01, 0., 0.5);
    formulaArgs.add(*cp);
    formula += Form("+@%d*@0^(%d)*(%d)/(%d^(%d)-@1^(%d))", formulaArgs.getSize() - 1, -3 + i, -2 + i, xmax, -2 + i, -2 + i);
    if (i < nhigher) {
      formula += "+";
    }
  }

  // // *************************************
  // // New implementation(NO fixed order)
  // // *************************************
  // for (int i = 0; i < order; i++) {
  //   RooRealVar *cp = new RooRealVar(Form("%s_cp%d", prefix.c_str(), i), Form("%s_cp%d", prefix.c_str(), i), 0.1*i, 0., 1.);
  //   RooRealVar *p = new RooRealVar(Form("%s_p%d", prefix.c_str(), i), Form("%s_p%d", prefix.c_str(), i), -2-i, -10, -0.9);
        
  //   formulaArgs.add(*p);
  //   formulaArgs.add(*cp);

  //   // Use p parameter only as a negative integer
  //   formula += Form("@%d*pow(@0,floor(@%d))*(1+floor(@%d))/(%d^(1+floor(@%d))-@1^(1+floor(@%d)))", formulaArgs.getSize() - 1, formulaArgs.getSize() - 2, formulaArgs.getSize() - 2, xmax, formulaArgs.getSize() - 2, formulaArgs.getSize() - 2);
  //   if (i + 1 < order) {
  //     formula += "+";
  //   }
  // }

  formula += ")";
  cout << "FORM -- " << formula << endl;
  RooGenericPdf *stepPdf = new RooGenericPdf(Form("%s_step_lau%d", prefix.c_str(), order), Form("%s_step_lau%d", prefix.c_str(), order), formula.c_str(), formulaArgs);
  RooGaussModel *gau = new RooGaussModel(Form("%s_gau_lau%d", prefix.c_str(), order), Form("%s_gau_lau%d", prefix.c_str(), order), *obs_var, *mean, *sigma);
  RooFFTConvPdf *gauxlau = new RooFFTConvPdf(Form("%s_gauxlau%d", prefix.c_str(), order), Form("%s_gauxlau%d", prefix.c_str(), order), *obs_var, *stepPdf, *gau);
  return gauxlau;
}

RooAbsPdf* PdfModelBuilder::getBernstein(string prefix, int order){
  
  RooArgList *coeffList = new RooArgList();
  //coeffList->add(RooConst(1.0)); // no need for cnstant in this interface
  for (int i=0; i<order; i++){
    string name = Form("%s_p%d",prefix.c_str(),i);
    //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.1*(i+1),-15.,15.);
    RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0*@0",RooArgList(*param));
    params.insert(pair<string,RooRealVar*>(name,param));
    prods.insert(pair<string,RooFormulaVar*>(name,form));
    coeffList->add(*prods[name]);
  }
  //RooBernstein *bern = new RooBernstein(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  if (order==1) {
	RooBernsteinFast<1> *bern = new RooBernsteinFast<1>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==2) {
	RooBernsteinFast<2> *bern = new RooBernsteinFast<2>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==3) {
	RooBernsteinFast<3> *bern = new RooBernsteinFast<3>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==4) {
	RooBernsteinFast<4> *bern = new RooBernsteinFast<4>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==5) {
	RooBernsteinFast<5> *bern = new RooBernsteinFast<5>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  // } else if (order==6) {
	// RooBernsteinFast<6> *bern = new RooBernsteinFast<6>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  // 	return bern;
//  } else if (order==7) {
//	RooBernsteinFast<7> *bern = new RooBernsteinFast<7>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
 // 	return bern;
  } else {
	return NULL;
  }
  //return bern;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(bern->GetName(),bern));

}

RooAbsPdf* PdfModelBuilder::getPowerLawGeneric(string prefix, int order){
  
  if (order%2==0){
    cerr << "ERROR -- addPowerLaw -- only odd number of params allowed" << endl;
    return NULL;
  }
  else {
    int nfracs=(order-1)/2;
    int npows=order-nfracs;
    assert(nfracs==npows-1);
    string formula="";
    RooArgList *dependents = new RooArgList();
    dependents->add(*obs_var);
    // first do recursive fraction
    if (order>1) {
      formula += "(1.-";
      for (int i=1; i<=nfracs; i++){
        if (i<nfracs) formula += Form("@%d-",i);
        else formula += Form("@%d)*",i);
        string name =  Form("%s_f%d",prefix.c_str(),i);
        params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.1,0.,1.)));
        dependents->add(*params[name]);
      }
    }
    for (int i=1; i<=npows; i++){
      string pname =  Form("%s_p%d",prefix.c_str(),i);
      string fname =  Form("%s_f%d",prefix.c_str(),i-1);
      params.insert(pair<string,RooRealVar*>(pname, new RooRealVar(pname.c_str(),pname.c_str(),TMath::Max(-10.,-2.*(i+1)),-10.,0.)));
      if (i==1) {
        formula += Form("TMath::Power(@0,@%d)",nfracs+i);
        dependents->add(*params[pname]);
      }
      else {
        formula += Form(" + @%d*TMath::Power(@0,@%d)",i-1,nfracs+i);
        dependents->add(*params[pname]);
      }
    }
    cout << "FORM -- " << formula << endl;
    dependents->Print("v");
    RooGenericPdf *pow = new RooGenericPdf(prefix.c_str(),prefix.c_str(),formula.c_str(),*dependents);
    pow->Print("v");
    return pow;
    //bkgPdfs.insert(pair<string,RooAbsPdf*>(pow->GetName(),pow));

  }
}

RooAbsPdf* PdfModelBuilder::getPowerLaw(string prefix, int order){
  
  RooArgList coefList;
  for (int i=0; i<order; i++){
    double start=-2.;
    double low=-10.;
    double high=0.;
    if (order>0){
      start=-0.001/double(i);
      low=-0.01;
      high=0.01;
    }
    RooRealVar *var = new RooRealVar(Form("%s_p%d",prefix.c_str(),i),Form("%s_p%d",prefix.c_str(),i),start,low,high);
    coefList.add(*var);
  }
  RooPowerLawSum *pow = new RooPowerLawSum(prefix.c_str(),prefix.c_str(),*obs_var,coefList);
  return pow;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(pow->GetName(),pow));

}

RooAbsPdf* PdfModelBuilder::getExponential(string prefix, int order){
  
  RooArgList coefList;
  for (int i=0; i<order; i++){
    double start=-1.;
    double low=-2.;
    double high=0.;
    if (order>0){
      start=-0.001/double(i);
      low=-0.01;
      high=0.01;
    }
    RooRealVar *var = new RooRealVar(Form("%s_p%d",prefix.c_str(),i),Form("%s_p%d",prefix.c_str(),i),start,low,high);
    coefList.add(*var);
  }
  RooPowerLawSum *exp = new RooPowerLawSum(prefix.c_str(),prefix.c_str(),*obs_var,coefList);
  return exp;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(exp->GetName(),exp));

}

RooAbsPdf* PdfModelBuilder::getPowerLawSingle(string prefix, int order){
  
  if (order%2==0){
    cerr << "ERROR -- addPowerLaw -- only odd number of params allowed" << endl;
    return NULL;
  }
  else {
    int nfracs=(order-1)/2;
    int npows=order-nfracs;
    assert(nfracs==npows-1);
    RooArgList *fracs = new RooArgList();
    RooArgList *pows = new RooArgList();
    for (int i=1; i<=nfracs; i++){
      string name =  Form("%s_f%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.9-float(i-1)*1./nfracs,0.,1.)));
      //params[name]->removeRange();
      fracs->add(*params[name]);
    }
    for (int i=1; i<=npows; i++){
      string name =  Form("%s_p%d",prefix.c_str(),i);
      string ename =  Form("%s_e%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),TMath::Max(-9.,-1.*(i+1)),-9.,1.)));
      //params[name]->removeRange();
      utilities.insert(pair<string,RooAbsPdf*>(ename, new RooPower(ename.c_str(),ename.c_str(),*obs_var,*params[name])));
      pows->add(*utilities[ename]);
    }
    //cout << "RooArgLists..." << endl;
    //fracs->Print("v");
    //pows->Print("v");
    //cout << "Function..." << endl;
    RooAbsPdf *pow = new RooAddPdf(prefix.c_str(),prefix.c_str(),*pows,*fracs,true); 
    //pow->Print("v");
    return pow;
    //bkgPdfs.insert(pair<string,RooAbsPdf*>(pow->GetName(),pow));
  }
}

RooAbsPdf* PdfModelBuilder::getLaurentSeries(string prefix, int order){
 
  int nlower=int(ceil(order/2.));
  int nhigher=order-nlower;
  // first do 0th order
  RooArgList *pows = new RooArgList();
  RooArgList *plist = new RooArgList();
  string pname =  Form("%s_pow0",prefix.c_str());
  utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(-4.))));
  pows->add(*utilities[pname]);

  // even terms
  for (int i=1; i<=nlower; i++){
    string name = Form("%s_l%d",prefix.c_str(),i);
    params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.25/order,0.000001,0.999999)));
    plist->add(*params[name]);
    string pname =  Form("%s_powl%d",prefix.c_str(),i);
    utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(-4.-i))));
    pows->add(*utilities[pname]);
  }
  // odd terms
  for (int i=1; i<=nhigher; i++){
    string name = Form("%s_h%d",prefix.c_str(),i);
    params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.25/order,0.000001,0.999999)));
    plist->add(*params[name]);
    string pname =  Form("%s_powh%d",prefix.c_str(),i);
    utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(-4.+i))));
    pows->add(*utilities[pname]);
  }
  RooAddPdf *pdf = new RooAddPdf(prefix.c_str(),prefix.c_str(),*pows,*plist,true);
  return pdf;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));
}

RooAbsPdf* PdfModelBuilder::getKeysPdf(string prefix){
  if (!keysPdfAttributesSet){
    cerr << "ERROR -- keysPdf attributes not set" << endl;
    exit(1);
  }
  return new RooKeysPdf(prefix.c_str(),prefix.c_str(),*obs_var,*keysPdfData,RooKeysPdf::MirrorBoth,keysPdfRho);
}

RooAbsPdf* PdfModelBuilder::getPdfFromFile(string &prefix){
  vector<string> details;
  split(details,prefix,boost::is_any_of(","));

  string fname = details[2];
  string wsname = details[1];
  string pdfname = details[0];

  TFile *tempFile = TFile::Open(fname.c_str());
  if (!tempFile){
    cerr << "PdfModelBuilder::getPdfFromFile -- file not found " << fname << endl;
    assert(0);
  }
  RooWorkspace *tempWS = (RooWorkspace*)tempFile->Get(wsname.c_str());
  if (!tempWS){
    cerr << "PdfModelBuilder::getPdfFromFile -- workspace not found " << wsname << endl;
    assert(0);
  }
  RooAbsPdf *tempPdf = (RooAbsPdf*)tempWS->pdf(pdfname.c_str());
  if (!tempPdf){
    cerr << "PdfModelBuilder::getPdfFromFile -- pdf not found " << pdfname << endl;
    assert(0);
  }
  prefix = pdfname;
  RooAbsPdf *pdf = (RooAbsPdf*)tempPdf->Clone(prefix.c_str());
  tempFile->Close();
  delete tempFile;
  return pdf;
}

RooAbsPdf* PdfModelBuilder::getExpModGaussian(string prefix){
  // Define parameters: mu, sigma and lambda
  RooRealVar *mu = new RooRealVar(Form("%s_mu",prefix.c_str()), Form("%s_mu",prefix.c_str()), 105., 95., 120.);
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma",prefix.c_str()), Form("%s_sigma",prefix.c_str()), 3., 1., 10.); // Adjusted lower limit
  RooRealVar *lambda = new RooRealVar(Form("%s_lambda",prefix.c_str()), Form("%s_lambda",prefix.c_str()), 0.1, 0.001, 0.5); // Keep upper limit limited for stability
  
  // Add parameters to the params map
  params.insert(pair<string,RooRealVar*>(Form("%s_mu",prefix.c_str()), mu));
  params.insert(pair<string,RooRealVar*>(Form("%s_sigma",prefix.c_str()), sigma));
  params.insert(pair<string,RooRealVar*>(Form("%s_lambda",prefix.c_str()), lambda));
  
  // Create list of dependent variables
  RooArgList formulaVars;
  formulaVars.add(*obs_var);    // x
  formulaVars.add(*mu);         // mu
  formulaVars.add(*sigma);      // sigma
  formulaVars.add(*lambda);     // lambda
  
  // Advanced formula with protection for both fitting and plotting
  // The formula has several safety features to prevent numerical issues:、、
  // 1. The exponent is constrained to prevent overflow
  string formula = "(@3/2) * exp(min(20.0, (@3/2) * (2*@1 + @3*@2*@2 - 2*@0))) * TMath::Erfc(((@1 + @3*@2*@2 - @0)/(sqrt(2)*@2)))";
  
  // Create RooGenericPdf with the stabilized formula
  RooGenericPdf *emg = new RooGenericPdf(prefix.c_str(), "Exponentially Modified Gaussian", formula.c_str(), formulaVars);
  
  return emg;
}

RooAbsPdf* PdfModelBuilder::getAsymGenGaussian(string prefix){
  // Define parameters: beta, alpha and mu
  // This function implements the formula: beta/(2*alpha*Gamma(1/beta)) * exp(-((|x-mu|/alpha)^beta))
  // where:
  // - beta is the shape parameter (controls the peakedness)
  // - alpha is the scale parameter (controls the width)
  // - mu is the location parameter (controls the center position)
  
  RooRealVar *beta = new RooRealVar(Form("%s_beta", prefix.c_str()), Form("%s_beta", prefix.c_str()), 2.0, 0.5, 5.0);
  RooRealVar *alpha = new RooRealVar(Form("%s_alpha", prefix.c_str()), Form("%s_alpha", prefix.c_str()), 1.0, 0.1, 10.0);
  RooRealVar *mu = new RooRealVar(Form("%s_mu", prefix.c_str()), Form("%s_mu", prefix.c_str()), 105.0, 95.0, 125.0);
  
  // Add parameters to the params map
  params.insert(pair<string,RooRealVar*>(Form("%s_beta", prefix.c_str()), beta));
  params.insert(pair<string,RooRealVar*>(Form("%s_alpha", prefix.c_str()), alpha));
  params.insert(pair<string,RooRealVar*>(Form("%s_mu", prefix.c_str()), mu));
  
  // Create a list of formula variables
  RooArgList formulaVars;
  formulaVars.add(*obs_var);    // x
  formulaVars.add(*mu);         // mu
  formulaVars.add(*alpha);      // alpha
  formulaVars.add(*beta);       // beta
  
  // Formula implementation of: beta/(2*alpha*Gamma(1/beta)) * exp(-((|x-mu|/alpha)^beta))
  // We use TMath::Gamma(1/beta) for the gamma function
  // We need to ensure numerical stability:
  // 1. Prevent division by zero in gamma function
  // 2. Prevent overflow in the exponent
  // 3. Ensure the absolute value is properly handled
  string formula = "@3/(2*@2*TMath::Gamma(1.0/@3)) * exp(-pow(abs(@0-@1)/@2, @3))";
  
  // Create the generalized Gaussian PDF with the specified formula
  RooGenericPdf *agg = new RooGenericPdf(prefix.c_str(), "Generalized Gaussian Distribution", formula.c_str(), formulaVars);
  
  return agg;
}

RooAbsPdf* PdfModelBuilder::getExponentialSingle(string prefix, int order){
  
  if (order%2==0){
    cerr << "ERROR -- addExponential -- only odd number of params allowed" << endl;
    return NULL;
  }
  else {
    int nfracs=(order-1)/2;
    int nexps=order-nfracs;
    assert(nfracs==nexps-1);
    RooArgList *fracs = new RooArgList();
    RooArgList *exps = new RooArgList();
    for (int i=1; i<=nfracs; i++){
      string name =  Form("%s_f%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.9-float(i-1)*1./nfracs,0.0001,0.9999)));
      fracs->add(*params[name]);
    }
    for (int i=1; i<=nexps; i++){
      string name =  Form("%s_p%d",prefix.c_str(),i);
      string ename =  Form("%s_e%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),TMath::Max(-1.,-0.04*(i+1)),-1.,0.)));
      utilities.insert(pair<string,RooAbsPdf*>(ename, new RooExponential(ename.c_str(),ename.c_str(),*obs_var,*params[name])));
      exps->add(*utilities[ename]);
    }
    //fracs->Print("v");
    //exps->Print("v");
    RooAbsPdf *exp = new RooAddPdf(prefix.c_str(),prefix.c_str(),*exps,*fracs,true);
    //exp->Print("v");
    cout << "--------------------------" << endl;
    return exp;
    //bkgPdfs.insert(pair<string,RooAbsPdf*>(exp->GetName(),exp));

  }
}


void PdfModelBuilder::addBkgPdf(string type, int nParams, string name, bool cache){
 
  if (!obs_var_set){
    cerr << "ERROR -- obs Var has not been set!" << endl;
    exit(1);
  }
  bool found=false;
  for (vector<string>::iterator it=recognisedPdfTypes.begin(); it!=recognisedPdfTypes.end(); it++){
    if (*it==type) found=true;
  }
  if (!found){
    cerr << "Pdf of type " << type << " is not recognised!" << endl;
    exit(1);
  }
  RooAbsPdf *pdf=0;// avoid uninitialised variable error in cmssw  

  if (type=="Bernstein") pdf = getBernstein(name,nParams);
  if (type=="Exponential") pdf = getExponentialSingle(name,nParams);
  if (type=="PowerLaw") pdf = getPowerLawSingle(name,nParams);
  if (type=="Laurent") pdf = getLaurentSeries(name,nParams);
  if (type=="KeysPdf") pdf = getKeysPdf(name);
  if (type=="File") pdf = getPdfFromFile(name);
  if (type=="ExpModGaussian") pdf = getExpModGaussian(name);
  if (type=="AsymGenGaussian") pdf = getAsymGenGaussian(name);

  if (cache) {
    wsCache->import(*pdf);
    RooAbsPdf *cachePdf = wsCache->pdf(pdf->GetName());
    bkgPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
  }
  else {
    bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));
  }

}

void PdfModelBuilder::setKeysPdfAttributes(RooDataSet *data, double rho){
  keysPdfData = data;
  keysPdfRho = rho;
  keysPdfAttributesSet=true;
}

void PdfModelBuilder::setSignalPdf(RooAbsPdf *pdf, RooRealVar *norm){
  sigPdf=pdf;
  sigNorm=norm;
  signal_set=true;
}

void PdfModelBuilder::setSignalPdfFromMC(RooDataSet *data){
  
  RooDataHist *sigMCBinned = new RooDataHist(Form("roohist_%s",data->GetName()),Form("roohist_%s",data->GetName()),RooArgSet(*obs_var),*data);
  sigPdf = new RooHistPdf(Form("pdf_%s",data->GetName()),Form("pdf_%s",data->GetName()),RooArgSet(*obs_var),*sigMCBinned);
  sigNorm = new RooConstVar(Form("sig_events_%s",data->GetName()),Form("sig_events_%s",data->GetName()),data->sumEntries());
  signal_set=true;
}

void PdfModelBuilder::makeSBPdfs(bool cache){
  
  if (!signal_set){
    cerr << "ERROR - no signal model set!" << endl;
    exit(1);
  }
  if (!signal_modifier_set){
    cerr << "ERROR - no signal modifier set!" << endl;
    exit(1);
  }
 
  if (sigNorm) {
    sigYield = new RooProduct("sig_yield","sig_yield",RooArgSet(*signalModifier,*sigNorm));
  }
  else {
    sigYield = signalModifier;
  }
  bkgYield = new RooRealVar("bkg_yield","bkg_yield",1000.,0.,1.e6);

  for (map<string,RooAbsPdf*>::iterator bkg=bkgPdfs.begin(); bkg!=bkgPdfs.end(); bkg++){
    RooAbsPdf *sbMod = new RooAddPdf(Form("sb_%s",bkg->first.c_str()),Form("sb_%s",bkg->first.c_str()),RooArgList(*(bkg->second),*sigPdf),RooArgList(*bkgYield,*sigYield));
    if (cache) {
      wsCache->import(*sbMod,RecycleConflictNodes());
      RooAbsPdf *cachePdf = (RooAbsPdf*)wsCache->pdf(sbMod->GetName());
      signalModifier = (RooRealVar*)wsCache->var(signalModifier->GetName());
      sbPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
    }
    else {
      sbPdfs.insert(pair<string,RooAbsPdf*>(sbMod->GetName(),sbMod));
    }
  }
}

map<string,RooAbsPdf*> PdfModelBuilder::getBkgPdfs(){
  return bkgPdfs;
}

map<string,RooAbsPdf*> PdfModelBuilder::getSBPdfs(){
  return sbPdfs;
}

RooAbsPdf* PdfModelBuilder::getSigPdf(){
  return sigPdf;
}

void PdfModelBuilder::plotPdfsToData(RooAbsData *data, int binning, string name, bool bkgOnly,string specificPdfName){
  
  TCanvas *canv = new TCanvas();
  bool specPdf=false;
  if (specificPdfName!="") specPdf=true;

  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) pdfSet = bkgPdfs;
  else pdfSet = sbPdfs;
  
  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    if (specPdf && it->first!=specificPdfName && specificPdfName!="NONE") continue;
    RooPlot *plot = obs_var->frame();
    data->plotOn(plot,Binning(binning));
    if (specificPdfName!="NONE") {
	 it->second->plotOn(plot);
	 it->second->paramOn(plot,RooFit::Layout(0.34,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
    }	
    plot->Draw();
    canv->Print(Form("%s_%s.pdf",name.c_str(),it->first.c_str()));
    canv->Print(Form("%s_%s.png",name.c_str(),it->first.c_str()));
  }
  delete canv;
}

void PdfModelBuilder::fitToData(RooAbsData *data, bool bkgOnly, bool cache, bool print){
  
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) pdfSet = bkgPdfs;
  else pdfSet = sbPdfs;

  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    RooFitResult *fit = (RooFitResult*)it->second->fitTo(*data,Save(true));
    if (print){
      cout << "Fit Res Before: " << endl;
      fit->floatParsInit().Print("v");
      cout << "Fit Res After: " << endl;
      fit->floatParsFinal().Print("v");
    }
    if (cache) {
      RooArgSet *fitargs = (RooArgSet*)it->second->getParameters(*obs_var);
      // remove the signal strength since this will be set AFTER fitting the background 
      fitargs->remove(*signalModifier); 
      wsCache->defineSet(Form("%s_params",it->first.c_str()),*fitargs);
      wsCache->defineSet(Form("%s_observs",it->first.c_str()),*obs_var);
      wsCache->saveSnapshot(it->first.c_str(),*fitargs,true);
      if (print) {
        cout << "Cached values: " << endl;
        fitargs->Print("v");
      }
    }
  }
  if (bkgOnly) bkgHasFit=true;
  else sbHasFit=true;
}

void PdfModelBuilder::setSeed(int seed){
  RooRandom::randomGenerator()->SetSeed(seed);
}

RooDataSet* PdfModelBuilder::makeHybridDataset(vector<float> switchOverMasses, vector<RooDataSet*> dataForHybrid){
  
  assert(switchOverMasses.size()==dataForHybrid.size()-1);

  vector<string> cut_strings;
  cut_strings.push_back("cutstring0");
  obs_var->setRange("cutstring0",obs_var->getMin(),switchOverMasses[0]);
  for (unsigned int i=1; i<switchOverMasses.size(); i++){
    cut_strings.push_back(Form("cutstring%d",i));
    obs_var->setRange(Form("cutstring%d",i),switchOverMasses[i-1],switchOverMasses[i]);
  }
  cut_strings.push_back(Form("cutstring%d",int(switchOverMasses.size())));
  obs_var->setRange(Form("cutstring%d",int(switchOverMasses.size())),switchOverMasses[switchOverMasses.size()-1],obs_var->getMax());
  
  obs_var->Print("v");
  assert(cut_strings.size()==dataForHybrid.size());
  
	RooDataSet *data=0;// avoid uninitialised variable error in cmssw  
  for (unsigned int i=0; i<dataForHybrid.size(); i++){
    RooDataSet *cutData = (RooDataSet*)dataForHybrid[i]->reduce(Name("hybridToy"),Title("hybridToy"),CutRange(cut_strings[i].c_str()));
    //RooDataSet *cutData = new RooDataSet("hybridToy","hybridToy",RooArgSet(*obs_var),Import(*dataForHybrid[i]),CutRange(cut_strings[i].c_str()));
    if (i==0) data=cutData;
    else data->append(*cutData);
  }
  return data;
}

void PdfModelBuilder::throwHybridToy(string postfix, int nEvents, vector<float> switchOverMasses, vector<string> functions, bool bkgOnly, bool binned, bool poisson, bool cache){
  
  assert(switchOverMasses.size()==functions.size()-1);
  toyHybridData.clear();

  // have to throw unbinned for the hybrid
  throwToy(postfix,nEvents,bkgOnly,false,poisson,cache);

  vector<RooDataSet*> dataForHybrid;
  string hybridName = "hybrid";
  for (vector<string>::iterator func=functions.begin(); func!=functions.end(); func++){
    hybridName += "_"+*func;
    for (map<string,RooDataSet*>::iterator it=toyDataSet.begin(); it!=toyDataSet.end(); it++){
      if (it->first.find(*func)!=string::npos){
        dataForHybrid.push_back(it->second);
      }
    }
  }
  if (dataForHybrid.size()!=functions.size()){
    cerr << "One of the requested hybrid functions has not been found" << endl;
    exit(1);
  }

  RooDataSet *hybridData = makeHybridDataset(switchOverMasses,dataForHybrid);
  toyHybridData.clear();
  if (binned) {
    RooDataHist *hybridDataHist = hybridData->binnedClone();
    hybridDataHist->SetName(Form("%s_%s",hybridName.c_str(),postfix.c_str()));
    toyHybridData.insert(pair<string,RooAbsData*>(hybridDataHist->GetName(),hybridDataHist));
  }
  else {
    hybridData->SetName(Form("%s_%s",hybridName.c_str(),postfix.c_str()));
    toyHybridData.insert(pair<string,RooAbsData*>(hybridData->GetName(),hybridData));
  }
}

void PdfModelBuilder::throwToy(string postfix, int nEvents, bool bkgOnly, bool binned, bool poisson, bool cache){

  toyData.clear();
  toyDataSet.clear();
  toyDataHist.clear();
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) {
    pdfSet = bkgPdfs;
    if (!bkgHasFit) cerr << "WARNING -- bkg has not been fit to data. Are you sure this is wise?" << endl; 
  }
  else {
    pdfSet = sbPdfs;
    if (!sbHasFit) cerr << "WARNING -- sb has not been fit to data. Are you sure this is wise?" << endl;
  }
  
  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    if (cache) {
      wsCache->loadSnapshot(it->first.c_str());
      cout << "Loaded snapshot, params at.." << endl;
      it->second->getVariables()->Print("v");
    }
    RooAbsData *toy;
    if (binned){
      RooDataHist *toyHist;
      if (poisson) toyHist = it->second->generateBinned(RooArgSet(*obs_var),nEvents,Extended(),Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      else toyHist = it->second->generateBinned(RooArgSet(*obs_var),nEvents,Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      toyDataHist.insert(pair<string,RooDataHist*>(toyHist->GetName(),toyHist));
      toy=toyHist;
    }
    else {
      RooDataSet *toySet;
      if (poisson) toySet = it->second->generate(RooArgSet(*obs_var),nEvents,Extended(),Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      else toySet = it->second->generate(RooArgSet(*obs_var),nEvents,Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      toyDataSet.insert(pair<string,RooDataSet*>(toySet->GetName(),toySet));
      toy=toySet;
    }
    toyData.insert(pair<string,RooAbsData*>(toy->GetName(),toy));
  }
  
}

map<string,RooAbsData*> PdfModelBuilder::getToyData(){
  return toyData;
}

map<string,RooAbsData*> PdfModelBuilder::getHybridToyData(){
  return toyHybridData;
}

void PdfModelBuilder::plotHybridToy(string prefix, int binning, vector<float> switchOverMasses, vector<string> functions, bool bkgOnly){

  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) {
    pdfSet = bkgPdfs;
  }
  else {
    pdfSet = sbPdfs;
  }
  
  int tempColors[4] = {kBlue,kRed,kGreen+2,kMagenta};

  vector<string> cut_strings;
  cut_strings.push_back("cutstring0");
  obs_var->setRange("cutstring0",obs_var->getMin(),switchOverMasses[0]);
  for (unsigned int i=1; i<switchOverMasses.size(); i++){
    cut_strings.push_back(Form("cutstring%d",i));
    obs_var->setRange(Form("cutstring%d",i),switchOverMasses[i-1],switchOverMasses[i]);
  }
  cut_strings.push_back(Form("cutstring%d",int(switchOverMasses.size())));
  obs_var->setRange(Form("cutstring%d",int(switchOverMasses.size())),switchOverMasses[switchOverMasses.size()-1],obs_var->getMax());
  
  RooPlot *plot = obs_var->frame();
  TCanvas *canv = new TCanvas();
  int i=0;
  for (vector<string>::iterator func=functions.begin(); func!=functions.end(); func++){
    for (map<string,RooAbsPdf*>::iterator pdfIt = pdfSet.begin(); pdfIt != pdfSet.end(); pdfIt++){
      // check if in list of hybrid functions
      if (pdfIt->first.find(*func)!=string::npos) {
        for (map<string,RooAbsData*>::iterator toyIt = toyData.begin(); toyIt != toyData.end(); toyIt++){
          //cout << "pdf: " << pdfIt->first << " - toy: " << toyIt->first << endl; 
          if (toyIt->first.find(pdfIt->first)!=string::npos){
            RooAbsData *data = toyIt->second->reduce(CutRange(cut_strings[i].c_str()));
            data->plotOn(plot,Binning(binning),MarkerColor(tempColors[i]),LineColor(tempColors[i]),CutRange(cut_strings[i].c_str()));
            pdfIt->second->plotOn(plot,LineColor(tempColors[i]),Range(cut_strings[i].c_str()));
            i++;
          }
        }
      }
    }
  }
  for (map<string,RooAbsData*>::iterator hybrid=toyHybridData.begin(); hybrid!=toyHybridData.end(); hybrid++){
    hybrid->second->plotOn(plot,Binning(binning),MarkerSize(0.8),MarkerStyle(kFullSquare));
    plot->SetMinimum(0.0001);
    plot->Draw();
    canv->Print(Form("%s_%s.pdf",prefix.c_str(),hybrid->first.c_str()));
  }
  delete canv;
}

void PdfModelBuilder::plotToysWithPdfs(string prefix, int binning, bool bkgOnly){
  
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) {
    pdfSet = bkgPdfs;
  }
  else {
    pdfSet = sbPdfs;
  }
  TCanvas *canv = new TCanvas();
  for (map<string,RooAbsPdf*>::iterator pdfIt = pdfSet.begin(); pdfIt != pdfSet.end(); pdfIt++){
    for (map<string,RooAbsData*>::iterator toyIt = toyData.begin(); toyIt != toyData.end(); toyIt++){
      //cout << "pdf: " << pdfIt->first << " - toy: " << toyIt->first << endl; 
      if (toyIt->first.find(pdfIt->first)!=string::npos){
        RooPlot *plot = obs_var->frame();
        toyIt->second->plotOn(plot,Binning(binning));
        pdfIt->second->plotOn(plot,LineColor(kRed));
        pdfIt->second->paramOn(plot,LineColor(kRed),RooFit::Layout(0.34,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
        plot->Draw();
        canv->Print(Form("%s_%s.pdf",prefix.c_str(),pdfIt->first.c_str()));
        canv->Print(Form("%s_%s.png",prefix.c_str(),pdfIt->first.c_str()));
      }
    }
  }
  delete canv;

}

void PdfModelBuilder::saveWorkspace(string filename){
  
  TFile *outFile = new TFile(filename.c_str(),"RECREATE");
  outFile->cd();
  wsCache->Write();
  outFile->Close();
  delete outFile;
}

void PdfModelBuilder::saveWorkspace(TFile *file){
  file->cd();
  wsCache->Write();
}
