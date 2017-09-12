#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <sstream>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TMultiGraph.h"

// lhapdf routines
#include "LHAPDF/LHAPDF.h"

using namespace std;
using namespace LHAPDF;


int main(int argc, char **argv)
{  
  cout << "\n ********************************************\n";
  cout <<   " * PDF reweighting with Lattice QCD moments *\n";
  cout <<   " *******************************************\n\n";
  cout << endl;

  //************************************************
  // Reading the moments
  //************************************************

  // Read the moments computed with Emanuele
  // Loop over the unpolarized and polarized cases in turn
  int const NrepUnpol=1000;
  int const NrepPol=100;
  int nrep=0;

  for(int ipol=1;ipol<2;ipol++){

    // Set number of replicas
    if(ipol==0) nrep=NrepUnpol;
    if(ipol==1) nrep=NrepPol;
  
    ifstream in1;
    // IPOL=0 => This is \int_0^1 x * ( d(x,Q) + dbar(x,Q) ), Eq (A.4) of whitepaper
    // IPOL=1 => This is \int_0^1 * ( Delta_d(x,Q) + Delta_dbar(x,Q) ), Eq (A.9) of whitepaper
    double unp_mom_rep_xdp[NrepUnpol+1]={0.0};
    if(ipol==0)in1.open("../code_moments/replicas/res/unp_mom_rep_xdp.res");
    if(ipol==1)in1.open("../code_moments/replicas/res/pol_mom_dp.res");
    for(int irep=1;irep<=nrep;irep++){
      int idum=0;
      in1>>idum>> unp_mom_rep_xdp[irep];
      cout<<irep<<" "<<unp_mom_rep_xdp[irep]<<endl;
      if(idum!=irep){
	cout<<"Error when reading replica file"<<endl;
	exit(-10);
      }
    }
    in1.close();

    // IPOL=0 => This is \int_0^1 x * ( u(x,Q) + ubar(x,Q) ), Eq (A.4) of whitepaper
    // IPOL=1 => This is \int_0^1 * ( Delta_u(x,Q) + Delta_ubar(x,Q) ), Eq (A.9) of whitepaper
    double unp_mom_rep_xup[NrepUnpol+1]={0.0};
    if(ipol==0)in1.open("../code_moments/replicas/res/unp_mom_rep_xup.res");
    if(ipol==1)in1.open("../code_moments/replicas/res/pol_mom_up.res");
    for(int irep=1;irep<=nrep;irep++){
      int idum=0;
      in1>>idum>> unp_mom_rep_xup[irep];
      cout<<irep<<" "<<unp_mom_rep_xup[irep]<<endl;
      if(idum!=irep){
	cout<<"Error when reading replica file"<<endl;
	exit(-10);
      }
    }
    in1.close();
    
    // IPOL=0 => This is \int_0^1 x * ( s(x,Q) + sbar(x,Q) ), Eq (A.4) of whitepaper
    // IPOL=1 => This is \int_0^1 * ( Delta_s(x,Q) + Delta_sbar(x,Q) ), Eq (A.9) of whitepaper
    double unp_mom_rep_xsp[NrepUnpol+1]={0.0};
    if(ipol==0)in1.open("../code_moments/replicas/res/unp_mom_rep_xsp.res");
    if(ipol==1)in1.open("../code_moments/replicas/res/pol_mom_sp.res");
    for(int irep=1;irep<=nrep;irep++){
      int idum=0;
      in1>>idum>> unp_mom_rep_xsp[irep];
      if(idum!=irep){
	cout<<"Error when reading replica file"<<endl;
	exit(-10);
      }
    }
    in1.close();

    //  IPOL=0 => This is \int_0^1 x * ( g(x,Q) ), Eq (A.4) of whitepaper
    //  IPOL=1 => This is \int_0^1 * x ( Delta_u(x,Q) - Delta_ubar(x,Q)
    // - Delta_d(x,Q) + Delta_dbar(x,Q) ), Eq (A.8) of whitepaper
    double unp_mom_rep_xg[NrepUnpol+1]={0.0};
    if(ipol==0)in1.open("../code_moments/replicas/res/unp_mom_rep_xg.res");
    if(ipol==1)in1.open("../code_moments/replicas/res/pol_mom_xummdm.res");
    for(int irep=1;irep<=nrep;irep++){
      int idum=0;
      in1>>idum>> unp_mom_rep_xg[irep];
      if(idum!=irep){
	cout<<"Error when reading replica file"<<endl;
	exit(-10);
      }
    }
    in1.close();

    // IPOL=0 => This is \int_0^1 x * ( u(x,Q) + ubar(x,Q) - d(x,Q) - dbar(x,Q) ), Eq (A.2) of whitepaper
    // IPOL=1 => This is \int_0^1 ( Delta_u(x,Q) + Delta_ubar(x,Q)
    // - Delta_d(x,Q) - Delta_dbar(x,Q) ), Eq (A.7) of whitepaper
    double unp_mom_rep_xupmdp[NrepUnpol+1]={0.0};
    if(ipol==0)in1.open("../code_moments/replicas/res/unp_mom_rep_xupmdp.res");
    if(ipol==1)in1.open("../code_moments/replicas/res/pol_mom_upmdp.res");
    for(int irep=1;irep<=nrep;irep++){
      int idum=0;
      in1>>idum>> unp_mom_rep_xupmdp[irep];
      if(idum!=irep){
	cout<<"Error when reading replica file"<<endl;
	exit(-10);
      }
    }
    in1.close();
    
  //************************************************
  // Compute means and variances
  //************************************************
  
  // Compute mean and variance for XDP
  double unp_mom_xdp_mean=0;
  double unp_mom_xdp_std=0;

  for(int irep=1;irep<=nrep;irep++){
    unp_mom_xdp_mean+=unp_mom_rep_xdp[irep]/nrep;
    unp_mom_xdp_std+=pow(unp_mom_rep_xdp[irep],2.0)/nrep;
  }
  unp_mom_xdp_std=pow(unp_mom_xdp_std - pow(unp_mom_xdp_mean,2.0),0.5);

  // Compute mean and variance for XUP
  double unp_mom_xup_mean=0;
  double unp_mom_xup_std=0;

  for(int irep=1;irep<=nrep;irep++){
    unp_mom_xup_mean+=unp_mom_rep_xup[irep]/nrep;
    unp_mom_xup_std+=pow(unp_mom_rep_xup[irep],2.0)/nrep;
  }
  unp_mom_xup_std=pow(unp_mom_xup_std - pow(unp_mom_xup_mean,2.0),0.5);

  // Compute mean and variance for XSP
  double unp_mom_xsp_mean=0;
  double unp_mom_xsp_std=0;

  for(int irep=1;irep<=nrep;irep++){
    unp_mom_xsp_mean+=unp_mom_rep_xsp[irep]/nrep;
    unp_mom_xsp_std+=pow(unp_mom_rep_xsp[irep],2.0)/nrep;
  }
  unp_mom_xsp_std=pow(unp_mom_xsp_std - pow(unp_mom_xsp_mean,2.0),0.5);

  // Compute mean and variance for XG
  double unp_mom_xg_mean=0;
  double unp_mom_xg_std=0;

  for(int irep=1;irep<=nrep;irep++){
    unp_mom_xg_mean+=unp_mom_rep_xg[irep]/nrep;
    unp_mom_xg_std+=pow(unp_mom_rep_xg[irep],2.0)/nrep;
  }
  unp_mom_xg_std=pow(unp_mom_xg_std - pow(unp_mom_xg_mean,2.0),0.5);

  // Compute mean and variance for XUPMDP
  double unp_mom_xupmdp_mean=0;
  double unp_mom_xupmdp_std=0;

  for(int irep=1;irep<=nrep;irep++){
    unp_mom_xupmdp_mean+=unp_mom_rep_xupmdp[irep]/nrep;
    unp_mom_xupmdp_std+=pow(unp_mom_rep_xupmdp[irep],2.0)/nrep;
  }
  unp_mom_xupmdp_std=pow(unp_mom_xupmdp_std - pow(unp_mom_xupmdp_mean,2.0),0.5);


  //*************************************
  // Prepare PDF comparison plots
  //*************************************
  
  // Now do some plots of the PDFs before and after reweighting
  double const Q = pow(4,0.5); // Scale where PDFs are compared - same as lattice results
  
  std::vector<Color_t> color;
  color.push_back(kGreen);
  color.push_back(kRed);
  color.push_back(kBlue);

  std::vector<Color_t> color2;
  color2.push_back(kGreen+3);
  color2.push_back(kRed+1);
  color2.push_back(kBlue+2);
  
  std::vector<std::string> pdfname;
  pdfname.push_back("xg");
  pdfname.push_back("xup");
  pdfname.push_back("xdp");
  pdfname.push_back("xsp");
  pdfname.push_back("xupmdp");

  std::vector<std::string> titleY;
  titleY.push_back("g ( x, Q^{2} = 4 GeV^{2} ) ");
  titleY.push_back("u^{+} ( x, Q^{2} = 4 GeV^{2} ) ");
  titleY.push_back("d^{+} ( x, Q^{2} = 4 GeV^{2} ) ");
  titleY.push_back("s^{+} ( x, Q^{2} = 4 GeV^{2} ) ");
  titleY.push_back("( u^{+}-d^{+} ) ( x, Q^{2} = 4 GeV^{2} ) ");
 
  double const xmin=1e-3;
  double const xmax=0.90;

  int const nx=200;

  std::vector<std::string> outfilename;
  std::vector<std::string> outfilename2;
  for(unsigned i=0;i< pdfname.size();i++){
    outfilename.push_back(pdfname.at(i)+"-lattice");
    outfilename2.push_back(pdfname.at(i)+"-lattice-relerr");
  }
  
  // Define the PDF prior
  string pdfset;
  if(ipol==0) pdfset="NNPDF31_nnlo_as_0118_1000";
  if(ipol==1) pdfset="NNPDFpol11_100";
  // Initialization
  LHAPDF::initPDFSet(pdfset);
  // Set low verbosity as usual
  LHAPDF::setVerbosity(0);

  // Loop over PDf flavors
  for (size_t l = 0; l < pdfname.size(); l++)
    {
           
      //*****************************
      // Prepare the plots
      //*****************************
      std::vector<TCanvas*> c;
      c.push_back(new TCanvas());
      c[0]->SetTickx();
      c[0]->SetTicky();
      c[0]->SetLogx();
      
      std::vector<TLegend*> leg;
      leg.push_back(new TLegend(0.11,0.65,0.49,0.89));
      leg[0]->SetLineStyle(1);
      leg[0]->SetBorderSize(1);
      leg[0]->SetFillColor(kWhite);
      
      //************************************************
      // Compute the values of the chi2
      //************************************************
      
      // Conservative: 5% total systematics
      // Optimistic: 1% total systematics
      int const n_scenarios=3;
      double total_latQCD_sysunc[n_scenarios]={0.0};
      if(ipol==0){
	total_latQCD_sysunc[0] = 0.01;
	total_latQCD_sysunc[1] = 0.03;
	total_latQCD_sysunc[2] = 0.05;
      }
      if(ipol==1){
	total_latQCD_sysunc[0] = 0.10;
	total_latQCD_sysunc[1] = 0.20;
	total_latQCD_sysunc[2] = 0.30;
      }
      
      
      // Loop over the number of scenarios
      for(int i_scen=0;i_scen <n_scenarios;i_scen++ ){
	
	// Compute the chi2 for each replica
	double chi2[NrepUnpol+1] = {0.0};
	
	// Fix here the number of moments that are included in the chi2 computation
	int const ndatfit=5;
	for(int irep=1;irep<=nrep;irep++){
	  
	  // xdp
	  chi2[irep] += pow( unp_mom_xdp_mean - unp_mom_rep_xdp[irep] ,2.0)/
	    pow(total_latQCD_sysunc[i_scen]*unp_mom_xdp_mean ,2.0);
	  
	  // xup
	  chi2[irep] += pow( unp_mom_xup_mean - unp_mom_rep_xup[irep] ,2.0)/
	    pow(total_latQCD_sysunc[i_scen]*unp_mom_xup_mean ,2.0);
	  
	  // xsp
	  chi2[irep] += pow( unp_mom_xsp_mean - unp_mom_rep_xsp[irep] ,2.0)/
	    pow(total_latQCD_sysunc[i_scen]*unp_mom_xsp_mean ,2.0);
	  
	  // xg
	  chi2[irep] += pow( unp_mom_xg_mean - unp_mom_rep_xg[irep] ,2.0)/
	    pow(total_latQCD_sysunc[i_scen]*unp_mom_xg_mean ,2.0);
	  
	  // xupmdp
	  chi2[irep] += pow( unp_mom_xupmdp_mean - unp_mom_rep_xupmdp[irep] ,2.0)/
	    pow(total_latQCD_sysunc[i_scen]*unp_mom_xupmdp_mean ,2.0);
	  
	  chi2[irep]/=ndatfit;
	  //cout<<irep<<" "<<chi2[irep]<<endl;
	}
	
	// Compute the weights
	double weights[NrepUnpol+1]={0.0};
	double sum=0.0;
	int const irep_ref=1;
	for(int irep=1;irep<=nrep;irep++){
	  double logw = ( double(ndatfit-1)/2.0 ) * log( chi2[irep] * ndatfit )
	    - ( chi2[irep] * ndatfit ) / 2.0;
	  //cout<<"logw"<<" "<<logw<<endl;
	  logw -= ( (ndatfit-1)/2 ) * log( chi2[irep_ref] * ndatfit )
	    - ( chi2[irep_ref] * ndatfit ) / 2.0;
	  weights[irep]=exp(logw);
	  //cout<<irep<<" "<<logw<<" "<<weights[irep]<<endl;
	  
	  sum += weights[irep] / nrep;
	}
	
	// Check normalization
	double check=0;
	for(int irep=1;irep<=nrep;irep++){
	  weights[irep]/=sum;
	  check+=weights[irep];
	  // cout<<"irep, weights =  "<<irep<<"  "<<weights[irep]<<endl;
	}
	std::cout<<"\n Sum of weights = "<<check<<"\n"<<std::endl;
	if(fabs(check-nrep)>1e-3){
	  cout<<"Problem with the computation of the weights"<<endl;
	  exit(-10);
	}
	
	// Now compute the effective number of replicas
	double shannon=0.0;
	for(int irep=1;irep<=nrep;irep++){
	  if(weights[irep]>1e-8){
	    shannon += weights[irep] * log(double(nrep)/weights[irep]) / nrep;
	  }
	}
	double const neff=exp(shannon);
	std::cout<<"Effective number of replicas = "<<neff<<"\n "<<std::endl;
	
	// ******************************************************
	// Now compute the reweighted moments
	// ******************************************************
	
	// -- xdp -------------------------------------------------------
	
	if(ipol==0){
	  cout<<"\nOriginal Moment: \\int_0^1 * x * ( d(x,Q) + dbar(x,Q) )"<<endl;
	}
	if(ipol==1){
	  cout<<"\nOriginal Moment: \\int_0^1 * ( Delta_d(x,Q) + Delta_dbar(x,Q) )"<<endl;
	}
	cout<<"Mean +- Std = "<<unp_mom_xdp_mean<<"  "<<unp_mom_xdp_std<<
	  "  ( "<<100*unp_mom_xdp_std/unp_mom_xdp_mean<<" % )"<<endl;
	
	cout<<"Reweighted Moment: "<<endl;
	
	// Compute mean and variance
	double unp_mom_xdp_mean_rw=0;
	double unp_mom_xdp_std_rw=0;
	
	for(int irep=1;irep<=nrep;irep++){
	  unp_mom_xdp_mean_rw+=weights[irep]*unp_mom_rep_xdp[irep]/nrep;
	  unp_mom_xdp_std_rw+=weights[irep]*pow(unp_mom_rep_xdp[irep],2.0)/nrep;
	}
	unp_mom_xdp_std_rw=pow(unp_mom_xdp_std_rw - pow(unp_mom_xdp_mean_rw,2.0),0.5);
	
	cout<<"Mean +- Std = "<<unp_mom_xdp_mean_rw<<"  "<<unp_mom_xdp_std_rw<<
	  "  ( "<<100*unp_mom_xdp_std_rw/unp_mom_xdp_mean_rw<<" % )"<<endl;

	// -- xup -------------------------------------------------------
	
	if(ipol==0){
	  cout<<"\nOriginal Moment: \\int_0^1 x * ( u(x,Q) + ubar(x,Q) )"<<endl;
	}
	if(ipol==1){
	  cout<<"\nOriginal Moment: \\int_0^1 * ( Delta_u(x,Q) + Delta_ubar(x,Q) )"<<endl;
	}
	cout<<"Mean +- Std = "<<unp_mom_xup_mean<<"  "<<unp_mom_xup_std<<
	  "  ( "<<100*unp_mom_xup_std/unp_mom_xup_mean<<" % )"<<endl;
	
	cout<<"\n Reweighted Moment: "<<endl;
	
	// Compute mean and variance
	double unp_mom_xup_mean_rw=0;
	double unp_mom_xup_std_rw=0;
	
	for(int irep=1;irep<=nrep;irep++){
	  unp_mom_xup_mean_rw+=weights[irep]*unp_mom_rep_xup[irep]/nrep;
	  unp_mom_xup_std_rw+=weights[irep]*pow(unp_mom_rep_xup[irep],2.0)/nrep;
	}
	unp_mom_xup_std_rw=pow(unp_mom_xup_std_rw - pow(unp_mom_xup_mean_rw,2.0),0.5);
	
	cout<<"Mean +- Std = "<<unp_mom_xup_mean_rw<<"  "<<unp_mom_xup_std_rw<<
	  "  ( "<<100*unp_mom_xup_std_rw/unp_mom_xup_mean_rw<<" % )"<<endl;
	
	// -- xsp -------------------------------------------------------
	
	if(ipol==0){
	  cout<<"\nOriginal Moment: \\int_0^1 * x * ( s(x,Q) + sbar(x,Q) )"<<endl;
	}
	if(ipol==1){
	  cout<<"\nOriginal Moment: \\int_0^1 * ( Delta_s(x,Q) + Delta_sbar(x,Q) )"<<endl;
	}
	cout<<"Mean +- Std = "<<unp_mom_xsp_mean<<"  "<<unp_mom_xsp_std<<
	  "  ( "<<100*unp_mom_xsp_std/unp_mom_xsp_mean<<" % )"<<endl;
	
	cout<<"Reweighted Moment: "<<endl;
	
	// Compute mean and variance
	double unp_mom_xsp_mean_rw=0;
	double unp_mom_xsp_std_rw=0;
	
	for(int irep=1;irep<=nrep;irep++){
	  unp_mom_xsp_mean_rw+=weights[irep]*unp_mom_rep_xsp[irep]/nrep;
	  unp_mom_xsp_std_rw+=weights[irep]*pow(unp_mom_rep_xsp[irep],2.0)/nrep;
	}
	unp_mom_xsp_std_rw=pow(unp_mom_xsp_std_rw - pow(unp_mom_xsp_mean_rw,2.0),0.5);
	
	cout<<"Mean +- Std = "<<unp_mom_xsp_mean_rw<<"  "<<unp_mom_xsp_std_rw<<
	  "  ( "<<100*unp_mom_xsp_std_rw/unp_mom_xsp_mean_rw<<" % )"<<endl;
	
	// -- xg -------------------------------------------------------
	
	if(ipol==0){
	  cout<<"\nOriginal Moment: \\int_0^1 x * ( g(x,Q) )"<<endl;
	}
	if(ipol==1){
	  cout<<"\nOriginal Moment: \\int_0^1 * x * ( Delta_u(x,Q) - Delta_ubar(x,Q) - Delta_d(x,Q) + Delta_dbar(x,Q) )"<<endl;
	}
	cout<<"Mean +- Std = "<<unp_mom_xg_mean<<"  "<<unp_mom_xg_std<<
	  "  ( "<<100*unp_mom_xg_std/unp_mom_xg_mean<<" % )"<<endl;
	
	cout<<"Reweighted Moment: "<<endl;
	
	// Compute mean and variance
	double unp_mom_xg_mean_rw=0;
	double unp_mom_xg_std_rw=0;
	
	for(int irep=1;irep<=nrep;irep++){
	  unp_mom_xg_mean_rw+=weights[irep]*unp_mom_rep_xg[irep]/nrep;
	  unp_mom_xg_std_rw+=weights[irep]*pow(unp_mom_rep_xg[irep],2.0)/nrep;
	}
	unp_mom_xg_std_rw=pow(unp_mom_xg_std_rw - pow(unp_mom_xg_mean_rw,2.0),0.5);
	
	cout<<"Mean +- Std = "<<unp_mom_xg_mean_rw<<"  "<<unp_mom_xg_std_rw<<
	  "  ( "<<100*unp_mom_xg_std_rw/unp_mom_xg_mean_rw<<" % )"<<endl;
	
	// -- xupmdp -------------------------------------------------------
	
	if(ipol==0){
	  cout<<"\nOriginal Moment: \\int_0^1 x * ( u(x,Q)+ubar(x,Q) - d(x,Q)+dbar(x,Q) )"<<endl;
	}
	if(ipol==1){
	  cout<<"\nOriginal Moment: \\int_0^1 *( Delta_u(x,Q) + Delta_ubar(x,Q) - Delta_d(x,Q) - Delta_dbar(x,Q) )" <<endl;
	}
	cout<<"Mean +- Std = "<<unp_mom_xupmdp_mean<<"  "<<unp_mom_xupmdp_std<<
	  "  ( "<<100*unp_mom_xupmdp_std/unp_mom_xupmdp_mean<<" % )"<<endl;
	
	cout<<"Reweighted Moment: "<<endl;
	
	// Compute mean and variance
	double unp_mom_xupmdp_mean_rw=0;
	double unp_mom_xupmdp_std_rw=0;
	
	for(int irep=1;irep<=nrep;irep++){
	  unp_mom_xupmdp_mean_rw+=weights[irep]*unp_mom_rep_xupmdp[irep]/nrep;
	  unp_mom_xupmdp_std_rw+=weights[irep]*pow(unp_mom_rep_xupmdp[irep],2.0)/nrep;
	}
	unp_mom_xupmdp_std_rw=pow(unp_mom_xupmdp_std_rw - pow(unp_mom_xupmdp_mean_rw,2.0),0.5);
	
	cout<<"Mean +- Std = "<<unp_mom_xupmdp_mean_rw<<"  "<<unp_mom_xupmdp_std_rw<<
	  "  ( "<<100*unp_mom_xupmdp_std_rw/unp_mom_xupmdp_mean_rw<<" % )"<<endl;
	
	//*******************************************************
	// Compute the reweighted PDFs
	//*******************************************************

	cout<<"\n\n******************************"<<endl;
	cout<<"   Computing Plots  "<<endl;
	cout<<"*******************************"<<endl;
	
	// Relative PDF errors
	std::vector<TGraph*> cvLuxRel;
	
	// Reference
	double *CVpdf = new double[nx];
	
	// Initialization
	for (int w = 0; w < nx; w++) CVpdf[w] = 1.0;
	
	double *x = new double[nx];
	
	// Relative uncertainties
	cvLuxRel.push_back(new TGraph(nx));
	cvLuxRel.push_back(new TGraph(nx));
	
	double *cv = new double[nx];
	double *err = new double[nx];
	double *cv_rw = new double[nx];
	double *err_rw = new double[nx];
	
	for (int ix = 0; ix < nx; ix++)
	  {
	    
	    double dexp=double(ix)/(nx-1);
	    
	    x[ix] = xmin*pow(xmax/xmin,dexp);
	    
	    cv[ix]  = 0;
	    err[ix]  = 0;
	    cv_rw[ix]  = 0;
	    err_rw[ix]  = 0;
	    for (int irep = 1; irep <= nrep; irep++)
	      {
		LHAPDF::initPDF(irep);
		double xpdflh;
		
		// Gluon
		if(l==0){
		  xpdflh=LHAPDF::xfx(x[ix],Q,0);
		}
		// Total up quark
		if(l==1){
		  xpdflh=LHAPDF::xfx(x[ix],Q,2) + LHAPDF::xfx(x[ix],Q,-2);
		}
		// Total down quark
		if(l==2){
		  xpdflh=LHAPDF::xfx(x[ix],Q,1) + LHAPDF::xfx(x[ix],Q,-1);
		}
		// Total strange quark
		if(l==3){
		  xpdflh=LHAPDF::xfx(x[ix],Q,3) + LHAPDF::xfx(x[ix],Q,-3);
		}
		// upmdp
		if(l==4){
		  xpdflh=LHAPDF::xfx(x[ix],Q,2) + LHAPDF::xfx(x[ix],Q,-2)
		    - ( LHAPDF::xfx(x[ix],Q,1) + LHAPDF::xfx(x[ix],Q,-1) );
		}
		
		cv[ix]+=xpdflh/nrep;
		err[ix]+=pow(xpdflh,2.0)/nrep;
		
		cv_rw[ix]+=weights[irep] * xpdflh/nrep;
		err_rw[ix]+=weights[irep] * pow(xpdflh,2.0)/nrep;
	      }

	    /*
	    cout<<"x = "<<x[ix]<<endl;
	    cout<<"cv, err = "<<cv[ix]<<" "<<err[ix]<<endl;
	    cout<<"cv_rw, err_rw = "<<cv_rw[ix]<<" "<<err_rw[ix]<<endl;
	    */
	    
	    err[ix] = pow(err[ix]-pow(cv[ix],2.0),0.5);
	    err_rw[ix] = pow(err_rw[ix]-pow(cv_rw[ix],2.0),0.5);
	    
	    for(int iPdf=0;iPdf<2;iPdf++){

	      // Take absolute value to avoid negative centrak values
	      // Specially relevant for the polarized case
	      // Maybe show instead the absolute PDF uncertainty
	      if(iPdf==0){
		//cvLuxRel[iPdf]->SetPoint(ix, x[ix], 100*err[ix]/fabs( cv[ix]) );
		cvLuxRel[iPdf]->SetPoint(ix, x[ix], err[ix] );
	      }
	      if(iPdf==1){
		//cvLuxRel[iPdf]->SetPoint(ix, x[ix], 100*err_rw[ix]/fabs( cv_rw[ix]) );
		cvLuxRel[iPdf]->SetPoint(ix, x[ix], err_rw[ix] );
	      }
	    }
	    
	  }
	
	for(int iPdf=0;iPdf<2;iPdf++){
	  if(ipol==0){
	    if(l==0)cvLuxRel[iPdf]->SetTitle("#delta( g ) for Q^{2}=4 GeV^{2}, NNPDF3.1 NNLO");
	    if(l==1)cvLuxRel[iPdf]->SetTitle("#delta( u^{+} ) ) for Q^{2}=4 GeV^{2}, NNPDF3.1 NNLO");
	    if(l==2)cvLuxRel[iPdf]->SetTitle("#delta( d^{+} ) for Q^{2}=4 GeV^{2}, NNPDF3.1 NNLO");
	    if(l==3)cvLuxRel[iPdf]->SetTitle("#delta( s^{+} ) for Q^{2}=4 GeV^{2}, NNPDF3.1 NNLO");
	    if(l==4)cvLuxRel[iPdf]->SetTitle("#delta( u^{+}-d^{+} ) for Q^{2}=4 GeV^{2}, NNPDF3.1 NNLO");
	  }
	  if(ipol==1){
	    if(l==0)cvLuxRel[iPdf]->SetTitle("#delta( #Delta g ) for Q^{2}=4 GeV^{2}, NNPDFpol1.1");
	    if(l==1)cvLuxRel[iPdf]->SetTitle("#delta( #Delta u^{+} ) ) for Q^{2}=4 GeV^{2}, NNPDFpol1.1");
	    if(l==2)cvLuxRel[iPdf]->SetTitle("#delta( #Delta d^{+} ) for Q^{2}=4 GeV^{2}, NNPDFpol1.1");
	    if(l==3)cvLuxRel[iPdf]->SetTitle("#delta( #Delta s^{+} ) for Q^{2}=4 GeV^{2}, NNPDFpol1.1");
	    if(l==4)cvLuxRel[iPdf]->SetTitle("#delta( #Delta u^{+} - #Delta d^{+} ) for Q^{2}=4 GeV^{2}, NNPDFpol1.1");
	  }
	  cvLuxRel[iPdf]->SetLineWidth(4);

	  cvLuxRel[0]->SetLineColor(1);
	  cvLuxRel[0]->SetLineStyle(1);
	  if(i_scen==0){
	    cvLuxRel[1]->SetLineColor(kGreen+2);
	    cvLuxRel[1]->SetLineStyle(3);
	  }
	  if(i_scen==1){
	    cvLuxRel[1]->SetLineColor(kRed+2);
	    cvLuxRel[1]->SetLineStyle(2);
	  }
	  if(i_scen==2){
	    cvLuxRel[1]->SetLineColor(kBlue+2);
	    cvLuxRel[1]->SetLineStyle(4);
	  }
		  
	  if(iPdf==0)cvLuxRel[iPdf]->SetLineStyle(1);
	  cvLuxRel[iPdf]->GetXaxis()->SetTitle("       x  ");
	  cvLuxRel[iPdf]->GetXaxis()->SetTitleSize(0.05);
	  cvLuxRel[iPdf]->GetXaxis()->SetTitleOffset(1.1);
	  cvLuxRel[iPdf]->GetXaxis()->SetLabelSize(0.045);
	  cvLuxRel[iPdf]->GetXaxis()->SetLimits(5e-3, 0.65);
	  cvLuxRel[iPdf]->GetXaxis()->CenterTitle(true);
	  cvLuxRel[iPdf]->GetYaxis()->SetTitle("Percentage PDF uncertainty");
	  cvLuxRel[iPdf]->GetYaxis()->CenterTitle(true);
	  if(ipol==0){
	    if(l==0)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,10);	
	    if(l==1)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,10);
	    if(l==2)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,10);
	    if(l==3)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,35);
	    if(l==4)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,35);
	  }
	  if(ipol==1){
	    // For percentage PDF uncertainty
	    /*
	    if(l==0)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,100);	
	    if(l==1)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,30);
	    if(l==2)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,100);
	    if(l==3)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,350);
	    if(l==4)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,20);
	    */
	    // for absolute PDF uncertainty
	    if(l==0)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,0.15);	
	    if(l==1)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,0.035);
	    if(l==2)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,0.1);
	    if(l==3)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,0.035);
	    if(l==4)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,0.1);
	  }
	  cvLuxRel[iPdf]->GetYaxis()->SetTitleSize(0.04);
	  cvLuxRel[iPdf]->GetYaxis()->SetLabelSize(0.03);
	}

	cout<<"Saving plot "<<endl;
	cout<<"i_scen, l = "<<i_scen<<" "<<l<<endl;
	
	int iPdf=0;
	c[0]->cd();
	if(i_scen==0)cvLuxRel[iPdf]->Draw("aC");
	
	iPdf=1;
	c[0]->cd();
	cvLuxRel[iPdf]->Draw("C,same");
	
	if(i_scen==0)leg[0]->AddEntry(cvLuxRel[0], "no LQCD moms","L");
	if(ipol==0){
	  if(i_scen==0)leg[0]->AddEntry(cvLuxRel[1], "with LQCD moms, #delta_{L} = 1%","L");
	  if(i_scen==1)leg[0]->AddEntry(cvLuxRel[1], "with LQCD moms, #delta_{L} = 3%","L");
	  if(i_scen==2)leg[0]->AddEntry(cvLuxRel[1], "with LQCD moms, #delta_{L} = 5%","L");
	}
	if(ipol==1){
	  if(i_scen==0)leg[0]->AddEntry(cvLuxRel[1], "with LQCD moms, #delta_{L} = 10%","L");
	  if(i_scen==1)leg[0]->AddEntry(cvLuxRel[1], "with LQCD moms, #delta_{L} = 20%","L");
	  if(i_scen==2)leg[0]->AddEntry(cvLuxRel[1], "with LQCD moms, #delta_{L} = 30%","L");
	}
	  
      }//  End loop over scenearios

      leg[0]->Draw();

      // Save the files
      string filename= outfilename2[l]+".eps";
      c[0]->SaveAs(filename.c_str());
      filename= outfilename2[l]+".pdf";
      c[0]->SaveAs(filename.c_str());
      
      
    } //End loop over PDF flavours

  } // End loop over unpolarized vs polarized


  return 0;
}

