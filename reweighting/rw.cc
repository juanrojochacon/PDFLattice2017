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
  // First the unpolarized moments - 1000 replicas with NNPDF3.1
  int const NrepUnpol=1000;

  ifstream in1;
  // This is \int_0^1 x * ( d(x,Q) + dbar(x,Q) ), Eq A.4 of whitepaper
  double unp_mom_rep_xdp[NrepUnpol+1]={0.0};
  in1.open("../code_moments/replicas/res/unp_mom_rep_xdp.res");
  for(int irep=1;irep<=NrepUnpol;irep++){
    int idum=0;
    in1>>idum>> unp_mom_rep_xdp[irep];
    if(idum!=irep){
      cout<<"Error when reading replica file";
      exit(-10);
    }
  }
  in1.close();

  // This is \int_0^1 x * ( u(x,Q) + ubar(x,Q) ), Eq A.4 of whitepaper
  double unp_mom_rep_xup[NrepUnpol+1]={0.0};
  in1.open("../code_moments/replicas/res/unp_mom_rep_xup.res");
  for(int irep=1;irep<=NrepUnpol;irep++){
    int idum=0;
    in1>>idum>> unp_mom_rep_xup[irep];
    if(idum!=irep){
      cout<<"Error when reading replica file";
      exit(-10);
    }
  }
  in1.close();

  // This is \int_0^1 x * ( s(x,Q) + sbar(x,Q) ), Eq A.4 of whitepaper
  double unp_mom_rep_xsp[NrepUnpol+1]={0.0};
  in1.open("../code_moments/replicas/res/unp_mom_rep_xsp.res");
  for(int irep=1;irep<=NrepUnpol;irep++){
    int idum=0;
    in1>>idum>> unp_mom_rep_xsp[irep];
    if(idum!=irep){
      cout<<"Error when reading replica file";
      exit(-10);
    }
  }
  in1.close();

  //************************************************
  // Compute means and variances
  //************************************************
  
  // Compute mean and variance
  double unp_mom_xdp_mean=0;
  double unp_mom_xdp_std=0;

  for(int irep=1;irep<=NrepUnpol;irep++){
    unp_mom_xdp_mean+=unp_mom_rep_xdp[irep]/NrepUnpol;
    unp_mom_xdp_std+=pow(unp_mom_rep_xdp[irep],2.0)/NrepUnpol;
  }
  unp_mom_xdp_std=pow(unp_mom_xdp_std - pow(unp_mom_xdp_mean,2.0),0.5);

  // Compute mean and variance
  double unp_mom_xup_mean=0;
  double unp_mom_xup_std=0;

  for(int irep=1;irep<=NrepUnpol;irep++){
    unp_mom_xup_mean+=unp_mom_rep_xup[irep]/NrepUnpol;
    unp_mom_xup_std+=pow(unp_mom_rep_xup[irep],2.0)/NrepUnpol;
  }
  unp_mom_xup_std=pow(unp_mom_xup_std - pow(unp_mom_xup_mean,2.0),0.5);

  // Compute mean and variance
  double unp_mom_xsp_mean=0;
  double unp_mom_xsp_std=0;

  for(int irep=1;irep<=NrepUnpol;irep++){
    unp_mom_xsp_mean+=unp_mom_rep_xsp[irep]/NrepUnpol;
    unp_mom_xsp_std+=pow(unp_mom_rep_xsp[irep],2.0)/NrepUnpol;
  }
  unp_mom_xsp_std=pow(unp_mom_xsp_std - pow(unp_mom_xsp_mean,2.0),0.5);

  
  //************************************************
  // Compute the values of the chi2
  //************************************************
  
  // Conservative: 5% total systematics
  // Optimistic: 1% total systematics
  int const n_scenarios=2;
  double const total_latQCD_sysunc[n_scenarios]={0.01, 0.05};
  
    // Compute the chi2 for each replica
  double chi2[NrepUnpol+1] = {0.0};

  // Fix here the number of moments that are included in the chi2 computation
  int const ndatfit=3;
  for(int irep=1;irep<=NrepUnpol;irep++){

    // xdp
    chi2[irep] += pow( unp_mom_xdp_mean - unp_mom_rep_xdp[irep] ,2.0)/
	     pow(total_latQCD_sysunc[0]*unp_mom_xdp_mean ,2.0);

    // xup
    chi2[irep] += pow( unp_mom_xup_mean - unp_mom_rep_xup[irep] ,2.0)/
	     pow(total_latQCD_sysunc[0]*unp_mom_xup_mean ,2.0);

    // xsp
    chi2[irep] += pow( unp_mom_xsp_mean - unp_mom_rep_xsp[irep] ,2.0)/
	     pow(total_latQCD_sysunc[0]*unp_mom_xsp_mean ,2.0);

    
    chi2[irep]/=ndatfit;
    cout<<irep<<" "<<chi2[irep]<<endl;
  }
     
  // Compute the weights
  double weights[NrepUnpol+1]={0.0};
  double sum=0.0;
  int const irep_ref=1;
  for(int irep=1;irep<=NrepUnpol;irep++){
    double logw = ( double(ndatfit-1)/2.0 ) * log( chi2[irep] * ndatfit )
      - ( chi2[irep] * ndatfit ) / 2.0;
    //cout<<"logw"<<" "<<logw<<endl;
    logw -= ( (ndatfit-1)/2 ) * log( chi2[irep_ref] * ndatfit )
      - ( chi2[irep_ref] * ndatfit ) / 2.0;
    weights[irep]=exp(logw);
    //cout<<irep<<" "<<logw<<" "<<weights[irep]<<endl;
        
    sum += weights[irep] / NrepUnpol;
  }
  
  // Check normalization
  double check=0;
  for(int irep=1;irep<=NrepUnpol;irep++){
    weights[irep]/=sum;
    check+=weights[irep];
    cout<<"irep, weights =  "<<irep<<"  "<<weights[irep]<<endl;
  }
  std::cout<<"\n Sum of weights = "<<check<<"\n"<<std::endl;
  if(fabs(check-NrepUnpol)>1e-3){
    cout<<"Problem with the computation of the weights"<<endl;
    exit(-10);
  }
  
  // Now compute the effective number of replicas
  double shannon=0.0;
  for(int irep=1;irep<=NrepUnpol;irep++){
    if(weights[irep]>1e-8){
      shannon += weights[irep] * log(double(NrepUnpol)/weights[irep]) / NrepUnpol;
    }
  }
  double const neff=exp(shannon);
  std::cout<<"Effective number of replicas = "<<neff<<"\n "<<std::endl;

  
  // ******************************************************
  // Now compute the reweighted moments
  // ******************************************************

  // -- xdp -------------------------------------------------------
  
  cout<<"\n Orginal Moment: \\int_0^1 x * ( d(x,Q) + dbar(x,Q) )"<<endl;
  cout<<"Mean +- Std = "<<unp_mom_xdp_mean<<"  "<<unp_mom_xdp_std<<
    "  ( "<<100*unp_mom_xdp_std/unp_mom_xdp_mean<<" % )"<<endl;

  cout<<"\n Reweighted Moment: "<<endl;

  // Compute mean and variance
  double unp_mom_xdp_mean_rw=0;
  double unp_mom_xdp_std_rw=0;

  for(int irep=1;irep<=NrepUnpol;irep++){
    unp_mom_xdp_mean_rw+=weights[irep]*unp_mom_rep_xdp[irep]/NrepUnpol;
    unp_mom_xdp_std_rw+=weights[irep]*pow(unp_mom_rep_xdp[irep],2.0)/NrepUnpol;
  }
  unp_mom_xdp_std_rw=pow(unp_mom_xdp_std_rw - pow(unp_mom_xdp_mean_rw,2.0),0.5);
  
  cout<<"Mean +- Std = "<<unp_mom_xdp_mean_rw<<"  "<<unp_mom_xdp_std_rw<<
    "  ( "<<100*unp_mom_xdp_std_rw/unp_mom_xdp_mean_rw<<" % )"<<endl;

  // -- xup -------------------------------------------------------
  
  cout<<"\n Original Moment: \\int_0^1 x * ( u(x,Q) + ubar(x,Q) )"<<endl;
  cout<<"Mean +- Std = "<<unp_mom_xup_mean<<"  "<<unp_mom_xup_std<<
    "  ( "<<100*unp_mom_xup_std/unp_mom_xup_mean<<" % )"<<endl;

  cout<<"\n Reweighted Moment: "<<endl;
  
  // Compute mean and variance
  double unp_mom_xup_mean_rw=0;
  double unp_mom_xup_std_rw=0;

  for(int irep=1;irep<=NrepUnpol;irep++){
    unp_mom_xup_mean_rw+=weights[irep]*unp_mom_rep_xup[irep]/NrepUnpol;
    unp_mom_xup_std_rw+=weights[irep]*pow(unp_mom_rep_xup[irep],2.0)/NrepUnpol;
  }
  unp_mom_xup_std_rw=pow(unp_mom_xup_std_rw - pow(unp_mom_xup_mean_rw,2.0),0.5);
  
  cout<<"Mean +- Std = "<<unp_mom_xup_mean_rw<<"  "<<unp_mom_xup_std_rw<<
    "  ( "<<100*unp_mom_xup_std_rw/unp_mom_xup_mean_rw<<" % )"<<endl;

  // -- xsp -------------------------------------------------------
    
  cout<<"\n Original Moment: \\int_0^1 x * ( s(x,Q) + sbar(x,Q) )"<<endl;
  cout<<"Mean +- Std = "<<unp_mom_xsp_mean<<"  "<<unp_mom_xsp_std<<
    "  ( "<<100*unp_mom_xsp_std/unp_mom_xsp_mean<<" % )"<<endl;

  cout<<"\n Reweighted Moment: "<<endl;
  
  // Compute mean and variance
  double unp_mom_xsp_mean_rw=0;
  double unp_mom_xsp_std_rw=0;

  for(int irep=1;irep<=NrepUnpol;irep++){
    unp_mom_xsp_mean_rw+=weights[irep]*unp_mom_rep_xsp[irep]/NrepUnpol;
    unp_mom_xsp_std_rw+=weights[irep]*pow(unp_mom_rep_xsp[irep],2.0)/NrepUnpol;
  }
  unp_mom_xsp_std_rw=pow(unp_mom_xsp_std_rw - pow(unp_mom_xsp_mean_rw,2.0),0.5);
  
  cout<<"Mean +- Std = "<<unp_mom_xsp_mean_rw<<"  "<<unp_mom_xsp_std_rw<<
    "  ( "<<100*unp_mom_xsp_std_rw/unp_mom_xsp_mean_rw<<" % )"<<endl;

  

  //*******************************************************
  // Compute the reweighted PDFs
  //*******************************************************

  std::cout<<"\n ******************************************************** \n"<<std::endl;
  std::cout<<"      Plotting orig vs rw vs unw  "<<std::endl;
  std::cout<<"      Also relative decrease in PDF uncertainties  "<<std::endl;
  std::cout<<"\n ******************************************************* \n"<<std::endl;
  
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

  std::vector<std::string> titleY;
  titleY.push_back("g ( x, Q^{2} = 4 GeV^{2} ) ");
  titleY.push_back("u^{+} ( x, Q^{2} = 4 GeV^{2} ) ");
  titleY.push_back("d^{+} ( x, Q^{2} = 4 GeV^{2} ) ");
  titleY.push_back("s^{+} ( x, Q^{2} = 4 GeV^{2} ) ");
 
  const int histoFillStyle[] = { 1001, 3005, 3004};

  double const xmin=1e-3;
  double const xmax=0.90;

  int const nx=100;

  std::vector<std::string> outfilename;
  std::vector<std::string> outfilename2;
  for(unsigned i=0;i< pdfname.size();i++){
    outfilename.push_back(pdfname.at(i)+"-lattice");
    outfilename2.push_back(pdfname.at(i)+"-lattice-relerr");
  }
  
  // Define the PDF prior
  string pdfset="NNPDF31_nnlo_as_0118_1000";
  // Initialization
  LHAPDF::initPDFSet(pdfset);
  // Set low verbosity as usual
  LHAPDF::setVerbosity(0);
  
  for (size_t l = 0; l < pdfname.size(); l++)
    {
      
      std::vector<TCanvas*> c;
      
      // Ratio plots
      std::vector<TGraphErrors*> cvLux;
      std::vector<TGraphErrors*> cvLuxUp;
      std::vector<TGraphErrors*> cvLuxDn;
      std::vector<TGraphErrors*> cvLuxCT;

      // Relative PDF errors
      std::vector<TGraph*> cvLuxRel;
         
      // Reference
      double *CVpdf = new double[nx];
      
      // Initialization
      for (int w = 0; w < nx; w++) CVpdf[w] = 1.0;

      double *x = new double[nx];

      cvLux.push_back(new TGraphErrors(nx));
      cvLuxUp.push_back(new TGraphErrors(nx));
      cvLuxDn.push_back(new TGraphErrors(nx));
      cvLuxCT.push_back(new TGraphErrors(nx));

      cvLux.push_back(new TGraphErrors(nx));
      cvLuxUp.push_back(new TGraphErrors(nx));
      cvLuxDn.push_back(new TGraphErrors(nx));
      cvLuxCT.push_back(new TGraphErrors(nx));

      cvLux.push_back(new TGraphErrors(nx));
      cvLuxUp.push_back(new TGraphErrors(nx));
      cvLuxDn.push_back(new TGraphErrors(nx));
      cvLuxCT.push_back(new TGraphErrors(nx));

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
	  for (int irep = 1; irep <= NrepUnpol; irep++)
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
	      
	      cv[ix]+=xpdflh/NrepUnpol;
	      err[ix]+=pow(xpdflh,2.0)/NrepUnpol;
	      
	      cv_rw[ix]+=weights[irep] * xpdflh/NrepUnpol;
	      err_rw[ix]+=weights[irep] * pow(xpdflh,2.0)/NrepUnpol;
	    }
	  
	  err[ix] = pow(err[ix]-pow(cv[ix],2.0),0.5);
	  err_rw[ix] = pow(err_rw[ix]-pow(cv_rw[ix],2.0),0.5);
	  
	  for(int iPdf=0;iPdf<3;iPdf++){
	    
	    if(iPdf==0){
	      cvLux[iPdf]->SetPoint(ix, x[ix], cv[ix]);
	      cvLux[iPdf]->SetPointError(ix, 0.0, err[ix]);
	      cvLuxCT[iPdf]->SetPoint(ix, x[ix], cv[ix]);
	      cvLuxUp[iPdf]->SetPoint(ix, x[ix], (cv[ix]+err[ix]));
	      cvLuxDn[iPdf]->SetPoint(ix, x[ix], (cv[ix]-err[ix]));

	      cvLuxRel[iPdf]->SetPoint(ix, x[ix], 100*err[ix]/cv[ix]);
	      
	    }
	    if(iPdf==1){
	      cvLux[iPdf]->SetPoint(ix, x[ix], cv_rw[ix]);
	      cvLux[iPdf]->SetPointError(ix, 0.0, err_rw[ix]);
	      cvLuxCT[iPdf]->SetPoint(ix, x[ix], cv_rw[ix]);
	      cvLuxUp[iPdf]->SetPoint(ix, x[ix], (cv_rw[ix]+err_rw[ix]));
	      cvLuxDn[iPdf]->SetPoint(ix, x[ix], (cv_rw[ix]-err_rw[ix]));

	      cvLuxRel[iPdf]->SetPoint(ix, x[ix], 100*err_rw[ix]/cv_rw[ix]);
	      
	    }
	  }
	  
	}
      
      c.push_back(new TCanvas());
      c[0]->SetTickx();
      c[0]->SetTicky();
      c[0]->SetLogx();
      
      for(int iPdf=0;iPdf<3;iPdf++){
	cvLux[iPdf]->SetTitle("NNPDF3.1 NNLO");
	cvLux[iPdf]->SetLineWidth(2);
	cvLux[iPdf]->SetLineColor(color2[iPdf]);
	cvLux[iPdf]->SetLineStyle(2);
	cvLux[iPdf]->SetFillColor(color[iPdf]);
	cvLux[iPdf]->SetFillStyle(histoFillStyle[iPdf]);
	cvLux[iPdf]->GetXaxis()->SetTitle("       x  ");
	cvLux[iPdf]->GetXaxis()->SetTitleSize(0.05);
	cvLux[iPdf]->GetXaxis()->SetTitleOffset(1.0);
	cvLux[iPdf]->GetXaxis()->SetLabelSize(0.035);
	cvLux[iPdf]->GetXaxis()->SetLimits(1e-3, 0.8);
	cvLux[iPdf]->GetXaxis()->CenterTitle(true);
	cvLux[iPdf]->GetYaxis()->SetTitle(titleY[l].c_str());
	cvLux[iPdf]->GetYaxis()->CenterTitle(true);
	if(l==0)cvLux[iPdf]->GetYaxis()->SetRangeUser(-2,20);
	if(l!=0)cvLux[iPdf]->GetYaxis()->SetRangeUser(0,20);
	cvLux[iPdf]->GetYaxis()->SetTitleSize(0.05);
	cvLux[iPdf]->GetYaxis()->SetLabelSize(0.05);
	cvLuxCT[iPdf]->SetLineColor(color2[iPdf]);
	cvLuxCT[iPdf]->SetLineWidth(2);
	cvLuxCT[iPdf]->SetLineStyle(2);
	cvLuxUp[iPdf]->SetLineColor(color2[iPdf]);
	cvLuxUp[iPdf]->SetLineWidth(2);
	cvLuxUp[iPdf]->SetLineStyle(1);
	cvLuxDn[iPdf]->SetLineColor(color2[iPdf]);
	cvLuxDn[iPdf]->SetLineWidth(2);
	cvLuxDn[iPdf]->SetLineStyle(1);
      }
      
      int iPdf=0;
      cvLux[iPdf]->Draw("a3");
      cvLuxCT[iPdf]->Draw("l");
      cvLuxUp[iPdf]->Draw("l");
      cvLuxDn[iPdf]->Draw("l");
      
      iPdf=1;
      cvLux[iPdf]->Draw("3,same");
      cvLuxCT[iPdf]->Draw("l,same");
      cvLuxUp[iPdf]->Draw("l,same");
      cvLuxDn[iPdf]->Draw("l,same");

      vector<TLegend*> leg;
      leg.push_back(new TLegend(0.43,0.70,0.89,0.89));
      leg[0]->SetLineStyle(1);
      leg[0]->SetBorderSize(1);
      leg[0]->SetFillColor(kWhite);
      leg[0]->AddEntry(cvLux[0], "no LQCD moms","FL");
      leg[0]->AddEntry(cvLux[1], "with LQCD moms","FL");
      leg[0]->Draw();
      
      std::string filename= outfilename[l]+".eps";
      c[0]->SaveAs(filename.c_str());
      filename= outfilename[l]+".pdf";
      c[0]->SaveAs(filename.c_str());

      //////////////////////////////////////////////
      // Relative PDF uncertainties
      
      c.push_back(new TCanvas());
      c[1]->SetTickx();
      c[1]->SetTicky();
      c[1]->SetLogx();
      
      for(int iPdf=0;iPdf<2;iPdf++){
	if(l==0)cvLuxRel[iPdf]->SetTitle("#Delta( g(x,Q^{2}) ) for Q^{2}=4 GeV^{2}, NNPDF3.1 NNLO");
	if(l==1)cvLuxRel[iPdf]->SetTitle("#Delta( u^{+}(x,Q^{2}) ) for Q^{2}=4 GeV^{2}, NNPDF3.1 NNLO");
	if(l==2)cvLuxRel[iPdf]->SetTitle("#Delta( d^{+}(x,Q^{2}) ) for Q^{2}=4 GeV^{2}, NNPDF3.1 NNLO");
	if(l==3)cvLuxRel[iPdf]->SetTitle("#Delta( s^{+}(x,Q^{2}) ) for Q^{2}=4 GeV^{2}, NNPDF3.1 NNLO");
	cvLuxRel[iPdf]->SetLineWidth(4);
	cvLuxRel[iPdf]->SetLineColor(color2[iPdf]);
	cvLuxRel[iPdf]->SetLineStyle(2);
	if(iPdf==0)cvLuxRel[iPdf]->SetLineStyle(1);
	cvLuxRel[iPdf]->GetXaxis()->SetTitle("       x  ");
	cvLuxRel[iPdf]->GetXaxis()->SetTitleSize(0.05);
	cvLuxRel[iPdf]->GetXaxis()->SetTitleOffset(1.0);
	cvLuxRel[iPdf]->GetXaxis()->SetLabelSize(0.035);
	cvLuxRel[iPdf]->GetXaxis()->SetLimits(5e-3, 0.8);
	cvLuxRel[iPdf]->GetXaxis()->CenterTitle(true);
	cvLuxRel[iPdf]->GetYaxis()->SetTitle("Percentage PDF uncertainty");
	cvLuxRel[iPdf]->GetYaxis()->CenterTitle(true);
	if(l==0)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,15);	
	if(l==1)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,15);
	if(l==2)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,15);
	if(l==3)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,35);
	cvLuxRel[iPdf]->GetYaxis()->SetTitleSize(0.05);
	cvLuxRel[iPdf]->GetYaxis()->SetLabelSize(0.05);
      }
      
      iPdf=0;
      cvLuxRel[iPdf]->Draw("aC");
           
      iPdf=1;
      cvLuxRel[iPdf]->Draw("C,same");

      leg.push_back(new TLegend(0.11,0.70,0.45,0.89));
      leg[1]->SetLineStyle(1);
      leg[1]->SetBorderSize(1);
      leg[1]->SetFillColor(kWhite);
      leg[1]->AddEntry(cvLuxRel[0], "no LQCD moms","L");
      leg[1]->AddEntry(cvLuxRel[1], "with LQCD moms","L");
      leg[1]->Draw();
      
      filename= outfilename2[l]+".eps";
      c[1]->SaveAs(filename.c_str());
      filename= outfilename2[l]+".pdf";
      c[1]->SaveAs(filename.c_str());
      
    }

  
  return 0;
}

