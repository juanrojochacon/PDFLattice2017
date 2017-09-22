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
  cout << "\n *********************************************************\n";
  cout <<   " * PDF reweighting with Lattice QCD x-space calculations *\n";
  cout <<   " *******************************************************\n\n";
  cout << endl;

  // Loop over the unpolarized and polarized cases in turn
  int const NrepUnpol=1000;
  int const NrepPol=100;
  int nrep=0;

  // Values of x where we assume that lattice-QCD will provide info
  int const nx=5;
  double const xl[nx]={0.65,0.70,0.75,0.80,0.85};
  double const Q=2; // that is 4 GeV2
  
  double xumxd_data[nx]={0.0};
  double xubmxdb_data[nx]={0.0};
  
  // Start loop between polarized and unpolarized
  for(int ipol=0;ipol<2;ipol++){

    // Set number of replicas
    if(ipol==0) nrep=NrepUnpol;
    if(ipol==1) nrep=NrepPol;

    // Define the pseudo-data: xu-xd and xubar-xdbar
    // for a bunch of values of x at large x
    // Define the PDF prior
    string pdfset;
    if(ipol==0) pdfset="NNPDF31_nnlo_as_0118_1000";
    if(ipol==1) pdfset="NNPDFpol11_100";
    // Initialization
    LHAPDF::initPDFSet(pdfset);
    // Set low verbosity as usual
    LHAPDF::setVerbosity(0);
    
    // Compute the pseudo-data
    for(int ix=0;ix<nx;ix++){
      
      // compute xu-xd
      double sum=0;
      for (int irep = 1; irep <= nrep; irep++)
	{
	  LHAPDF::initPDF(irep);
	  double xpdflh = LHAPDF::xfx(xl[ix],Q,2) + LHAPDF::xfx(xl[ix],Q,1) ;
	  sum += xpdflh/nrep;
	}
      xumxd_data[ix] = sum;

      // compute xub-xdb
      sum=0;
      for (int irep = 1; irep <= nrep; irep++)
	{
	  LHAPDF::initPDF(irep);
	  double xpdflh = LHAPDF::xfx(xl[ix],Q,-2) + LHAPDF::xfx(xl[ix],Q,-1) ;
	  sum += xpdflh/nrep;
	}
      xubmxdb_data[ix] = sum;
    }
    
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
    pdfname.push_back("xu");
    pdfname.push_back("xd");
    pdfname.push_back("xubar");
    pdfname.push_back("xdbar");
    
    std::vector<std::string> titleY;
    titleY.push_back("u ( x, Q^{2} = 4 GeV^{2} ) ");
    titleY.push_back("d ( x, Q^{2} = 4 GeV^{2} ) ");
    titleY.push_back("ubar ( x, Q^{2} = 4 GeV^{2} ) ");
    titleY.push_back("dbar ( x, Q^{2} = 4 GeV^{2} ) ");
    
      //double const xmin=1e-3;
      double const xmin=0.1;
      double const xmax=0.9;
      
      int const nx_p=50;
      
      std::vector<std::string> outfilename;
      std::vector<std::string> outfilename2;
      if(ipol==0){
	for(unsigned i=0;i< pdfname.size();i++){
	  outfilename.push_back(pdfname.at(i)+"-unpol-lattice-xdata");
	  outfilename2.push_back(pdfname.at(i)+"-unpol-lattice-relerr-xdata");
	}
      }
      if(ipol==1){
	for(unsigned i=0;i< pdfname.size();i++){
	  outfilename.push_back(pdfname.at(i)+"-pol-lattice-xdata");
	  outfilename2.push_back(pdfname.at(i)+"-pol-lattice-relerr-xdata");
	}
      }
      
      // Open file to save information on the effective number of replicas
      ofstream neffout;
      if(ipol==0) neffout.open("neff_unpol.res");
      if(ipol==1) neffout.open("neff_pol.res");
      
      
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
	  //c[0]->SetLogx();
	  
	  std::vector<TLegend*> leg;
	  leg.push_back(new TLegend(0.11,0.11,0.45,0.42));
	  leg[0]->SetLineStyle(1);
	  leg[0]->SetBorderSize(1);
	  leg[0]->SetFillColor(kWhite);
	  
	  //************************************************
	  // Compute the values of the chi2
	  //************************************************
	  
	  int const n_scenarios=3;
	  double total_latQCD_sysunc[n_scenarios]={0.0};
	  
	  // Unpolarized moments
	  if(ipol==0){
	    total_latQCD_sysunc[0]=12;
	    total_latQCD_sysunc[1]=6;
	    total_latQCD_sysunc[2]=3;
	  }

	  // Polarized moments (typically affected by larger uncertainties)
	  if(ipol==1){
	    total_latQCD_sysunc[0]=12;
	    total_latQCD_sysunc[1]=6;
	    total_latQCD_sysunc[2]=3;
	  }

	  // Loop over the number of scenarios
	  for(int i_scen=0;i_scen <n_scenarios;i_scen++ ){
	    
	    // Compute the chi2 for each replica
	    double chi2[NrepUnpol+1] = {0.0};

	    int const ndatfit=2*nx;
	    
	    for(int irep=1;irep<=nrep;irep++){

	      for(int ix=0;ix<nx;ix++){

		// Compute xu-xd for this replica
		LHAPDF::initPDF(irep);
		double xpdflh = LHAPDF::xfx(xl[ix],Q,2) + LHAPDF::xfx(xl[ix],Q,1) ;
		// cout<<irep<<" "<<xl[ix]<<" "<<xpdflh<<" "<<xumxd_data[ix]<<endl;
		chi2[irep] += pow( xpdflh - xumxd_data[ix],2.0)/
		  pow(total_latQCD_sysunc[i_scen]*xumxd_data[ix] ,2.0);
		
		// Compute xub-xdb for this replica
		LHAPDF::initPDF(irep);
		xpdflh = LHAPDF::xfx(xl[ix],Q,-2) + LHAPDF::xfx(xl[ix],Q,-1) ;
		
		chi2[irep] += pow( xpdflh - xubmxdb_data[ix],2.0)/
		  pow(total_latQCD_sysunc[i_scen]*xubmxdb_data[ix] ,2.0);

	      } // end loop over ix
	      // Normalization to the number of data points
	      chi2[irep]=chi2[irep] / ndatfit;
	      
	    }// end loop over replicas
	    
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
	    
	    if(l==0){
	      neffout << i_scen << " \t "<< neff<<endl;
	    }


	    
	    //*******************************************************
	    // Compute the reweighted PDFs
	    //*******************************************************
	    
	    cout<<"\n\n******************************"<<endl;
	    cout<<"   Computing Plots  "<<endl;
	    cout<<"*******************************"<<endl;
	    
	    // Relative PDF errors
	    std::vector<TGraph*> cvLuxRel;
	    
	    // Reference
	    double *CVpdf = new double[nx_p];
	    
	    // Initialization
	    for (int w = 0; w < nx_p; w++) CVpdf[w] = 1.0;
	    
	    double *x = new double[nx_p];
	    
	    // Relative uncertainties
	    cvLuxRel.push_back(new TGraph(nx_p));
	    cvLuxRel.push_back(new TGraph(nx_p));
	    
	    double *cv = new double[nx_p];
	    double *err = new double[nx_p];
	    double *cv_rw = new double[nx_p];
	    double *err_rw = new double[nx_p];
	    
	    for (int ix = 0; ix < nx_p; ix++)
	      {
		
		double dexp=double(ix)/(nx_p-1);
		
		x[ix] = xmin*pow(xmax/xmin,dexp);
		
		cv[ix]  = 0;
		err[ix]  = 0;
		cv_rw[ix]  = 0;
		err_rw[ix]  = 0;
		for (int irep = 1; irep <= nrep; irep++)
		  {
		    LHAPDF::initPDF(irep);
		    double xpdflh;
		    
		    // up quark
		    if(l==0){
		      xpdflh=LHAPDF::xfx(x[ix],Q,2); 
		    }
		    // down quark
		    if(l==1){
		      xpdflh=LHAPDF::xfx(x[ix],Q,1) ;
		    }
		    // up antiquark
		    if(l==2){
		      xpdflh=LHAPDF::xfx(x[ix],Q,-2);
		    }
		    // down antiquark
		    if(l==3){
		      xpdflh=LHAPDF::xfx(x[ix],Q,-1);
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
		    if(ipol==0){
		      cvLuxRel[iPdf]->SetPoint(ix, x[ix], 1.0 );
		    }
		    if(ipol==1){
		      // cvLuxRel[iPdf]->SetPoint(ix, x[ix], err[ix] );
		      cvLuxRel[iPdf]->SetPoint(ix, x[ix], 1.0 );
		    }
		    
		  }
		  if(iPdf==1){
		    if(ipol==0){
		      //		      cout<<ix<<" "<< ( 100*err_rw[ix]/fabs( cv_rw[ix])  )/( 100*err[ix]/fabs( cv[ix]) )  <<endl; 
		      cvLuxRel[iPdf]->SetPoint(ix, x[ix], ( err_rw[ix] / err[ix] )  );
		    }
		    if(ipol==1){
		      // cvLuxRel[iPdf]->SetPoint(ix, x[ix], err_rw[ix] );
		      cvLuxRel[iPdf]->SetPoint(ix, x[ix], ( err_rw[ix] / err[ix] )  );
		    }
		  }
		}
		
	      }
	    
	    for(int iPdf=0;iPdf<2;iPdf++){
	      if(ipol==0){
		if(l==0)cvLuxRel[iPdf]->SetTitle("#delta( u ) @ Q^{2}=4 GeV^{2}, NNPDF3.1");
		if(l==1)cvLuxRel[iPdf]->SetTitle("#delta( d ) @ Q^{2}=4 GeV^{2}, NNPDF3.1");
		if(l==2)cvLuxRel[iPdf]->SetTitle("#delta( #bar{u} ) @ Q^{2}=4 GeV^{2}, NNPDF3.1");
		if(l==3)cvLuxRel[iPdf]->SetTitle("#delta( #bar{d} ) @ Q^{2}=4 GeV^{2}, NNPDF3.1");
	      }
	      if(ipol==1){
		if(l==0)cvLuxRel[iPdf]->SetTitle("#delta( #Delta u ) @ Q^{2}=4 GeV^{2}, NNPDFpol1.1");
		if(l==1)cvLuxRel[iPdf]->SetTitle("#delta( #Delta d ) @ Q^{2}=4 GeV^{2}, NNPDFpol1.1");
		if(l==2)cvLuxRel[iPdf]->SetTitle("#delta( #Delta #bar{u} ) @ Q^{2}=4 GeV^{2}, NNPDFpol1.1");
		if(l==3)cvLuxRel[iPdf]->SetTitle("#delta( #Delta #bar{d} ) @ Q^{2}=4 GeV^{2}, NNPDFpol1.1");
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
	      cvLuxRel[iPdf]->GetXaxis()->SetTitleOffset(0.9);
	      cvLuxRel[iPdf]->GetXaxis()->SetLabelSize(0.045);
	      cvLuxRel[iPdf]->GetXaxis()->SetLimits(0.3, 0.9);
	      cvLuxRel[iPdf]->GetXaxis()->CenterTitle(true);
	      if(ipol==0)cvLuxRel[iPdf]->GetYaxis()->SetTitle("PDF uncertainty (ratio to NNPDF3.1)");
	      if(ipol==1)cvLuxRel[iPdf]->GetYaxis()->SetTitle("PDF uncertainty (ratio to NNPDFpol1.1)");
	      cvLuxRel[iPdf]->GetYaxis()->CenterTitle(true);
	      if(ipol==0){
		// if(l==0)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,10);	
		//if(l==1)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,5);
		//if(l==2)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,8);
		//if(l==3)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,35);
		// if(l==4)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,35);
		
		if(l==0)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0.4,1.12);	
		if(l==1)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0.4,1.12);
		if(l==2)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0.4,1.12);
		if(l==3)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0.4,1.12);
		if(l==4)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0.4,1.12);
		if(l==5)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0.4,1.12);
	      }
	      if(ipol==1){

			
		if(l==0)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0.4,1.12);	
		if(l==1)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0.4,1.12);
		if(l==2)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0.4,1.12);
		if(l==3)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0.4,1.12);
		if(l==4)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0.4,1.12);
		if(l==5)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0.4,1.12);

		
		// For percentage PDF uncertainty
		/*
		  if(l==0)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,100);	
		  if(l==1)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,30);
		  if(l==2)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,100);
		  if(l==3)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,350);
		  if(l==4)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,20);
		*/
		// for absolute PDF uncertainty
		/*
		if(l==0)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,0.20);	
		if(l==1)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,0.025);
		if(l==2)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,0.03);
		if(l==3)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,0.035);
		if(l==4)cvLuxRel[iPdf]->GetYaxis()->SetRangeUser(0,0.1);
		*/
	      }
	      cvLuxRel[iPdf]->GetYaxis()->SetTitleSize(0.044);
	      cvLuxRel[iPdf]->GetYaxis()->SetTitleOffset(0.9);
	      cvLuxRel[iPdf]->GetYaxis()->SetLabelSize(0.034);
	    }
	    
	    int iPdf=0;
	    c[0]->cd();
	    if(i_scen==0)cvLuxRel[iPdf]->Draw("aL");
	    
	    iPdf=1;
	    c[0]->cd();
	    cvLuxRel[iPdf]->Draw("C,same");
	    
	    if(i_scen==0)leg[0]->AddEntry(cvLuxRel[0], "no LQCD moms","L");
	    if(ipol==0){
	      if(i_scen==0)leg[0]->AddEntry(cvLuxRel[1], "w LattQCD x-sp, Scen D","L");
	      if(i_scen==1)leg[0]->AddEntry(cvLuxRel[1], "w LattQCD x-sp, Scen E","L");
	      if(i_scen==2)leg[0]->AddEntry(cvLuxRel[1], "w LattQCD x-sp, Scen F","L");
	    }
	    if(ipol==1){
	      if(i_scen==0)leg[0]->AddEntry(cvLuxRel[1], "w LattQCD x-sp, Scen D","L");
	      if(i_scen==1)leg[0]->AddEntry(cvLuxRel[1], "w LattQCD x-sp, Scen E","L");
	      if(i_scen==2)leg[0]->AddEntry(cvLuxRel[1], "w LattQCD x-sp, Scen F","L");
	    }
	    
	  }//  End loop over scenearios
	  
	  leg[0]->Draw();
	  
	  // Save the files
	  string filename= outfilename2[l]+"-xspace.eps";
	  c[0]->SaveAs(filename.c_str());
	  filename= outfilename2[l]+"-xspace.pdf";
	  c[0]->SaveAs(filename.c_str());
	  
	  
	} //End loop over PDF flavours
      
      neffout.close();
      
    } // End loop over unpolarized vs polarized

 
    
    return 0;
}

