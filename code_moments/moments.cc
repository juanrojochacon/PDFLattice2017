
#include "momentsint.h"

int main(int argc, char **argv)
{  
  cout << "\n ****************************************\n";
  cout <<   " * Computation of PDF moments *\n";
  cout <<   " ****************************************\n\n";
  cout << endl;

  //! Relative error integration
  double const eps = 1e-6;

  // Max number of PDF error members
  int const nrepmax=101;

  // Scale at which moments are evaluated
  double const Q = 2; 

  LHAPDF::setVerbosity(0);
  MomIntegral *lum = new MomIntegral(eps);

  // Array of PDF sets
  vector<string> pdfname;

  // Type of PDF errors
  vector<string> pdfertype;


  ////////////////////////////////////
  // List of PDF sets for which we
  // want to compute moments
  
  // PDF4LHC2015
  pdfname.push_back("PDF4LHC15_nnlo_mc");
  // This is a Monte Carlo set
  pdfertype.push_back("mc");

  // NNPDF3.0 NNLO
  pdfname.push_back("NNPDF30_nnlo_as_0118");
  // This is a Monte Carlo set
  pdfertype.push_back("mc");

  // NNPDF3.1 NNLO
  pdfname.push_back("NNPDF31_nnlo_as_0118");
  // This is a Monte Carlo set
  pdfertype.push_back("mc");

  // CT14NNLO
  pdfname.push_back("CT14nnlo");
  // This is a Hessian set
  pdfertype.push_back("hessianCT");

  // MMHT2014NNLO
  pdfname.push_back("MMHT2014nnlo68cl");
  // This is also a Hessian set
  pdfertype.push_back("hessian");
  
  // ABMP16 NNLO
  pdfname.push_back("ABMP16_4_nnlo");
  // This is a symmetric Hessian set
  pdfertype.push_back("symhessian");

  double cv_ref=0;

  for (size_t iPdf = 0; iPdf < pdfname.size(); iPdf++){

    // Initialize PDF set
    LHAPDF::initPDFSet(pdfname[iPdf]);

    // Get number of replicas
    int repfinal = LHAPDF::numberPDF();

      // Begin loop over replicas
      double pdfmom[nrepmax]={0.0};
      for (int irep = 0; irep <= repfinal; irep++){
	
	// Initialize replica
	LHAPDF::initPDF(irep);
	
	// Get the momentum fraction
	pdfmom[irep] = lum->getMom(Q);
		
      } // end loop over replicas

      double cv=0;
      double err=0;

      
      if(pdfertype[iPdf]=="mc"){
	// Now compute 68% CL
	
	vector<double> values;
	for (int nrep = 1; nrep <= repfinal; nrep++)
	  {
	    values.push_back(pdfmom[nrep]);
	  }
	// Now sort the array
	sort(values.begin(),values.end());
	
	// Compute the central value as the midpoint of the 68% CL interval
	int rep_min=repfinal*0.16;
	int rep_max=repfinal*0.84;
	cv= 0.5*(values.at(rep_min)+values.at(rep_max));
	err= 0.5*(values.at(rep_max)-values.at(rep_min));
	cv_ref=cv;

      }
      
      if(pdfertype[iPdf]=="hessianCT"){
	cv=pdfmom[0];
	for (int nrep = 1; nrep <= repfinal/2; nrep++)
	  {
	    err += pow(pdfmom[2*nrep-1]-pdfmom[2*nrep],2.0 );
	  }
	err=0.5 *pow(err,0.5)/1.642;

      }

      if(pdfertype[iPdf]=="hessian"){
	cv=pdfmom[0];
	for (int nrep = 1; nrep <= repfinal/2; nrep++)
	  {
	    err += pow(pdfmom[2*nrep-1]-pdfmom[2*nrep],2.0 );
	  }
	err=0.5 *pow(err,0.5);

      }

      if(pdfertype[iPdf]=="symhessian"){
	cv=pdfmom[0];
	for (int nrep = 1; nrep <= repfinal; nrep++)
	  {
	    err += pow(pdfmom[nrep]-pdfmom[0],2.0 );
	  }
	err=pow(err,0.5);

      }
	
      cout<<"pdfset = "<<pdfname[iPdf]<<endl;
      
      cout<<" <x(u-d)> = "<<cv<<" +- "<<err<<" ( "<< 100*err/cv  <<"% )  (shift from cv ref = "<<
	100*(cv-cv_ref)/cv_ref<<"% )"<<endl;
      
      
  } // End loop over PDF sets
  
  
  return 0;
}

