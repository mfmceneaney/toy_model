#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooCategory.h"
using namespace RooFit;

std::vector<double> getFitParams(int fit_idx, int nevents) {

   // S e t u p   m o d e l
   // ---------------------

   // Define discrete variables and set state labels
   RooCategory h("h", "helicity");
   h.defineType("minus", -1);
   h.defineType("plus", 1);
 
   // Declare continuous variables x,mean,sigma with associated name, title, initial value and allowed range
   RooRealVar x("x", "x fit variable", -1.0, 1.0);
   RooRealVar y("y", "y fit variable", 0.0, 2.0*TMath::Pi());
   RooRealVar a0("a0", "a0", 0.00, -1.0, 1.0);
   RooRealVar a1("a1", "a1", 0.1, -1.0, 1.0);
   RooRealVar a2("a2", "a2", 0.00, -1.0, 1.0);
   
   // Build generic pdf in terms of fit variables and parameters
   RooGenericPdf gen("gen", "1.0+0.747*1.0*x[0]*x[1]*(cos(x[2])*x[3]+x[4]+cos(2.0*x[2])*x[5])", RooArgList(h,x,y,a0,a1,a2));
   
   // G e n e r a t e   e v e n t s
   // -----------------------------
 
   // Generate a dataset of 1000 events in x from gen
   std::unique_ptr<RooDataSet> data{gen.generate(RooArgSet(h,x,y), nevents)};
   // std::unique_ptr<RooDataSet> data1Dx{gen_x_1d.generate(RooArgSet(x), nevents)};
   // std::unique_ptr<RooDataSet> data1Dy{gen_y_1d.generate(RooArgSet(y), nevents)};

   // // Create binned clones
   // std::unique_ptr<RooDataHist> data_hist{data->binnedClone()};
   // std::unique_ptr<RooDataHist> data1Dx_hist{data1Dx->binnedClone()};
   // std::unique_ptr<RooDataHist> data1Dy_hist{data1Dy->binnedClone()};

   RooDataSet data_plus("data_hplus", "data_hplus", &*data, RooArgSet(h,x,y), "h>0");
   RooDataSet data_minus("data_hminus", "data_hminus", &*data, RooArgSet(h,x,y), "h<0");

   // U s e   B a s i c   R O O T   F i t t i n g
   // -------------------------------------------

   // Set parametres
   std::string fitformula = "1.0/(0.747*1.0)*x*(TMath::Cos(y)*[0]+[1]+TMath::Cos(2.0*y)*[2])";
   int nbinsx = 8;
   int nbinsy = 8;
   double xmin = -1.0;
   double xmax = 1.0;
   double ymin = 0.0;
   double ymax = 2.0*TMath::Pi();
   int nparams = 3;
   double params[3] = {0.00,0.1,0.00}; // initial parameters of fit function before fitting
   std::string fitopt = "S";

   // Create asymmetry histograms
   TH2D* h_plus  = new TH2D("h_plus","h_plus",nbinsx,xmin,xmax,nbinsy,ymin,ymax);
   TH2D* h_minus = new TH2D("h_minus","h_minus",nbinsx,xmin,xmax,nbinsy,ymin,ymax);
   data_plus.fillHistogram(h_plus,RooArgList(x,y));
   data_minus.fillHistogram(h_minus,RooArgList(x,y));
   TH2D* hasym = (TH2D*)h_plus->GetAsymmetry(h_minus);
   TCanvas *c1 = new TCanvas("c1");
   c1->cd();
   hasym->Draw("COLZ");
   
   // Set fit function
    TF2 *f1 = new TF2("f1",fitformula.c_str(),xmin,xmax,ymin,ymax);
    for (int idx=0; idx<nparams; idx++) {
        f1->SetParameter(idx,params[idx]);
        f1->SetParName(idx,Form("A%d",idx));
    }
    
    // Fit and get covariance matrix
    TFitResultPtr fr = hasym->Fit("f1",fitopt.c_str(),"S"); // IMPORTANT THAT YOU JUST FIT TO WHERE YOU STOPPED PLOTTING THE FIT VARIABLE.
    TMatrixDSym *covMat = new TMatrixDSym(fr->GetCovarianceMatrix());

    // Get fit parameters
    double * pars   = (double *)f1->GetParameters();
    double * Epars  = (double *)f1->GetParErrors();
    double  chi2    = f1->GetChisquare();
    double  ndf     = f1->GetNDF();
    double  chi2ndf = (double)chi2/ndf;

    std::cout << " params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        std::cout << pars[idx] << "±" << Epars[idx];
        if (idx<nparams-1) { std::cout << " , "; }
    }
    std::cout << "]" << std::endl;
    std::cout << " chi2/ndf = " << chi2ndf << std::endl;

    covMat->Print();//DEBUGGING

    std::cout<<" END BASIC ROOT FIT PROCEDURE"<<std::endl;//DEBUGGING
   
   // P l o t   d a t a   a n d   m o d e l s
   // ---------------------------------------

   // Construct plot frame in 'x' and 'y'
   RooPlot *xframe_2d = x.frame(Title("2D pdf fit x."));
   RooPlot *yframe_2d = y.frame(Title("2D pdf fit y."));

   // Plot data in 2D frames
   data->plotOn(xframe_2d);
   data->plotOn(yframe_2d);
 
   // Plot pdfs in 2D frames
   gen.plotOn(xframe_2d, LineColor(kRed));
   gen.plotOn(yframe_2d, LineColor(kRed));

   // F i t   m o d e l   t o   d a t a
   // -----------------------------
 
   // Fit pdf to data
   std::unique_ptr<RooFitResult> r{gen.fitTo(*data, Save(), PrintLevel(-1))};

   //--------------------------------------------------//
   // Print fit result
   r->Print("v");

   // Extract covariance and correlation matrix as TMatrixDSym
   const TMatrixDSym &cor = r->correlationMatrix();
   const TMatrixDSym &cov = r->covarianceMatrix();
 
   // Print correlation, covariance matrix
   std::cout << "correlation matrix" << endl;
   cor.Print();
   std::cout << "covariance matrix" << endl;
   cov.Print();

   // Print values of parameters (that now reflect fitted values and errors)
   a0.Print();
   a1.Print();
   a2.Print();
   
   // Draw all frames on a canvas
   std::string cname = Form("c_%d_nevents_%d",fit_idx,nevents);
   std::string ctitle = Form("");
   TCanvas *c = new TCanvas(cname.c_str(), cname.c_str(), 1600, 800); c->Divide(2,1);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   xframe_2d->GetYaxis()->SetTitleOffset(1.6);
   xframe_2d->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   yframe_2d->GetYaxis()->SetTitleOffset(1.6);
   yframe_2d->Draw();
   c->SaveAs(Form("%s.pdf",c->GetName()));

   // Extract fitted values
   double a0_2d_val = a0.getVal();
   double a1_2d_val = a1.getVal();
   double a2_2d_val = a2.getVal();

   // Extract fitted errors
   double a0_2d_err = a0.getError();
   double a1_2d_err = a1.getError();
   double a2_2d_err = a2.getError();

   std::vector<double> a;
   a.push_back(a0_2d_val);
   a.push_back(a1_2d_val);
   a.push_back(a2_2d_val);
   a.push_back(a0_2d_err);
   a.push_back(a1_2d_err);
   a.push_back(a2_2d_err);

    // TODO: Add basic ROOT Fitting results here
    for (int idx=0; idx<nparams; idx++) {
        a.push_back(pars[idx]);
    }
    for (int idx=0; idx<nparams; idx++) {
        a.push_back(Epars[idx]);
    }

   return a;

} // std::vector<double> getFitParams(int nevents) {


void ToyModelRepeatAsymmetryPdf() {

   // Loop nevents and do 1D and 2D fits returning fit parameters and errors
   std::vector<int> fit_idcs;
   std::vector<std::vector<double>> as;
   int nreps = 1;
   for (int fit_idx=0; fit_idx<nreps; fit_idx++) {
      int nevents = 10000;
      fit_idcs.push_back(fit_idx);
      std::vector<double> a = getFitParams(fit_idx,nevents);
      as.push_back(a);
    //   std::cout<<"a = [ ";
    //   for (int i=0; i<a.size(); i++) {
    //      std::cout<<a[i];
    //      if (i<a.size()-1) std::cout<<" , ";
    //      else std::cout<<" ]"<<std::endl;
    //   }
    //   std::cout<<"DONE"<<std::endl;//DEBUGGING
   } // for (int fit_idx=2; fit_idx<9; fit_idx++) {

   // Print out all info at the very end for plotting in python
   for (int idx=0; idx<fit_idcs.size(); idx++) {
      int nevents = fit_idcs[idx];
      // std::cout<<"----------------------------------------------------------------------"<<std::endl;
      // std::cout<<"nevents = "<<nevents<<std::endl;
      std::vector<double> a = as[idx];
      std::cout<<"[ ";
      for (int i=0; i<a.size(); i++) {
            std::cout<<a[i];
            if (i<a.size()-1) std::cout<<" , ";
            else std::cout<<" ],"<<std::endl;
         }
         // std::cout<<"DONE"<<std::endl;//DEBUGGING
   } // for (int idx=0; idx<fit_idcs.size(); idx++) {

   std::cout<<"[ ";
   for (int i=0; i<fit_idcs.size(); i++) {
         std::cout<<fit_idcs[i];
         if (i<fit_idcs.size()-1) std::cout<<" , ";
         else std::cout<<" ],"<<std::endl;
   }

} // void ToyModelRepeatAsymmetryPdf()
