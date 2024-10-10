#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
using namespace RooFit;


std::vector<double> getFitParams(int nevents) {

   // S e t u p   m o d e l
   // ---------------------
 
   // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
   RooRealVar x("x", "x", 0.0, 1.0);
   RooRealVar y("y", "y", 0.0, 1.0);
   RooRealVar a0x_2d("a0x_2d", "constant term 2D x", 0.5, 0.0, 1.0);
   RooRealVar a1x_2d("a1x_2d", "linear term 2D x", 0.1, 0.0, 1.0);
   RooRealVar a0y_2d("a0y_2d", "constant term 2D y", 0.5, 0.0, 1.0);
   RooRealVar a1y_2d("a1y_2d", "linear term 2D y", 0.1, 0.0, 0.5);
   RooRealVar a0x_1d("a0x_1d", "constant term 1D x", 0.5, 0.0, 1.0);
   RooRealVar a1x_1d("a1x_1d", "linear term 1D x", 0.1, 0.0, 1.0);
   RooRealVar a0y_1d("a0y_1d", "constant term 1D y", 0.5, 0.0, 1.0);
   RooRealVar a1y_1d("a1y_1d", "linear term 1D y", 0.1, 0.0, 0.5);
   
   // Build generic pdf in terms of fit variables and parameters
   RooGenericPdf gen_x_2d("gen_x_2d", "x[1]+x[2]*x[0]", RooArgList(x,a0x_2d,a1x_2d));
   RooGenericPdf gen_y_2d("gen_y_2d", "x[1]+x[2]*x[0]", RooArgList(y,a0y_2d,a1y_2d));
   RooGenericPdf gen_x_1d("gen_x_1d", "x[1]+x[2]*x[0]", RooArgList(x,a0x_1d,a1x_1d));
   RooGenericPdf gen_y_1d("gen_y_1d", "x[1]+x[2]*x[0]", RooArgList(y,a0y_1d,a1y_1d));
   RooProdPdf gen("gen","gen",RooArgSet(gen_x_2d,gen_y_2d));

   // G e n e r a t e   e v e n t s
   // -----------------------------
 
   // Generate a dataset of 1000 events in x from gen
   std::unique_ptr<RooDataSet> data{gen.generate(RooArgSet(x,y), nevents)};
   std::unique_ptr<RooDataSet> data1Dx{gen_x_1d.generate(RooArgSet(x), nevents)};
   std::unique_ptr<RooDataSet> data1Dy{gen_y_1d.generate(RooArgSet(y), nevents)};

   // // Create binned clones
   // std::unique_ptr<RooDataHist> data_hist{data->binnedClone()};
   // std::unique_ptr<RooDataHist> data1Dx_hist{data1Dx->binnedClone()};
   // std::unique_ptr<RooDataHist> data1Dy_hist{data1Dy->binnedClone()};
   
   // P l o t   d a t a   a n d   m o d e l s
   // ---------------------------------------

   // Construct plot frame in 'x' and 'y'
   RooPlot *xframe_1d = x.frame(Title("1D pdf fit x."));
   RooPlot *yframe_1d = y.frame(Title("1D pdf fit y."));
   RooPlot *xframe_2d = x.frame(Title("2D pdf fit x."));
   RooPlot *yframe_2d = y.frame(Title("2D pdf fit y."));

   // Plot data in 1D frames
   data->plotOn(xframe_1d);
   data->plotOn(yframe_1d);

   // Plot data in 2D frames
   data->plotOn(xframe_2d);
   data->plotOn(yframe_2d);
 
   // Plot pdfs in 1D frames
   gen_x_1d.plotOn(xframe_1d, LineColor(kRed));
   gen_y_1d.plotOn(yframe_1d, LineColor(kRed));
 
   // Plot pdfs in 2D frames
   gen.plotOn(xframe_2d, LineColor(kRed));
   gen.plotOn(yframe_2d, LineColor(kRed));

   // F i t   m o d e l   t o   d a t a
   // -----------------------------
 
   // Fit pdf to data
   gen.fitTo(*data, PrintLevel(-1));
   gen_x_1d.fitTo(*data, PrintLevel(-1));
   gen_y_1d.fitTo(*data, PrintLevel(-1));
 
   // Print values of parameters (that now reflect fitted values and errors)
   a0x_2d.Print();
   a1x_2d.Print();
   a0y_2d.Print();
   a1y_2d.Print();
   a0x_1d.Print();
   a1x_1d.Print();
   a0y_1d.Print();
   a1y_1d.Print();
   
   // Draw all frames on a canvas
   std::string cname = Form("c_nevents_%d",nevents);
   std::string ctitle = Form("");
   TCanvas *c = new TCanvas(cname.c_str(), cname.c_str(), 1600, 800); c->Divide(2,2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   xframe_1d->GetYaxis()->SetTitleOffset(1.6);
   xframe_1d->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   yframe_1d->GetYaxis()->SetTitleOffset(1.6);
   yframe_1d->Draw();
   c->cd(3);
   gPad->SetLeftMargin(0.15);
   xframe_2d->GetYaxis()->SetTitleOffset(1.6);
   xframe_2d->Draw();
   c->cd(4);
   gPad->SetLeftMargin(0.15);
   yframe_2d->GetYaxis()->SetTitleOffset(1.6);
   yframe_2d->Draw();
   c->SaveAs(Form("%s.pdf",c->GetName()));

   // Extract fitted values
   double a0x_2d_val = a0x_2d.getVal();
   double a1x_2d_val = a1x_2d.getVal();
   double a0y_2d_val = a0y_2d.getVal();
   double a1y_2d_val = a1y_2d.getVal();
   double a0x_1d_val = a0x_1d.getVal();
   double a1x_1d_val = a1x_1d.getVal();
   double a0y_1d_val = a0y_1d.getVal();
   double a1y_1d_val = a1y_1d.getVal();

   // Extract fitted errors
   double a0x_2d_err = a0x_2d.getError();
   double a1x_2d_err = a1x_2d.getError();
   double a0y_2d_err = a0y_2d.getError();
   double a1y_2d_err = a1y_2d.getError();
   double a0x_1d_err = a0x_1d.getError();
   double a1x_1d_err = a1x_1d.getError();
   double a0y_1d_err = a0y_1d.getError();
   double a1y_1d_err = a1y_1d.getError();

   std::vector<double> a;
   a.push_back(a0x_2d_val);
   a.push_back(a1x_2d_val);
   a.push_back(a0y_2d_val);
   a.push_back(a1y_2d_val);
   a.push_back(a0x_1d_val);
   a.push_back(a1x_1d_val);
   a.push_back(a0y_1d_val);
   a.push_back(a1y_1d_val);
   a.push_back(a0x_2d_err);
   a.push_back(a1x_2d_err);
   a.push_back(a0y_2d_err);
   a.push_back(a1y_2d_err);
   a.push_back(a0x_1d_err);
   a.push_back(a1x_1d_err);
   a.push_back(a0y_1d_err);
   a.push_back(a1y_1d_err);

   return a;

} // std::vector<double> getFitParams(int nevents) {


void ToyModel() {

   // Loop nevents and do 1D and 2D fits returning fit parameters and errors
   std::vector<int> neventss;
   std::vector<std::vector<double>> as;
   for (int idx=4; idx<7; idx++) {
      int nevents = (int)std::pow(10,idx);
      neventss.push_back(nevents);
      std::vector<double> a = getFitParams(nevents);
      as.push_back(a);
      std::cout<<"a = [ ";
      for (int i=0; i<a.size(); i++) {
         std::cout<<a[i];
         if (i<a.size()-1) std::cout<<" , ";
         else std::cout<<" ]"<<std::endl;
      }
      std::cout<<"DONE"<<std::endl;//DEBUGGING
   } // for (int idx=2; idx<9; idx++) {

   // Print out all info at the very end for plotting in python
   for (int idx=0; idx<neventss.size(); idx++) {
      int nevents = neventss[idx];
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
   } // for (int idx=0; idx<neventss.size(); idx++) {

} // void ToyModel()
