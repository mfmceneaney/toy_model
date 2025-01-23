#include <iostream>
#include <memory>
#include <string>

// ROOT Includes
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TH1.h>
#include <ROOT/RDataFrame.hxx>
#include <Fit/Fitter.h>
#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Math/WrappedMultiTF1.h>
#include <HFitInterface.h>
#include <TGraphErrors.h>
#include <TRandom.h>
#include <TF2.h>
#include <TLatex.h>

// RooFit Includes
#include <RooCategory.h>
#include <RooRealVar.h>
#include <RooProduct.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooAbsDataHelper.h>
#include <RooDataHist.h>
#include <RooArgList.h>
#include <RooAddPdf.h>
#include <RooGenericPdf.h>
#include <RooExtendPdf.h>
#include <RooSimultaneous.h>
#include <RooFFTConvPdf.h>
#include <RooCrystalBall.h>
#include <RooLandau.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooFitResult.h>
#include <RooWorkspace.h>

// RooStats includes
#include "RooStats/SPlot.h"

using namespace RooFit;

std::vector<std::vector<std::vector<double>>> getGraphs(int fit_idx, int nevents) {

    //----------------------------------------------------------------------------------------------------//
    // Build datasets
    std::cout<<"[INFO]: Build datasets"<<std::endl;

    // Define discrete variables and set state labels
    RooCategory h("h", "helicity");
    h.defineType("minus", -1);
    h.defineType("plus", 1);

    // Declare independent fit variables
    RooRealVar m("m", "m fit variable", 1.08, 1.24);
    RooRealVar x("x", "x fit variable", -1.0, 1.0);
    RooRealVar y("y", "y fit variable", 0.0, 2.0*TMath::Pi());
    RooRealVar y_kin("y_kin", "y_kin", 0.2, 0.8);
    RooRealVar z("z", "z bin variable", 0.0, 1.0);

    // Construct signal mass variable pdf and parameters
    RooRealVar mu("mu","#mu",1.1157,0.0,2.0);
    RooRealVar s("s","#sigma",0.008,0.0,1.0);
    RooRealVar a_left("a_left","alpha_left",10.0);
    RooRealVar n_left("n_left","n_left",10.0);
    RooRealVar a("a","#alpha",1.0,0.0,3.0);
    RooRealVar n("n","n",2.0,2.0,10.0);
    RooCrystalBall pdf_gen_massvars_sg("pdf_gen_massvars_sg", "crystal_ball", m, mu, s, a_left, n_left, a, n); //NOTE: Include model name for uniqueness within bin.

    // Consruct background mass variable pdf and parameters
    RooRealVar b1("b1","b_{1}",  0.72,-10.0,10.0);
    RooRealVar b2("b2","b_{2}", -0.17,-10.0,10.0);
    RooRealVar b3("b3","b_{3}",  0.05,-10.0,10.0);
    RooRealVar b4("b4","b_{4}", -0.01,-10.0,10.0);
    RooChebychev pdf_gen_massvars_bg("pdf_gen_massvars_bg", "pdf_gen_massvars_bg", m, RooArgList(b1,b2,b3,b4));

    // Add signal and background functions for full mass variable pdf
    double sgfrac = 0.1;
    double sgYield_init = sgfrac * nevents;
    double bgYield_init = (1.0-sgfrac) * nevents;
    RooRealVar sgYield("sgYield", "fitted yield for signal", sgYield_init, 0., 2.0*nevents);
    RooRealVar bgYield("bgYield", "fitted yield for background", bgYield_init, 0., 2.0*nevents);
    RooAddPdf model("model", "pdf_gen_massvars_sg+pdf_gen_massvars_bg", RooArgList(pdf_gen_massvars_sg,pdf_gen_massvars_bg), RooArgList(sgYield,bgYield)); //NOTE: N-1 Coefficients!  Unless you want extended ML Fit

    // Build generation asymmetry pdf formulae
    std::string pdf_gen_fitvars_asymmetry_formula = "0.747*1.0*x[0]*x[1]*(x[3]*(x[6]*(1-0.5*x[6])/(1-x[6]+0.5*x[6]*x[6]))*cos(x[2])+x[4]*(x[6]*sqrt(1-x[6])/(1-x[6]+0.5*x[6]*x[6]))+x[5]*(x[6]*sqrt(1-x[6])/(1-x[6]+0.5*x[6]*x[6]))*cos(2.0*x[2]))";
    std::string pdf_gen_fitvars_formula = Form("1.0+%s",pdf_gen_fitvars_asymmetry_formula.c_str());

    // Build signal generation asymmetry pdfs in terms of fit variables and parameters
    double parLim = 0.5;
    RooRealVar a0_sg("a0_sg", "a0_sg", 0.0, -parLim, parLim);
    RooRealVar a1_sg("a1_sg", "a1_sg", 0.1, -parLim, parLim);
    RooRealVar a2_sg("a2_sg", "a2_sg", 0.0, -parLim, parLim);
    RooGenericPdf pdf_gen_fitvars_sg("pdf_gen_fitvars_sg", pdf_gen_fitvars_formula.c_str(), RooArgList(h,x,y,a0_sg,a1_sg,a2_sg,y_kin));

    // Build background generation asymmetry pdfs in terms of fit variables and parameters
    RooRealVar a0_bg("a0_bg", "a0_bg", 0.0, -parLim, parLim);
    RooRealVar a1_bg("a1_bg", "a1_bg", 0.0, -parLim, parLim);
    RooRealVar a2_bg("a2_bg", "a2_bg", 0.0, -parLim, parLim);
    RooGenericPdf pdf_gen_fitvars_bg("pdf_gen_fitvars_bg", pdf_gen_fitvars_formula.c_str(), RooArgList(h,x,y,a0_bg,a1_bg,a2_bg,y_kin));

    // Build fit asymmetry pdf formulae
    std::string pdf_fit_fitvars_asymmetry_formula = "0.747*1.0*x[0]*(x[2]*(x[5]*(1-0.5*x[5])/(1-x[5]+0.5*x[5]*x[5]))*cos(x[1])+x[3]*(x[5]*sqrt(1-x[5])/(1-x[5]+0.5*x[5]*x[5]))+x[4]*(x[5]*sqrt(1-x[5])/(1-x[5]+0.5*x[5]*x[5]))*cos(2.0*x[1]))";
    std::string pdf_fit_fitvars_pos_formula = Form("1.0+%s",pdf_fit_fitvars_asymmetry_formula.c_str());
    std::string pdf_fit_fitvars_neg_formula = Form("1.0-%s",pdf_fit_fitvars_asymmetry_formula.c_str());

    // Build fit asymmetry pdfs in terms of fit variables and parameters
    int nparams = 3;
    RooRealVar * params_fit[nparams];
    std::vector<double> params_init = {0.0, 0.1, 0.0};
    std::vector<std::string> params_fit_names;
    RooArgSet *pdf_fit_fitvars_fit_args = new RooArgSet();
    pdf_fit_fitvars_fit_args->add(x); //NOTE: ORDER IS IMPORTANT HERE!
    pdf_fit_fitvars_fit_args->add(y);
    for (int aa=0; aa<nparams; aa++) {
        params_fit_names.push_back((std::string)Form("A%d",aa));
        params_fit[aa] = new RooRealVar(params_fit_names[aa].c_str(), params_fit_names[aa].c_str(), params_init[aa], -parLim, parLim);
        pdf_fit_fitvars_fit_args->add(*params_fit[aa]);
    }
    pdf_fit_fitvars_fit_args->add(y_kin);
    RooGenericPdf pdf_fit_fitvars_fit_pos("pdf_fit_fitvars_fit_pos", pdf_fit_fitvars_pos_formula.c_str(), *pdf_fit_fitvars_fit_args);
    RooGenericPdf pdf_fit_fitvars_fit_neg("pdf_fit_fitvars_fit_neg", pdf_fit_fitvars_neg_formula.c_str(), *pdf_fit_fitvars_fit_args);
    RooSimultaneous pdf_fit_fitvars("pdf_fit_fitvars", "pdf_fit_fitvars", {{"plus", &pdf_fit_fitvars_fit_pos}, {"minus", &pdf_fit_fitvars_fit_neg}}, h);

    // Build independent variable pdf
    RooRealVar mean_gauss("mean_gauss", "mean_gauss", 0.5, 0.0, 1.0);
    RooRealVar sigma_gauss("sigma_gauss", "sigma_gauss", 0.5, 0.0, 2.0);
    RooGaussian pdf_gen_indepvars("pdf_gen_indepvars", "pdf_gen_indepvars", m, mean_gauss, sigma_gauss);

    // Generate signal datasets and merege
    std::unique_ptr<RooDataSet> data{pdf_gen_fitvars_sg.generate(RooArgSet(h,x,y,y_kin), sgYield_init)};
    std::unique_ptr<RooDataSet> data_massvar_sg{pdf_gen_massvars_sg.generate(RooArgSet(m), sgYield_init)};
    static_cast<RooDataSet&>(*data).merge(&static_cast<RooDataSet&>(*data_massvar_sg));

    // Generate background datasets and merge
    std::unique_ptr<RooDataSet> data_bg{pdf_gen_fitvars_bg.generate(RooArgSet(h,x,y,y_kin), bgYield_init)};
    std::unique_ptr<RooDataSet> data_massvar_bg{pdf_gen_massvars_bg.generate(RooArgSet(m), bgYield_init)};
    static_cast<RooDataSet&>(*data_bg).merge(&static_cast<RooDataSet&>(*data_massvar_bg));

    // Concatenate signal and background datasets
    data->append(*data_bg);

    // Generate independent variable dataset and merge
    std::unique_ptr<RooDataSet> data_indep_var{pdf_gen_indepvars.generate(RooArgSet(z), nevents)};
    static_cast<RooDataSet&>(*data).merge(&static_cast<RooDataSet&>(*data_indep_var));

    //TODO: Add option to just generate signal dataset and pass to below and not apply mass fit or sPlot

    //----------------------------------------------------------------------------------------------------//
    // Fit mass distribution
    std::cout<<"[INFO]: Fit mass distribution"<<std::endl;

    // Fit mass distribution
    std::unique_ptr<RooFitResult> fit_result_massvar{model.fitTo(*data, Save(), PrintLevel(-1))};

    // Plot invariant mass fit from RooFit
    RooPlot *mframe_1d = m.frame(Title("1D pdf fit mass_ppim."));
    data->plotOn(mframe_1d);
    model.plotOn(mframe_1d);
    model.plotOn(mframe_1d, Components(pdf_gen_massvars_sg), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(mframe_1d, Components(pdf_gen_massvars_bg), LineStyle(kDashed), LineColor(kBlue));
    TCanvas *c_massfit = new TCanvas("c_massfit");
    c_massfit->cd();
    gPad->SetLeftMargin(0.15);
    mframe_1d->GetYaxis()->SetTitleOffset(1.6);
    mframe_1d->Draw();

    // Define a range named "signal" in x from -5,5
    double sg_region_min = 1.11;
    double sg_region_max = 1.13;
    m.setRange("signal", sg_region_min, sg_region_max);
 
    // Integrate the signal PDF
    std::unique_ptr<RooAbsReal> igm_sig{pdf_gen_massvars_sg.createIntegral(m, NormSet(m), Range("signal"))};
    RooProduct i_sg_yield{"i_sg_yield", "i_sg_yield", {*igm_sig, sgYield}};
    double integral_sg_value = i_sg_yield.getVal();
    double integral_sg_value_error = i_sg_yield.getPropagatedError(*fit_result_massvar, m); // Note Need fit result saved for this to work!

    // Integrate the background PDF
    std::unique_ptr<RooAbsReal> igm_bg{pdf_gen_massvars_bg.createIntegral(m, NormSet(m), Range("signal"))};
    RooProduct i_bg_yield{"i_bg_yield", "i_bg_yield", {*igm_bg, bgYield}};
    double integral_bg_value = i_bg_yield.getVal();
    double integral_bg_value_error = i_bg_yield.getPropagatedError(*fit_result_massvar, m); // Note Need fit result saved for this to work!                                                                                                                

    // Get Full Model PDF integral
    std::unique_ptr<RooAbsReal> igm_model{model.createIntegral(m, NormSet(m), Range("signal"))};
    RooRealVar model_yield("model_yield", "model_yield",sgYield.getVal()+bgYield.getVal());
    double integral_model_value = igm_model->getVal() * model_yield.getVal();
    double integral_model_value_error = igm_model->getPropagatedError(*fit_result_massvar, m) * model_yield.getVal(); // Note Need fit result saved for this to work!                                                                                                                 

    // Sum on dataset
    std::string signal_cut = Form("%.8f<%s && %s<%.8f",(double)sg_region_min,m.GetName(),m.GetName(),(double)sg_region_max);
    double i_ds = data->sumEntries(signal_cut.c_str());
    double i_ds_err = TMath::Sqrt(i_ds);

    // Compute epsilon in a couple different ways
    double eps_pdf_int = integral_bg_value / i_ds;
    double eps_sg_pdf_int = 1.0 - integral_sg_value / i_ds ;

    // Compute epsilon errors in a couple different ways
    double eps_pdf_int_err = integral_bg_value_error / i_ds;
    double eps_sg_pdf_int_err = integral_sg_value_error / i_ds ;

    // Create histogram for mass fit variable
    int nbins_mass = 100;
    TH1D* h_mass = new TH1D("h_mass","h_mass",nbins_mass,m.getMin(),m.getMax());
    data->fillHistogram(h_mass,RooArgList(m));
    RooDataHist dh_mass("dh_mass", "dh_mass", m, Import(*h_mass));

    // Compute chi2 value
    RooFit::OwningPtr<RooAbsReal> chi2 = model.createChi2(dh_mass, Range("fullRange"),
                 Extended(true), DataError(RooAbsData::Poisson));
    int nparameters = (int) model.getParameters(RooArgSet(m))->size();
    int ndf = nbins_mass - nparams; //NOTE: ASSUME ALL BINS NONZERO
    double chi2ndf = (double) chi2->getVal()/ndf;

    // Create Legend Entries
    TString s_chi2, s_ntot, s_nbg, s_epsilon;
    s_chi2.Form("#chi^{2}/NDF = %.2f",chi2ndf);
    s_ntot.Form("N_{Tot} = %.2e #pm %.0f",i_ds,i_ds_err);
    s_nbg.Form("N_{bg} = %.2e #pm %.0f",integral_bg_value,integral_bg_value_error);
    s_epsilon.Form("#varepsilon = %.3f #pm %.3f",eps_pdf_int,eps_pdf_int_err);
    TString s_alpha, s_n, s_sigma, s_mu, s_c1;
    s_alpha.Form("#alpha = %.3f #pm %.3f",a.getVal(),a.getError());
    s_n.Form("n = %.2f #pm %.2f",n.getVal(),n.getError());
    s_sigma.Form("#sigma = %.5f #pm %.5f GeV",s.getVal(),s.getError());
    s_mu.Form("#mu = %.5f #pm %.2f GeV",mu.getVal(),mu.getError());
    s_c1.Form("C = %.5f #pm %.5f GeV",sgYield.getVal(),sgYield.getError());

    // Draw Legend
    TLegend *legend=new TLegend(0.45,0.2,0.875,0.625); //NOTE: FOR WITHOUT MC DECOMP
    legend->SetTextSize(0.04);
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, s_chi2, Form(" %g ",i_ds));
    legend->AddEntry((TObject*)0, s_alpha, Form(" %g ",i_ds));
    legend->AddEntry((TObject*)0, s_n, Form(" %g ",i_ds));
    legend->AddEntry((TObject*)0, s_sigma, Form(" %g ",i_ds));
    legend->AddEntry((TObject*)0, s_mu, Form(" %g ",i_ds));
    legend->AddEntry((TObject*)0, s_ntot, Form(" %g ",i_ds));
    legend->AddEntry((TObject*)0, s_nbg, Form(" %g ",i_ds));
    legend->AddEntry((TObject*)0, s_epsilon, Form(" %g ",i_ds));
    legend->Draw();

    // Save Canvas
    c_massfit->SaveAs(Form("%s__sigpdf_%s__fitidx_%d.pdf",c_massfit->GetName(),"cb",fit_idx));

    //----------------------------------------------------------------------------------------------------//
    // Compute sWeights
    std::cout<<"[INFO]: Compute sWeights"<<std::endl;

    // Create sweight dataset names
    std::string dataset_sg_name = Form("%s_sg",data->GetName());
    std::string dataset_bg_name = Form("%s_bg",data->GetName());

    // Run SPlot and create weighted datasets
    RooStats::SPlot sData{"sData", "SPlot Data", *data, &model, RooArgList(sgYield, bgYield)};
    auto& data1 = static_cast<RooDataSet&>(*data);
    RooDataSet data_sg_sw{dataset_sg_name.c_str(), data1.GetTitle(), &data1, *data1.get(), nullptr, "sgYield_sw"};
    RooDataSet data_bg_sw{dataset_bg_name.c_str(), data1.GetTitle(), &data1, *data1.get(), nullptr, "bgYield_sw"};

    //----------------------------------------------------------------------------------------------------//
    // Fit asymmetry
    std::cout<<"[INFO]: Fit asymmmetries"<<std::endl;

    // Initialize fit results arrays
    std::vector<std::vector<double>> asyms;
    std::vector<std::vector<double>> asym_errs;
    std::vector<double> binvars;
    std::vector<double> binvar_errs;

    // Set bins
    std::vector<double> bins = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    int nbins = bins.size()-1;

    // Loop bins
    std::cout<<"[INFO]: Looping bins"<<std::endl;
    for (int bin_idx=0; bin_idx<bins.size()-1; bin_idx++) {
        std::cout<<"[INFO]: bin_idx = "<<bin_idx<<std::endl;

        // Get bin limits
        double bin_min = bins[bin_idx];
        double bin_max = bins[bin_idx+1];

        // Apply bin cuts
        std::string bincut = Form("z>=%.2f && z<%.2f",bin_min,bin_max);
        RooDataSet *bin_ds = (RooDataSet*)data_sg_sw.reduce(bincut.c_str()); //NOTE: IMPORTANT USE THE sWEIGHTED DATASET!

        // Get bin count
        int count = bin_ds->sumEntries();

        // Get bin variable means
        double mean   = bin_ds->mean(z);
        double stddev = TMath::Sqrt(bin_ds->moment(z,2.0));
        binvars.push_back(mean);
        binvar_errs.push_back(stddev);

        // Set binvar outdir name
        std::string outdir = Form("z_bin_%d",bin_idx);

        // Fit the pdf to data
        std::unique_ptr<RooFitResult> r{pdf_fit_fitvars.fitTo(*bin_ds, RooFit::Save(), RooFit::SumW2Error(true), RooFit::PrintLevel(-1))};

        // Print fit result
        r->Print("v");

        // Extract covariance and correlation matrix as TMatrixDSym
        const TMatrixDSym &corMat = r->correlationMatrix();
        const TMatrixDSym &covMat = r->covarianceMatrix();

        // Print correlation, covariance matrix
        std::cout << "correlation matrix" << std::endl;
        corMat.Print();
        std::cout << "covariance matrix" << std::endl;
        covMat.Print();


        //TODO: Get asymmetry histograms in fit variables and plot together with the fitted asymmetry function

        // Get positive and negative helicity datasets
        RooDataSet bin_ds_pos("bin_ds_pos", "bin_ds_pos", &*bin_ds, RooArgSet(h,x,y), "h>0");
        RooDataSet bin_ds_neg("bin_ds_neg", "bin_ds_neg", &*bin_ds, RooArgSet(h,x,y), "h<0");

        // Set fit variable bins
        int nbinsx = 16;
        int nbinsy = 16;

        // // Create asymmetry histograms for x fit variable
        // TH1D* h_pos_x  = new TH1D("h_pos_x","h_pos_x",nbinsx,x.getMin(),x.getMax());
        // TH1D* h_neg_x = new TH1D("h_neg_x","h_neg_x",nbinsx,x.getMin(),x.getMax());
        // bin_ds_pos.fillHistogram(h_pos_x,RooArgList(x));
        // bin_ds_neg.fillHistogram(h_neg_x,RooArgList(x));
        // TH1D* hasym_x = (TH1D*)h_pos_x->GetAsymmetry(h_neg_x);
        // RooDataHist dh_x("dh_x", "dh_x", x, Import(*hasym_x));

        // // Create asymmetry histograms for y fit variable
        // TH1D* h_pos_y  = new TH1D("h_pos_y","h_pos_y",nbinsx,y.getMin(),y.getMax());
        // TH1D* h_neg_y = new TH1D("h_neg_y","h_neg_y",nbinsx,y.getMin(),y.getMax());
        // bin_ds_pos.fillHistogram(h_pos_y,RooArgList(y));
        // bin_ds_neg.fillHistogram(h_neg_y,RooArgList(y));
        // TH1D* hasym_y = (TH1D*)h_pos_y->GetAsymmetry(h_neg_y);
        // RooDataHist dh_y("dh_y", "dh_y", y, Import(*hasym_y));

        // Plot projection of fitted distribution in x.
        RooPlot *xframe = x.frame(RooFit::Title(Form("%s Projection, Bin: %s",x.GetTitle(),bincut.c_str())));
        bin_ds->plotOn(xframe);
        pdf_fit_fitvars.plotOn(xframe, ProjWData(h, *bin_ds), DataError(RooAbsData::SumW2));
        pdf_fit_fitvars.plotOn(xframe, Components("plus,minus"), ProjWData(h, *bin_ds), LineStyle(kDashed), DataError(RooAbsData::SumW2));

        // Draw the frame on the canvas
        std::string c1_x_name = Form("c1_%s__fitvarx_%s__fitidx_%d",outdir.c_str(),x.GetName(),fit_idx);
        TCanvas *c1_x = new TCanvas(c1_x_name.c_str(), c1_x_name.c_str());
        gPad->SetLeftMargin(0.15);
        xframe->GetYaxis()->SetTitleOffset(1.4);
        xframe->Draw();
        c1_x->Print(Form("%s.pdf",c1_x_name.c_str()));

        // Plot projection of fitted distribution in y.
        RooPlot *yframe = y.frame(RooFit::Title(Form("%s Projection, Bin: %s",y.GetTitle(),bincut.c_str())));
        bin_ds->plotOn(yframe);
        pdf_fit_fitvars.plotOn(yframe, ProjWData(h, *bin_ds), DataError(RooAbsData::SumW2));
        pdf_fit_fitvars.plotOn(yframe, Components("plus,minus"), ProjWData(h, *bin_ds), LineStyle(kDashed), DataError(RooAbsData::SumW2));

        // Draw the frame on the canvas
        std::string c1_y_name = Form("c1_%s__fitvary_%s__fitidx_%d",outdir.c_str(),y.GetName(),fit_idx);
        TCanvas *c1_y = new TCanvas(c1_y_name.c_str(), c1_y_name.c_str());
        gPad->SetLeftMargin(0.15);
        yframe->GetYaxis()->SetTitleOffset(1.4);
        yframe->Draw();
        c1_y->Print(Form("%s.pdf",c1_y_name.c_str()));

        // Get fit parameters
        std::vector<double> pars;
        std::vector<double> Epars;
        for (int aa=0; aa<nparams; aa++) {
            pars.push_back((double)params_fit[aa]->getVal());
            Epars.push_back((double)params_fit[aa]->getError());
        }
        asyms.push_back(pars);
        asym_errs.push_back(Epars);

    } // for (int bin_idx=0; bin_idx<bin_lims.size()-1; bin_idx++) {

    //----------------------------------------------------------------------------------------------------//
    // Plot asymmetries
    std::cout<<"[INFO]: Plot asymmetries"<<std::endl;

    // Loop results and plot
    std::vector<std::vector<std::vector<double>>> graphs;
    for (int idx=0; idx<nparams; idx++) {
        std::cout<<"[INFO]: Parameter idx = "<<idx<<std::endl;

        // Create arrays to contain graph data
        double xs[nbins];
        double exs[nbins];
        double ys[nbins];
        double eys[nbins];
        std::vector<std::vector<double>> graph;
        std::vector<double> graph_x;
        std::vector<double> graph_xerr;
        std::vector<double> graph_y;
        std::vector<double> graph_yerr;

        // Fill arrays with graph data
        for (int bin_idx=0; bin_idx<nbins; bin_idx++) {
            // std::cout<<"[INFO]: asyms["<<bin_idx<<"].size() = "<<asyms[bin_idx].size()<<std::endl;
            // std::cout<<"[INFO]: asym_errs["<<bin_idx<<"].size() = "<<asym_errs[bin_idx].size()<<std::endl;
            // std::cout<<"[INFO]: binvars["<<bin_idx<<"] = "<<binvars[bin_idx]<<std::endl;
            // std::cout<<"[INFO]: binvar_errs["<<bin_idx<<"] = "<<binvar_errs[bin_idx]<<std::endl;
            // std::cout<<"[INFO]: asyms["<<bin_idx<<"] = "<<asyms[bin_idx][idx]<<std::endl;
            // std::cout<<"[INFO]: asym_errs["<<bin_idx<<"] = "<<asym_errs[bin_idx][idx]<<std::endl;

            xs[bin_idx] = binvars[bin_idx];
            exs[bin_idx] = binvar_errs[bin_idx];
            ys[bin_idx] = asyms[bin_idx][idx];
            eys[bin_idx] = asym_errs[bin_idx][idx];

            graph_x.push_back(binvars[bin_idx]);
            graph_xerr.push_back(binvar_errs[bin_idx]);
            graph_y.push_back(asyms[bin_idx][idx]);
            graph_yerr.push_back(asym_errs[bin_idx][idx]);
        }
        graph.push_back(graph_x);
        graph.push_back(graph_xerr);
        graph.push_back(graph_y);
        graph.push_back(graph_yerr);
        graphs.push_back(graph);

        // Create graph of results binned in binvar
        TGraphErrors *gr = new TGraphErrors(nbins,xs,ys,exs,eys);
        gr->Write("gr");

        // Plot results graph
        TCanvas *c1 = new TCanvas();
        c1->SetBottomMargin(0.125);
        c1->cd(0);

        // Stylistic choices that aren't really necessary
        gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
        gr->SetMarkerSize(1.25);
        gr->GetXaxis()->SetTitleSize(0.05);
        gr->GetXaxis()->SetTitleOffset(0.9);
        gr->GetYaxis()->SetTitleSize(0.05);
        gr->GetYaxis()->SetTitleOffset(0.9);

        // More necessary stylistic choices
        gr->SetTitle("");
        gr->SetMarkerColor(4); // 4  blue
        gr->SetMarkerStyle(20); // 20 circle
        gr->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
        gr->GetXaxis()->SetTitle(z.GetTitle());
        std::string ytitle = Form("Asymmetry A%d",idx);
        gr->GetYaxis()->SetTitle(ytitle.c_str());
        gr->Draw("PA");

        // Add zero line
        TF1 *f2 = new TF1("f2","0",bins[0],bins[nbins]);
        f2->SetLineColor(1); // 1 black
        f2->SetLineWidth(1); // 1 thinnest
        f2->Draw("SAME");

        // Set outname and save
        TString fname;
        fname.Form("%s_%s_%s_%s_%.3f_%.3f_A%d__fitidx%d","MLFit",x.GetName(),y.GetName(),z.GetName(),bins[0],bins[nbins],idx,fit_idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

    }

    return graphs;

} //  std::vector<std::vector<std::vector<double>>> getGraphs(int fit_idx, int nevents)

void printGraphs(std::vector<std::vector<std::vector<double>>> graphs) {

    std::map<int,std::string> graph_idx_map;
    graph_idx_map[0] = "'x'   ";
    graph_idx_map[1] = "'xerr'";
    graph_idx_map[2] = "'y'   ";
    graph_idx_map[3] = "'yerr'";

    // Loop List of graphs
    std::cout<<"[ ";
    for (int graphs_idx=0; graphs_idx<graphs.size(); graphs_idx++) {
        std::vector<std::vector<double>> graph = graphs[graphs_idx];

        // Loop graph
        std::cout<<"{ ";
        for (int graph_idx=0; graph_idx<graph.size(); graph_idx++) {
            std::vector<double> data = graph[graph_idx];

            // Loop data in graph
            std::cout<<graph_idx_map[graph_idx].c_str()<<" : [ ";
            for (int i=0; i<data.size(); i++) {
                std::cout<<data[i];
                if (i<data.size()-1) std::cout<<" , ";
                else std::cout<<" ],"<<std::endl;
            }

            // Print out graph data separators
            if (graph_idx<graph.size()-1) std::cout<<"    ";
            else std::cout<<" },"<<std::endl;
        }

        // Print out graph separators
        if (graphs_idx<graphs.size()-1) std::cout<<"  ";
        else std::cout<<" ],"<<std::endl;
    }
}

void ToyModelSGBG() {

   // Loop nevents and do 1D and 2D fits returning fit parameters and errors
   std::vector<int> fit_idcs;
   std::vector<std::vector<std::vector<std::vector<double>>>> graphs_reps;
   int nreps = 1;
   for (int fit_idx=0; fit_idx<nreps; fit_idx++) {
      int nevents = 10000; //NOTE: Lambda dataset is actually ~3.7M events.
      fit_idcs.push_back(fit_idx);
      std::vector<std::vector<std::vector<double>>> gs = getGraphs(fit_idx, nevents);
      graphs_reps.push_back(gs);
   }

   // Print out all info at the very end for plotting in python
   for (int idx=0; idx<fit_idcs.size(); idx++) {
      int nevents = fit_idcs[idx];
      std::vector<std::vector<std::vector<double>>> graphs = graphs_reps[idx];
      printGraphs(graphs);
   }

   std::cout<<"fit_idcs = [ ";
   for (int i=0; i<fit_idcs.size(); i++) {
         std::cout<<fit_idcs[i];
         if (i<fit_idcs.size()-1) std::cout<<" , ";
         else std::cout<<" ],"<<std::endl;
   }

} // void ToyModelSGBG()
