#ifndef plotting_functions_H
#define plotting_functions_H
#include <iostream>

#include "PlotUtils/HistogramUtils.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/MnvVertErrorBand.h"
//#include "StackedHistogram.h"
#include "TCanvas.h"

namespace plotting {
const double xmin = 0.;
const double xmax = 20.e3;
const int nbins = 30;

const bool do_fractional_uncertainty = true;
const bool do_cov_area_norm = false;
const bool include_stat_error = false;
}  // namespace plotting

const std::string do_fractional_uncertainty_str =
    plotting::do_fractional_uncertainty ? std::string("Frac")
                                        : std::string("Abs");
const std::string do_cov_area_norm_str =
    plotting::do_cov_area_norm ? std::string("CovAreaNorm") : std::string("");

void PlotErrorSummary(PlotUtils::MnvH1D* hist, std::string label);
void PlotVertBand(std::string band, std::string method_str,
                  PlotUtils::MnvH1D* hist);
void PlotLatBand(std::string band, std::string method_str,
                 PlotUtils::MnvH1D* hist);
void PlotVertUniverse(std::string band, unsigned int universe,
                      std::string method_str, PlotUtils::MnvH1D* hist);
void PlotLatUniverse(std::string band, unsigned int universe,
                     std::string method_str, PlotUtils::MnvH1D* hist);
void PlotCVAndError(PlotUtils::MnvH1D* hist, std::string label);
void PlotTotalError(PlotUtils::MnvH1D* hist, std::string method_str);

void PlotTotalError(PlotUtils::MnvH1D* hist, std::string method_str) {
  TH1D* hTotalErr = (TH1D*)hist
                        ->GetTotalError(plotting::include_stat_error,
                                        plotting::do_fractional_uncertainty,
                                        plotting::do_cov_area_norm)
                        .Clone(Form("h_total_err_errSum_%d", __LINE__));
  TCanvas cF("c4", "c4");
  hTotalErr->SetTitle(
      Form("Total Uncertainty (%s); E_{#nu} (MeV)", method_str.c_str()));
  hTotalErr->Draw();
  cF.Print(Form("Enu_TotalUncertainty_%s_%s_%s.png",
                do_fractional_uncertainty_str.c_str(),
                do_cov_area_norm_str.c_str(), method_str.c_str()));
}

void PlotErrorSummary(PlotUtils::MnvH1D* hist, std::string label) {
  PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
  TCanvas cE("c1", "c1");
  hist->GetXaxis()->SetTitle(label.c_str());

  mnvPlotter.error_summary_group_map.clear();
  mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("genie_FrAbs_N");
  mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("genie_FrAbs_pi");
  mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("genie_FrCEx_N");
  mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("genie_FrCEx_pi");
  mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("genie_FrElas_N");
  mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("genie_FrElas_pi");
  mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("genie_FrInel_N");
  mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("genie_FrInel_pi");
  mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("genie_FrPiProd_N");
  mnvPlotter.error_summary_group_map["Genie_FSI"].push_back(
      "genie_FrPiProd_pi");
  mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("genie_MFP_N");
  mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("genie_MFP_pi");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_AGKYxF1pi");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_AhtBY");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_BhtBY");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_CCQEPauliSupViaKF");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_CV1uBY");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_CV2uBY");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_EtaNCEL");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_MaCCQE");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_MaCCQEshape");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_MaNCEL");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_MaRES");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_MvRES");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_NormCCQE");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_NormCCRES");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_NormDISCC");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_NormNCRES");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_RDecBR1gamma");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_Rvn1pi");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_Rvn2pi");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_Rvn3pi");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_Rvp1pi");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_Rvp2pi");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_Theta_Delta2Npi");
  mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
      "genie_VecFFCCQEshape");

  mnvPlotter.DrawErrorSummary(hist, "TR", plotting::include_stat_error, true,
                              0.0, plotting::do_cov_area_norm, "",
                              plotting::do_fractional_uncertainty);
  std::string plotname =
      Form("Enu_ErrorSummary_%s_%s_%s", do_fractional_uncertainty_str.c_str(),
           do_cov_area_norm_str.c_str(), label.c_str());
  mnvPlotter.MultiPrint(&cE, plotname, "png");
}

void PlotVertBand(std::string band, std::string method_str,
                  PlotUtils::MnvH1D* hist) {
  TH1* h1 = (TH1*)hist->GetVertErrorBand(band.c_str())
                ->GetErrorBand(plotting::do_fractional_uncertainty,
                               plotting::do_cov_area_norm)
                .Clone(Form("Enu_%s_%s", band.c_str(), method_str.c_str()));
  // TH1* h1 = (TH1*)hist->GetVertErrorBand(band.c_str());
  TCanvas cF("c4", "c4");
  h1->SetTitle(Form("%s Uncertainty (%s); E_{#nu} (MeV)", band.c_str(),
                    method_str.c_str()));
  h1->Draw("h");
  cF.Print(Form("Enu_%s_band_%s.png", band.c_str(), method_str.c_str()));
}

void PlotLatBand(std::string band, std::string method_str,
                 PlotUtils::MnvH1D* hist) {
  TH1* h1 = (TH1*)hist->GetLatErrorBand(band.c_str())
                ->GetErrorBand(plotting::do_fractional_uncertainty,
                               plotting::do_cov_area_norm)
                .Clone(Form("Enu_%s_%s", band.c_str(), method_str.c_str()));
  // TH1* h1 = (TH1*)hist->GetLatErrorBand(band.c_str());
  TCanvas cF("c1", "c1");
  h1->SetTitle(Form("%s Uncertainty (%s); E_{#nu} (MeV)", band.c_str(),
                    method_str.c_str()));
  h1->Draw("h");
  cF.Print(Form("Enu_%s_band_%s.png", band.c_str(), method_str.c_str()));
}

void PlotVertUniverse(std::string band, unsigned int universe,
                      std::string method_str, PlotUtils::MnvH1D* hist) {
  TH1D* h1 = hist->GetVertErrorBand(band.c_str())->GetHist(universe);

  TCanvas cF("c1", "c1");
  h1->SetLineColor(kBlack);
  h1->SetLineStyle(1);
  h1->Draw("hist");
  cF.Print(Form("Enu_%s_band_universe%i_%s.png", band.c_str(), universe + 1,
                method_str.c_str()));
}

void PlotLatUniverse(std::string band, unsigned int universe,
                     std::string method_str, PlotUtils::MnvH1D* hist) {
  TH1D* h1 = hist->GetLatErrorBand(band.c_str())->GetHist(universe);
  TCanvas cF("c1", "c1");
  h1->SetLineColor(kBlack);
  h1->SetLineStyle(1);
  h1->Draw("hist");
  cF.Print(Form("Enu_%s_band_universe%i_%s.png", band.c_str(), universe + 1,
                method_str.c_str()));
}

void PlotCVAndError(PlotUtils::MnvH1D* datahist,PlotUtils::MnvH1D* hist, std::string label,double mcScale) {
  PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
  TCanvas cE("c1", "c1");
  hist->GetXaxis()->SetTitle(Form("%s", label.c_str()));
//  PlotUtils::MnvH1D* datahist = new PlotUtils::MnvH1D(
//      "adsf", " ", plotting::nbins, plotting::xmin, plotting::xmax);
  bool statPlusSys = true;
 // int mcScale = 1.;
  bool useHistTitles = false;
  const PlotUtils::MnvH1D* bkgdHist = NULL;
  const PlotUtils::MnvH1D* dataBkgdHist = NULL;
  mnvPlotter.DrawDataMCWithErrorBand(datahist, hist, mcScale, "TL",
                                     useHistTitles, NULL, NULL, false,
                                     statPlusSys);
  // mnvPlotter.DrawMCWithErrorBand(hist); //I think that this call only shows
  // stat errors.
  std::string plotname = Form("CV_w_err_%s", label.c_str());
  mnvPlotter.MultiPrint(&cE, plotname, "png");
  delete datahist;
}



void PlotStacked(PlotUtils::MnvH1D* data,const TObjArray& array_mc, double mcScale, const string var,  std::string outfile_tag = "", std::string plot_title = "", double ymax = -1) {
  // Never don't clone when plotting
  TObjArray array = *(TObjArray*)array_mc.Clone("mc");

  PlotUtils::MnvPlotter mnv_plotter;
  if (ymax > 0) mnv_plotter.axis_maximum = ymax;
  std::string outfile_name = Form("StackedBreakdown_%s", outfile_tag.c_str());
  TCanvas cE("c1", "c1");
  double infoX = .68;
  double infoY = .78;
  
  //MnvH1D *h = (MnvH1D*) array.At(0); 
  //mnv_plotter.axis_maximum = (h->GetMaximum()*1.3);
  std::string x_label = var.c_str();
  //std::string x_label = ((PlotUtils::MnvH1D*)array.At(0))->GetXaxis()->GetTitle();
  std::string y_label = ((PlotUtils::MnvH1D*)array.At(0))->GetYaxis()->GetTitle();

//  double mcScale = 1.;
//  PlotUtils::MnvH1D* data = new PlotUtils::MnvH1D(
  //    "dummy", "dummy", plotting::nbins, plotting::xmin, plotting::xmax);

  // By default, this function uses the data x and y labels, instead of the mc
  // x and y labels. Given that this is much more commonly used with MC, this
  // seems like a bad choice, but who am I to break backwards compatibility?
  // Anyway, it means that we need to provide the x and y labels manually.
  mnv_plotter.DrawDataStackedMC(data, &array, mcScale, "TR", "Data", -1, -1, 1001, x_label.c_str(), y_label.c_str());
  //mnv_plotter.DrawDataStackedMC(data, &array, mcScale, "TR", "Data", 2, 1, 1001, x_label.c_str(), y_label.c_str());

  mnv_plotter.WriteNorm("Abs-Normalized", infoX, infoY);
  mnv_plotter.WritePreliminary("TL");
  mnv_plotter.AddHistoTitle(plot_title.c_str());
  mnv_plotter.MultiPrint(&cE, outfile_name, "png");
}

void Plot2D(PlotUtils::MnvH2D* hist, std::string name, std::string label_xaxis,
            std::string label_yaxis) {

    TCanvas c ("c1","c1"); 
    bool draw_as_matrix = false;
    hist->Draw("scat");
    hist->GetXaxis()->SetTitle(label_xaxis.c_str());
    //hist->GetYaxis()->SetRangeUser(0.0,5.0);
   // hist->GetXaxis()->SetRangeUser(0.0,3.0);
    hist->GetYaxis()->SetTitle(label_yaxis.c_str());
    c.Update();
    c.Print(Form("%s_2D.png", name.c_str()));
}

void PlotRatio(PlotUtils::MnvH1D* dataHist, PlotUtils::MnvH1D* mcHist, double mcScale, double dataPOT, double mcPOT, const string var,  std::string outfile_tag = "", std::string plot_title = "", double ymax = -1){
//void PlotRatio(PlotUtils::MnvH1D* dataHist, TObject hlist, double mcScale, double dataPOT, double mcPOT, const string var,  std::string outfile_tag = "", std::string plot_title = "", double ymax = -1){

  PlotUtils::MnvPlotter mnv_plotter;
  if (ymax > 0) mnv_plotter.axis_maximum = ymax;
  std::string outfile_name = Form("RatioPlot_%s", outfile_tag.c_str());
  TCanvas cE("c1", "c1");

 // hist->GetYaxis()->SetRangeUser(0.0,5.0);
 // hist->GetXaxis()->SetRangeUser(0.0,3.0);
  std::string x_label = var.c_str();
  std::string y_label = "Ratio";
  //mnv_plotter.DrawDataMCRatio(dataHist, mcHist, mcScale, true, 0.95, 1.05, y_label.c_str());
  //For plastic this range works
  mnv_plotter.DrawDataMCRatio(dataHist, mcHist, mcScale, true, true, 0.0, 2.0, x_label.c_str(), y_label.c_str());
  //mnv_plotter.DrawDataMCRatio(dataHist, hlist, mcScale, true, true, 0.0, 2.0, x_label.c_str(), y_label.c_str());
  //mnv_plotter.WriteNorm("Abs-Normalized", infoX, infoY);
  mnv_plotter.WritePreliminary("TR");
  //mnv_plotter.AddPOTNormBox(dataPOT, mcPOT, 0.3, 1.2, 0.03);
  mnv_plotter.AddChi2Label( dataHist, mcHist, 1., "TL" );
  mnv_plotter.AddHistoTitle(plot_title.c_str());
  mnv_plotter.MultiPrint(&cE, outfile_name, "png");

}


#endif

