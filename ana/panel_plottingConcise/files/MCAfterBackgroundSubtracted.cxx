//#include "myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
// #include "../util/plot/plot.h"

#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TEnv.h"
#include "TChain.h"
#include "TF2.h"
#include "Math/DistFunc.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TFile.h"
//#include "localColor.h"
#include "PlotUtils/MnvColors.h"
//#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"

using namespace PlotUtils;
using namespace std;

void makePlots(bool doMultipliers,bool doRatio,string location)
{
//  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  gStyle->SetLabelSize(0.04);

  //TFile f1(Form("%s_CV/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_pzmu.root",location.c_str()));//Final result
  TFile f1("Hists_EventSelection.root");
  TFile f2("BKGMethod1.root");
  MnvH2D* dataMnv=(MnvH2D*)f2.Get("h_background_subtracted_data");
  MnvH2D* mcMnv=(MnvH2D*)f2.Get("h_background_subtracted_mc");


  dataMnv->GetXaxis()->SetLabelSize(0.04);
  dataMnv->GetYaxis()->SetLabelSize(0.04);



  dataMnv->GetXaxis()->SetTitle("p_{||} (GeV/c)");

  double MCPOT = 1.049056685617040850944e+21;
  double DataPOT = 2.53761506093854621696e+20;
  double scale = DataPOT/MCPOT;


  dataMnv->Scale(1e-5, "width");
  mcMnv->Scale(scale*1e-5, "width");


  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data=new TH2D(dataMnv->GetCVHistoWithError(true));
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());
  TH2* mcTotalError=new TH2D(mcMnv->GetCVHistoWithError(true));


  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = MnvColors::GetColors(9);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);

  mcTotalError->SetLineColor(kRed);
  mcTotalError->SetLineWidth(2);
  mcTotalError->SetFillColor(kRed);
  mcTotalError->SetFillStyle(3002);

  data->SetLineColor(kBlack);
  data->SetLineWidth(2);
  data->SetFillColor(kBlack);
  data->SetFillStyle(3002);

  // These line and marker styles will be propagated to the 1D plots
  dataStat->SetMarkerStyle(kFullCircle);
  dataStat->SetMarkerSize(0.5);
  dataStat->SetLineColor(kBlack);
  dataStat->SetLineWidth(2);

  dataStat->SetLineColor(kBlack);


    

  // Make a list of the histograms we want to draw, along with the
  // draw options we want to use for them. You can add "graph" to the
  // draw options if you want the histogram to be converted to a graph
  // and then drawn. In that case the draw options are interpreted as
  // options to TGraphErrors::Draw().
  //
  // I don't know what happens if you put a "graph" first in the list,
  // so don't do that. Make sure the first item doesn't have "graph"
  // in its options
  std::vector<std::pair<TH2*, const char*> > histAndOpts;
  histAndOpts.push_back(std::make_pair(dataStat, "histpe1"));
  histAndOpts.push_back(std::make_pair(mcTotalError,       "graphe3"));
  histAndOpts.push_back(std::make_pair(mc,       "graph0LX"));
//  histAndOpts.push_back(std::make_pair(mc_qe,       "graph0LX"));
//  histAndOpts.push_back(std::make_pair(mc_res,       "graph0LX"));
//  histAndOpts.push_back(std::make_pair(mc_dis_dis,       "graph0LX"));
//  histAndOpts.push_back(std::make_pair(mc_dis_sis,       "graph0LX"));
//  histAndOpts.push_back(std::make_pair(mc_oth,       "graph0LX"));
//  histAndOpts.push_back(std::make_pair(mc_bkg,       "graph0LX"));
  histAndOpts.push_back(std::make_pair(data,     "histpe1"));



  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  vector<double> multipliers = GetScales(histAndOpts, true, 5,0.75);

  GridCanvas* gc=plotXAxis1DRebinPz(histAndOpts, "Muon Longitudinal Momentum (GeV/c)", "p_{t}", 4,4,800,500,doMultipliers ? &multipliers[0] : NULL);

  // Set the y range manually. Can also use gc->Remax() to guess automatically
  if(doRatio) gc->SetYLimits(0,1.99);
  else  gc->SetYLimits(0, 4.59);
  if(doRatio) gc->SetYTitle("Ratio bkgType/TotalBkg");
  else gc->SetYTitle("Events (x10^{5}) per (GeV/c)^{2}");
  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0.6, 0.05, 0.95, 0.32);
  leg->SetNColumns(2);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(dataStat, "MINERvA data", "l");
  leg->AddEntry(mc, "MINERvA Tune v1", "l");
//  leg->AddEntry(mc_qe,"QE+2p2h","l");
//  leg->AddEntry(mc_res,"Resonant","l");
//  leg->AddEntry(mc_dis_dis,"True DIS","l");
//  leg->AddEntry(mc_dis_sis,"Soft DIS","l");
//  leg->AddEntry(mc_oth,"Other CC","l");
//  leg->AddEntry(mc_bkg,"Background","l");

  TLegend* leg2 = new TLegend(0.6, 0.2525, 0.95, 0.32);
  leg2->SetNColumns(2);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(dataStat, "MINERvA data", "lpe");
  leg2->AddEntry(mc, "MINERvA Tune v1", "l");

  TLegend* leg3=new TLegend(0.6, 0.05, 0.95, 0.23);
  leg3->SetNColumns(2);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.03);
  TLatex *mytex = new TLatex();
  TLine *myline = new TLine();
  mytex->SetTextFont(42);
  mytex->SetTextSize(0.035);
  leg->Draw("SAME");
  //leg2->Draw("SAME");
  //  leg3->Draw("SAME");
  //  mytex->DrawLatex(0.67,0.23,"GENIE Components:");
  //  myline->DrawLine(0.67,0.225,0.847,0.225);
    
  if(doRatio){
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier_ratio.eps" : "nu-2d-xsec-comps-pt_ratio.eps");
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier_ratio.png" : "nu-2d-xsec-comps-pt_ratio.png");
//    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier_ratio.C" : "nu-2d-xsec-comps-pt_ratio.C");
  }
  else{
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier.eps" : "nu-2d-xsec-comps-pt.eps");
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier.png" : "nu-2d-xsec-comps-pt.png");
//    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier.C" : "nu-2d-xsec-comps-pt.C");
  }
  

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  vector<double> multipliers2 = GetScales(histAndOpts, false, 2.5,0.75);

  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  //GridCanvas* gc2=plotYAxis1D(histAndOpts, "Muon Transverse Momentum (GeV/c)", "p_{||}",4,4,800,500,doMultipliers ? &multipliers2[0] : NULL);
  GridCanvas* gc2=plotYAxis1DRebinPt(histAndOpts, "Muon Transverse Momentum (GeV/c)", "p_{||}",4,4,800,500,doMultipliers ? &multipliers2[0] : NULL);
  if(doRatio) gc2->SetYLimits(0,1.99);
  else  gc2->SetYLimits(0, 2.49);
  if(doRatio) gc2->SetYTitle("Ratio data/MINERvA Tune v1");
  else gc2->SetYTitle("Events (x10^{5}) per (GeV/c)^{2}");
  gc2->Modified();
  if(doRatio){
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier_ratio.eps" : "nu-2d-xsec-comps-pz_ratio.eps");
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier_ratio.png" : "nu-2d-xsec-comps-pz_ratio.png");
//    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier_ratio.C" : "nu-2d-xsec-comps-pz_ratio.C");
  }
  else{
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.eps" : "nu-2d-xsec-comps-pz.eps");
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.png" : "nu-2d-xsec-comps-pz.png");
//    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.C" : "nu-2d-xsec-comps-pz.C");
  }

}

int main(int argc, char* argv[])
{
  //  makePlots(true,true,argv[1]); //multipliers on ratio don't make sense
  makePlots(true,false,argv[1]);
  makePlots(false,true,argv[1]);
  makePlots(false,false,argv[1]);
  return 0;
}
