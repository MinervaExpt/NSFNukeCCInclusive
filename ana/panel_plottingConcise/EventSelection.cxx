#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"

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
#include "TParameter.h"
#include "PlotUtils/MnvColors.h"

#include "myPlotStyle.h"

#include "plot.h"

using namespace PlotUtils;
using namespace std;

void makePlots(bool doMultipliers,bool doRatio,string location)
{
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  gStyle->SetLabelSize(0.04);

  TFile f1("/minerva/app/users/anezkak/MAT_GitHub/NSFNukeCCInclusive/ana/panel_plottingConcise/EventSelection_minervame6A_t99_z99_nosys.root");
  MnvH2D* dataMnv=(MnvH2D*)f1.Get("selected_data2d_reco_pZmu_pTmu");
  MnvH2D* mcMnv=(MnvH2D*)f1.Get("selected_mc_reco2d_pZmu_pTmu");


  dataMnv->GetXaxis()->SetLabelSize(0.04);
  dataMnv->GetYaxis()->SetLabelSize(0.04);


  //MnvH2D* mcMnv_qe = (MnvH2D*)f1.Get("reco_QE_pZmu_pTmu");//Get from N track
  //MnvH2D* mcMnv_res = (MnvH2D*)f1.Get("reco_resonant_pZmu_pTmu");//Get from N track
  //MnvH2D* mcMnv_dis_dis = (MnvH2D*)f1.Get("reco_trueDIS_pZmu_pTmu");//Get from N track
  //MnvH2D* mcMnv_dis_sis = (MnvH2D*)f1.Get("reco_softDIS_pZmu_pTmu");//Get from N track
  //MnvH2D* mcMnv_2p2h = (MnvH2D*)f1.Get("reco_2p2h_pZmu_pTmu");//Get from N track
  //MnvH2D* mcMnv_oth = (MnvH2D*)f1.Get("reco_other_pZmu_pTmu");//Get from N track
  //MnvH2D* mcMnv_bkg = (MnvH2D*)f1.Get("bkg_total_pZmu_pTmu");

  dataMnv->GetXaxis()->SetTitle("p_{||} (GeV/c)");

  TParameter<double> *pot_mc = (TParameter<double>*)f1.Get("MCPOT");
  TParameter<double> *pot_data = (TParameter<double>*)f1.Get("DataPOT");

  double MCPOT = pot_mc->GetVal();
  double DataPOT = pot_data->GetVal();
  cout << "MC POT: " << MCPOT << endl;
  cout << "Data POT: " << DataPOT << endl;

  double scale = DataPOT/MCPOT;

  dataMnv->Scale(1e-5, "width");
  mcMnv->Scale(scale*1e-5, "width");

  //mcMnv_qe->Scale(scale*1e-5,"width");
  //mcMnv_res->Scale(scale*1e-5,"width");
  //mcMnv_dis_dis->Scale(scale*1e-5,"width");
  //mcMnv_dis_sis->Scale(scale*1e-5,"width");
  //mcMnv_2p2h->Scale(scale*1e-5,"width");
  //mcMnv_oth->Scale(scale*1e-5,"width");
  //mcMnv_bkg->Scale(scale*1e-5,"width");

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());
  TH2* mcTotalError=new TH2D(mcMnv->GetCVHistoWithError());

  //TH2* mc_qe = new TH2D(mcMnv_qe->GetCVHistoWithStatError());
  //TH2* mc_res = new TH2D(mcMnv_res->GetCVHistoWithStatError());
  //TH2* mc_dis = new TH2D(mcMnv_dis->GetCVHistoWithStatError());
  //TH2* mc_dis_dis = new TH2D(mcMnv_dis_dis->GetCVHistoWithStatError());
  //TH2* mc_dis_sis = new TH2D(mcMnv_dis_sis->GetCVHistoWithStatError());
  //TH2* mc_2p2h = new TH2D(mcMnv_2p2h->GetCVHistoWithStatError());
  //TH2* mc_oth = new TH2D(mcMnv_oth->GetCVHistoWithStatError());
  //TH2* mc_bkg = new TH2D(mcMnv_bkg->GetCVHistoWithStatError());

  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = MnvColors::GetColors(9);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);

  mcTotalError->SetLineColor(kRed);
  mcTotalError->SetLineWidth(2);
  mcTotalError->SetFillColor(kRed);
  mcTotalError->SetFillStyle(3002);

  //mc_qe->SetLineColor(mycolors[2]);
  //mc_res->SetLineColor(mycolors[0]);
  //mc_dis_dis->SetLineColor(kViolet-3);
  //mc_dis_sis->SetLineColor(mycolors[1]);
  //mc_oth->SetLineColor(mycolors[4]);
  //mc_bkg->SetLineColor(mycolors[8]);

  //Add 2p2h with qe to reduce number of categories
  //mc_qe->Add(mc_2p2h);

  //mc_qe->SetLineWidth(1.5);
  //mc_res->SetLineWidth(1.5);
  //mc_dis_dis->SetLineWidth(1.5);
  //mc_dis_sis->SetLineWidth(1.5);
  //mc_oth->SetLineWidth(1.5);
  //mc_bkg->SetLineWidth(1.5);

  // These line and marker styles will be propagated to the 1D plots
  dataStat->SetMarkerStyle(kFullCircle);
  dataStat->SetMarkerSize(0.5);
  dataStat->SetLineColor(kBlack);
  dataStat->SetLineWidth(2);

  dataStat->SetLineColor(kBlack);


  if(doRatio){
    //TH2 *tmpden = (TH2*)mc->Clone("tmpden");
    //tmpden->Sumw2(false);
    data->Divide(mc);
    dataStat->Divide(mc);
    //mc_qe->Divide(mc);
    //mc_res->Divide(mc);
    //mc_dis_dis->Divide(mc);
    //mc_dis_sis->Divide(mc);
    //mc_oth->Divide(mc);
    //mc_bkg->Divide(mc);
    mc->Divide(mc);
  }
    

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
  //histAndOpts.push_back(std::make_pair(mc_qe,       "graph0LX"));
  //histAndOpts.push_back(std::make_pair(mc_res,       "graph0LX"));
  //histAndOpts.push_back(std::make_pair(mc_dis_dis,       "graph0LX"));
  //histAndOpts.push_back(std::make_pair(mc_dis_sis,       "graph0LX"));
  //histAndOpts.push_back(std::make_pair(mc_oth,       "graph0LX"));
  //histAndOpts.push_back(std::make_pair(mc_bkg,       "graph0LX"));
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
  leg->AddEntry(dataStat, "MINERvA data", "lpe");
  leg->AddEntry(mc, "MINERvA Tune v1", "l");
  //leg->AddEntry(mc_qe,"QE+2p2h","l");
  //leg->AddEntry(mc_res,"Resonant","l");
  //leg->AddEntry(mc_dis_dis,"True DIS","l");
  //leg->AddEntry(mc_dis_sis,"Soft DIS","l");
  //leg->AddEntry(mc_oth,"Other CC","l");
  //leg->AddEntry(mc_bkg,"Background","l");

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
  //leg3->AddEntry(mc_qe,"QE+2p2h","l");
  //leg3->AddEntry(mc_res,"Resonant","l");
  //leg3->AddEntry(mc_dis_dis,"True DIS","l");
  //leg3->AddEntry(mc_dis_sis,"Soft DIS","l");
  //leg3->AddEntry(mc_oth,"Other CC","l");
  //leg3->AddEntry(mc_bkg,"Background","l");
  
  TLatex *mytex = new TLatex();
  TLine *myline = new TLine();
  mytex->SetTextFont(42);
  mytex->SetTextSize(0.035);
  leg->Draw("SAME");
  //leg2->Draw("SAME");
  //leg3->Draw("SAME");
  //mytex->DrawLatex(0.67,0.23,"GENIE Components:");
  //myline->DrawLine(0.67,0.225,0.847,0.225);
    
  if(doRatio){
    //gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier_ratio.eps" : "nu-2d-xsec-comps-pt_ratio.eps");
    gc->Print(doMultipliers ? "antinu-2d-eventsel-pt-multiplier_ratio.png" : "antinu-2d-eventsel-pt.png");
    //gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier_ratio.C" : "nu-2d-xsec-comps-pt_ratio.C");
  }
  
  else{
    //gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier.eps" : "nu-2d-xsec-comps-pt.eps");
    gc->Print(doMultipliers ? "antinu-2d-eventsel-pt-multiplier.png" : "antinu-2d-eventsel-pt.png");
    //gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier.C" : "nu-2d-xsec-comps-pt.C");
  }
  

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  vector<double> multipliers2 = GetScales(histAndOpts, false, 2.5,0.75);

  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc2=plotYAxis1DRebinPt(histAndOpts, "Muon Transverse Momentum (GeV/c)", "p_{||}",4,4,800,500,doMultipliers ? &multipliers2[0] : NULL);
  if(doRatio) gc2->SetYLimits(0,1.99);
  else  gc2->SetYLimits(0, 2.49);
  if(doRatio) gc2->SetYTitle("Ratio data/MINERvA Tune v1");
  else gc2->SetYTitle("Events (x10^{5}) per (GeV/c)^{2}");
  gc2->Modified();
  if(doRatio){
    //gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier_ratio.eps" : "nu-2d-xsec-comps-pz_ratio.eps");
    gc2->Print(doMultipliers ? "antinu-2d-eventsel-pz-multiplier_ratio.png" : "antinu-2d-eventsel-pz_ratio.png");
    //gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier_ratio.C" : "nu-2d-xsec-comps-pz_ratio.C");
  }

  else{
    //gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.eps" : "nu-2d-xsec-comps-pz.eps");
    gc2->Print(doMultipliers ? "antinu-2d-eventsel-pz-multiplier.png" : "antinu-2d-eventsel-pz.png");
    //gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.C" : "nu-2d-xsec-comps-pz.C");
  }

}

int main(int argc, char* argv[])
{
  //makePlots(true,true,argv[1]); //multipliers on ratio don't make sense
  makePlots(true,false,argv[1]);
  makePlots(false,true,argv[1]);
  makePlots(false,false,argv[1]);
  return 0;
}
