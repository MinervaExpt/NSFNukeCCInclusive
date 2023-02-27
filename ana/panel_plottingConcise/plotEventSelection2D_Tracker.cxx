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
#include <iostream>
#include <string>
#include <stdlib.h>

using namespace PlotUtils;
using namespace std;

void makePlots(bool doMultipliers,bool doRatio, string indir, string outdir, int targetID, int targetZ, string plist)
{
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  gStyle->SetLabelSize(0.04);

  TString histFileName;
  histFileName = Form("%s/EventSelection_%s_t%d_z%02d_sys.root", indir.c_str(), plist.c_str(), targetID, targetZ);
 
  TFile *f1 = new TFile( histFileName,"read" );

  MnvH2D* dataMnv=(MnvH2D*)f1->Get("selected_data_reco2d_pZmu_pTmu");
  MnvH2D* mcMnv=(MnvH2D*)f1->Get("selected_mc_reco2d_pZmu_pTmu");
  MnvH2D* mcMnv_NotTracker =(MnvH2D*)f1->Get("selected_mc_reco2d_NotTracker_pZmu_pTmu");
  MnvH2D* mcMnv_WrongSign =(MnvH2D*)f1->Get("selected_mc_reco2d_WrongSign_pZmu_pTmu");
  MnvH2D* mcMnv_NC =(MnvH2D*)f1->Get("selected_mc_reco2d_NC_pZmu_pTmu");

  string trueZ;
  string mat;

  if (targetZ == 26){
    trueZ = "Iron";
    mat = "Fe";
  };

  if (targetZ == 82){
    trueZ = "Lead";
    mat = "Pb";
  };

  if (targetZ == 6){
    trueZ = "Carbon";
    mat = "C";
  };

  if (targetZ == 99){
    trueZ = "Tracker";
    mat = "CH";
  };
  

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

  TParameter<double> *pot_mc = (TParameter<double>*)f1->Get("MCPOT");
  TParameter<double> *pot_data = (TParameter<double>*)f1->Get("DataPOT");

  double MCPOT = pot_mc->GetVal();
  double DataPOT = pot_data->GetVal();
  cout << "MC POT: " << MCPOT << endl;
  cout << "Data POT: " << DataPOT << endl;

  double scale = DataPOT/MCPOT;

  if(targetID==99){
    dataMnv->Scale(1e-4, "width");
    mcMnv->Scale(scale*1e-4, "width");
    mcMnv_NotTracker->Scale(scale*1e-4, "width");
    mcMnv_WrongSign->Scale(scale*1e-4, "width");
    mcMnv_NC->Scale(scale*1e-4, "width");
  }
  else{
    dataMnv->Scale(1e-3, "width");
    mcMnv->Scale(scale*1e-3, "width");
  }

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
  TH2* mc_NotTracker=new TH2D(mcMnv_NotTracker->GetCVHistoWithStatError());
  TH2* mc_WrongSign =new TH2D(mcMnv_WrongSign->GetCVHistoWithStatError());
  TH2* mc_NC =new TH2D(mcMnv_NC->GetCVHistoWithStatError());


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

  mc_NotTracker->SetLineColor(38);
  mc_WrongSign->SetLineColor(kMagenta+2);
  mc_NC->SetLineColor(kGray);

  mc_NotTracker->SetLineWidth(2);
  mc_WrongSign->SetLineWidth(2);
  mc_NC->SetLineWidth(2);

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

  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.5);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);

  data->SetLineColor(kBlack);


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
  histAndOpts.push_back(std::make_pair(mc_NotTracker,       "graph0LX"));
  histAndOpts.push_back(std::make_pair(mc_WrongSign,       "graph0LX"));
  histAndOpts.push_back(std::make_pair(mc_NC,       "graph0LX"));
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
  else{
    if(targetID==99){
      gc->SetYTitle("Events (x10^{4}) per (GeV/c)^{2}");
    }
    else{
      gc->SetYTitle("Events per (x10^{3}) (GeV/c)^{2}");
    }  
  }

  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right

  // Adding title!
  TLegend* title = new TLegend(0.05, 0.95, 0.95, 1);
  if(targetID==99){
    title->SetHeader(Form("%s", trueZ.c_str()),"C"); // option "C" allows to center the header
  }
   else {
    title->SetHeader(Form("Target %d %s", targetID, trueZ.c_str()),"C"); 
  };
  title->SetBorderSize(0);
  title->SetFillStyle(0);
  title->SetTextSize(0.04);
  title->Draw();


  TLegend* leg=new TLegend(0.58, 0.1, 1.0, 0.32);
  leg->SetNColumns(2);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(dataStat, "MINERvA data", "lpe");
  leg->AddEntry(mc, "MINERvA Tune v1", "fl");
  leg->AddEntry(mc_NotTracker,"Not tracker","l");
  leg->AddEntry(mc_WrongSign,"Wrong sign","l");
  leg->AddEntry(mc_NC,"Neutral current","l");


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
    //gc->Print(doMultipliers ? "antinu-2d-eventsel-pt-multiplier_ratio.png" : "antinu-2d-eventsel-pt.png");
    //gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier_ratio.C" : "nu-2d-xsec-comps-pt_ratio.C");
  }
  
  else{
    //gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier.eps" : "nu-2d-xsec-comps-pt.eps");
    gc->Print(doMultipliers ? Form("%s/EventSelection2D_t%d_z%02d_%s_pt_multiplier.png", outdir.c_str(), targetID, targetZ, plist.c_str()) : Form("%s/EventSelection2D_t%d_z%02d_%s_pt.png", outdir.c_str(), targetID, targetZ, plist.c_str()));
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
  else{
    if(targetID==99){
      gc2->SetYTitle("Events (x10^{4}) per (GeV/c)^{2}");
    }
    else{
      gc2->SetYTitle("Events per (x10^{3}) (GeV/c)^{2}");
    }  
  }

  // Adding title!
  TLegend* title2 = new TLegend(0.05, 0.95, 0.95, 1);
  if(targetID==99){
    title2->SetHeader(Form("%s", trueZ.c_str()),"C"); // option "C" allows to center the header
  }
   else {
    title2->SetHeader(Form("Target %d %s", targetID, trueZ.c_str()),"C"); 
  };
  title2->SetBorderSize(0);
  title2->SetFillStyle(0);
  title2->SetTextSize(0.04);
  title2->Draw();
  
  gc2->Modified();
  if(doRatio){
    //gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier_ratio.eps" : "nu-2d-xsec-comps-pz_ratio.eps");
    //gc2->Print(doMultipliers ? "antinu-2d-eventsel-pz-multiplier_ratio.png" : "antinu-2d-eventsel-pz_ratio.png");
    //gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier_ratio.C" : "nu-2d-xsec-comps-pz_ratio.C");
  }

  else{
    //gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.eps" : "nu-2d-xsec-comps-pz.eps");
    gc2->Print(doMultipliers ? Form("%s/EventSelection2D_t%d_z%02d_%s_pz_multiplier.png", outdir.c_str(), targetID, targetZ, plist.c_str()) : Form("%s/EventSelection2D_t%d_z%02d_%s_pz.png", outdir.c_str(),  targetID, targetZ, plist.c_str()));
    //gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.C" : "nu-2d-xsec-comps-pz.C");
  }

}

int main(int argc, char* argv[])
{
  //makePlots(true,true,argv[1]); //multipliers on ratio don't make sense
  
  string indir = argv[1];
  string outdir = argv[2];
  int targetID = atoi(argv[3]);
  int targetZ = atoi(argv[4]);
  const string playlist= argv[5];

  const std::string plist(playlist);

  makePlots(true,false, indir, outdir, targetID, targetZ, plist);
  // do multipliers, do ratios

  return 0;
}
