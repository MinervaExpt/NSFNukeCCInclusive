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

void makePlots(bool doMultipliers,bool doRatio, string indir, string outdir, string plist)
{
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  gStyle->SetLabelSize(0.04);

  TString histFileName_carbon, histFileName_iron, histFileName_lead, histFileName_daisy;

  histFileName_carbon = Form("%s/CrossSection2D_Daisy_t3_z06_%s.root", indir.c_str(), plist.c_str());
  histFileName_iron = Form("%s/CrossSection2D_Daisy_t235_z26_%s.root", indir.c_str(), plist.c_str());
  histFileName_lead = Form("%s/CrossSection2D_Daisy_t2345_z82_%s.root", indir.c_str(), plist.c_str());
  histFileName_daisy = Form("%s/CrossSection2D_Daisy_t99_z99_%s.root", indir.c_str(), plist.c_str());

 
  TFile *fcarbon = new TFile( histFileName_carbon,"read" );
  TFile *firon = new TFile( histFileName_iron,"read" );
  TFile *flead = new TFile( histFileName_lead,"read" );
  TFile *fdaisy = new TFile( histFileName_daisy,"read" );

  // Target cross-sections
  MnvH2D* xsec_data_carbon = (MnvH2D*)fcarbon->Get("crossSection_data_pZmu_pTmu");
  MnvH2D* xsec_mc_carbon = (MnvH2D*)fcarbon->Get("simEventRate_crossSection_mc_pZmu_pTmu");

  MnvH2D* xsec_data_iron = (MnvH2D*)firon->Get("crossSection_data_pZmu_pTmu");
  MnvH2D* xsec_mc_iron = (MnvH2D*)firon->Get("simEventRate_crossSection_mc_pZmu_pTmu");

  MnvH2D* xsec_data_lead = (MnvH2D*)flead->Get("crossSection_data_pZmu_pTmu");
  MnvH2D* xsec_mc_lead = (MnvH2D*)flead->Get("simEventRate_crossSection_mc_pZmu_pTmu");


  // Daisy tracker cross-sections
  MnvH2D* xsec_data_tracker_carbon = (MnvH2D*)fdaisy->Get("crossSection_carbon_data_pZmu_pTmu");
  MnvH2D* xsec_mc_tracker_carbon = (MnvH2D*)fdaisy->Get("simEventRate_crossSection_carbon_mc_pZmu_pTmu");

  MnvH2D* xsec_data_tracker_iron = (MnvH2D*)fdaisy->Get("crossSection_iron_data_pZmu_pTmu");
  MnvH2D* xsec_mc_tracker_iron = (MnvH2D*)fdaisy->Get("simEventRate_crossSection_iron_mc_pZmu_pTmu");

  MnvH2D* xsec_data_tracker_lead = (MnvH2D*)fdaisy->Get("crossSection_lead_data_pZmu_pTmu");
  MnvH2D* xsec_mc_tracker_lead = (MnvH2D*)fdaisy->Get("simEventRate_crossSection_lead_mc_pZmu_pTmu");

  
  // Calculate ratios

  MnvH2D* dataMnv_carbon = xsec_data_carbon->Clone();
  dataMnv_carbon->Divide(dataMnv_carbon, xsec_data_tracker_carbon);
  MnvH2D* mcMnv_carbon= xsec_mc_carbon->Clone();
  mcMnv_carbon->Divide(mcMnv_carbon, xsec_mc_tracker_carbon );

  MnvH2D* dataMnv_iron= xsec_data_iron->Clone();
  dataMnv_iron->Divide(dataMnv_iron, xsec_data_tracker_iron);
  MnvH2D* mcMnv_iron= xsec_mc_iron->Clone();
  mcMnv_iron->Divide(mcMnv_iron, xsec_mc_tracker_iron );

  MnvH2D* dataMnv_lead = xsec_data_lead->Clone();
  dataMnv_lead->Divide(dataMnv_lead, xsec_data_tracker_lead);
  MnvH2D* mcMnv_lead = xsec_mc_lead->Clone();
  mcMnv_lead->Divide(mcMnv_lead, xsec_mc_tracker_lead );


  dataMnv_carbon->GetXaxis()->SetLabelSize(0.04);
  dataMnv_carbon->GetYaxis()->SetLabelSize(0.04);

  dataMnv_carbon->GetXaxis()->SetTitle("p_{||} (GeV/c)");

  TParameter<double> *pot_mc = (TParameter<double>*)fcarbon->Get("MCPOT");
  TParameter<double> *pot_data = (TParameter<double>*)fcarbon->Get("DataPOT");

  double MCPOT = pot_mc->GetVal();
  double DataPOT = pot_data->GetVal();
  cout << "MC POT: " << MCPOT << endl;
  cout << "Data POT: " << DataPOT << endl;

  double scale = DataPOT/MCPOT;


  //mcMnv_qe->Scale(scale*1e-5,"width");
  //mcMnv_res->Scale(scale*1e-5,"width");
  //mcMnv_dis_dis->Scale(scale*1e-5,"width");
  //mcMnv_dis_sis->Scale(scale*1e-5,"width");
  //mcMnv_2p2h->Scale(scale*1e-5,"width");
  //mcMnv_oth->Scale(scale*1e-5,"width");
  //mcMnv_bkg->Scale(scale*1e-5,"width");

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat_carbon=new TH2D(dataMnv_carbon->GetCVHistoWithStatError());
  TH2* data_carbon=new TH2D(dataMnv_carbon->GetCVHistoWithError()); // total error
  TH2* mc_carbon=new TH2D(mcMnv_carbon->GetCVHistoWithStatError());
  TH2* mcStat_carbon=new TH2D(mcMnv_carbon->GetCVHistoWithStatError());

  TH2* dataStat_iron=new TH2D(dataMnv_iron->GetCVHistoWithStatError());
  TH2* data_iron=new TH2D(dataMnv_iron->GetCVHistoWithError()); // total error
  TH2* mc_iron=new TH2D(mcMnv_iron->GetCVHistoWithStatError());
  TH2* mcStat_iron=new TH2D(mcMnv_iron->GetCVHistoWithStatError());

  TH2* dataStat_lead=new TH2D(dataMnv_lead->GetCVHistoWithStatError());
  TH2* data_lead=new TH2D(dataMnv_lead->GetCVHistoWithError()); // total error
  TH2* mc_lead=new TH2D(mcMnv_lead->GetCVHistoWithStatError());
  TH2* mcStat_lead=new TH2D(mcMnv_lead->GetCVHistoWithStatError());

  TH2* lineone=mcMnv_lead->Clone();


  // These line and marker styles will be propagated to the 1D plots
  lineone->Divide(lineone, lineone);
  lineone->SetLineColor(12);
  lineone->SetLineWidth(1.0);
  vector<int> mycolors = MnvColors::GetColors(9);
  mc_carbon->SetLineColor(46);
  mc_carbon->SetLineWidth(1.5);

  mcStat_carbon->SetLineColor(46);
  mcStat_carbon->SetLineWidth(1.5);
  mcStat_carbon->SetFillColor(kRed-10);
  mcStat_carbon->SetFillStyle(3002);

  mc_iron->SetLineColor(38);
  mc_iron->SetLineWidth(2);
  mc_iron->SetLineStyle(7);


  mcStat_iron->SetLineColor(38);
  mcStat_iron->SetLineWidth(2);
  mcStat_iron->SetFillColor(kBlue-10);
  mcStat_iron->SetFillStyle(3002);
  mcStat_iron->SetLineStyle(7);


  mc_lead->SetLineColor(30);
  mc_lead->SetLineWidth(2);
  mc_lead->SetLineStyle(3);

  mcStat_lead->SetLineColor(30);
  mcStat_lead->SetLineWidth(2);
  mcStat_lead->SetFillColor(kGreen-10);
  mcStat_lead->SetFillStyle(3002);
  mcStat_lead->SetLineStyle(7);



  // These line and marker styles will be propagated to the 1D plots
  dataStat_carbon->SetMarkerStyle(kFullCircle);
  dataStat_carbon->SetMarkerSize(0.5);
  dataStat_carbon->SetLineColor(46);
  dataStat_carbon->SetLineWidth(1);
  dataStat_carbon->SetLineColor(46);
  dataStat_carbon->SetMarkerColor(46);

  data_carbon->SetMarkerStyle(kFullCircle);
  data_carbon->SetMarkerSize(0.5);
  data_carbon->SetLineColor(46);
  data_carbon->SetLineWidth(1);
  data_carbon->SetLineColor(46);
  data_carbon->SetMarkerColor(46);

  dataStat_iron->SetMarkerStyle(22);
  dataStat_iron->SetMarkerSize(0.9);
  dataStat_iron->SetLineColor(38);
  dataStat_iron->SetLineWidth(1);
  dataStat_iron->SetLineColor(38);
  dataStat_iron->SetMarkerColor(38);
  dataStat_iron->SetLineStyle(7);


  data_iron->SetMarkerStyle(22);
  data_iron->SetMarkerSize(0.9);
  data_iron->SetLineColor(38);
  data_iron->SetLineWidth(1);
  data_iron->SetLineColor(38);
  data_iron->SetMarkerColor(38);
  data_iron->SetLineStyle(7);

  dataStat_lead->SetMarkerStyle(23);
  dataStat_lead->SetMarkerSize(0.9);
  dataStat_lead->SetLineColor(30);
  dataStat_lead->SetLineWidth(1.5);
  dataStat_lead->SetLineColor(30);
  dataStat_lead->SetMarkerColor(30);
  dataStat_lead->SetLineStyle(3);

  data_lead->SetMarkerStyle(23);
  data_lead->SetMarkerSize(0.9);
  data_lead->SetLineColor(30);
  data_lead->SetLineWidth(1.5);
  data_lead->SetLineColor(30);
  data_lead->SetMarkerColor(30);
  data_lead->SetLineStyle(3);
  
  if(doRatio){
    //TH2 *tmpden = (TH2*)mc->Clone("tmpden");
    //tmpden->Sumw2(false);
    data_carbon->Divide(mcStat_carbon);
    dataStat_carbon->Divide(mcStat_carbon);
    //mc_qe->Divide(mc);
    //mc_res->Divide(mc);
    //mc_dis_dis->Divide(mc);
    //mc_dis_sis->Divide(mc);
    //mc_oth->Divide(mc);
    //mc_bkg->Divide(mc);
    mcStat_carbon->Divide(mc_carbon);
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

  histAndOpts.push_back(std::make_pair(dataStat_carbon, "histpe1"));
  histAndOpts.push_back(std::make_pair(mcStat_carbon,       "graphe3"));
  histAndOpts.push_back(std::make_pair(mc_carbon,       "graph0LX"));
  histAndOpts.push_back(std::make_pair(data_carbon,     "histpe1"));
  histAndOpts.push_back(std::make_pair(dataStat_iron, "histpe1"));
  histAndOpts.push_back(std::make_pair(mcStat_iron,       "graphe3"));
  histAndOpts.push_back(std::make_pair(mc_iron,       "graph0LX"));
  histAndOpts.push_back(std::make_pair(data_iron,     "histpe1"));
  histAndOpts.push_back(std::make_pair(dataStat_lead, "histpe1"));
  histAndOpts.push_back(std::make_pair(mcStat_lead,       "graphe3"));
  histAndOpts.push_back(std::make_pair(mc_lead,       "graph0LX"));
  histAndOpts.push_back(std::make_pair(data_lead,     "histpe1"));
  histAndOpts.push_back(std::make_pair(lineone,     "graph0LX"));
  




  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  vector<double> multipliers = GetScales(histAndOpts, true, 5,0.75);
  GridCanvas* gc=plotXAxis1DRebinPz(histAndOpts, "Muon Longitudinal Momentum (GeV/c)", "p_{t}", 4,4,800,500,doMultipliers ? &multipliers[0] : NULL);
  // Set the y range manually. Can also use gc->Remax() to guess automatically

  if(doRatio) gc->SetYLimits(0,1.99);
  else  gc->SetYLimits(0.49, 1.49);
  //else  gc->SetYLimits(-1., 3.59);
  if(doRatio) gc->SetYTitle("Ratio bkgType/TotalBkg");
  else{
    gc->SetYTitle("d^{2}#sigma_{A}/dp_{t}dp_{||} / d^{2}#sigma_{CH}/dp_{t}dp_{||}");
  }

  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right

  // Adding title!
  TLegend* title = new TLegend(0.05, 0.95, 0.95, 1);
  title->SetHeader("Cross-Section Ratios", "C"); 
  title->SetBorderSize(0);
  title->SetFillStyle(0);
  title->SetTextSize(0.04);
  title->Draw();

  TLegend* leg=new TLegend(0.7, 0.1, 0.9, 0.32);
  //leg->SetNColumns(3);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  
  leg->AddEntry(dataStat_carbon, "Carbon", "lpe");
  //leg->AddEntry(mcStat_carbon, "Carbon Tune v430", "l");
  leg->AddEntry(dataStat_iron, "Iron", "lpe");
  //leg->AddEntry(mcStat_iron, "Iron Tune v430", "l");
  leg->AddEntry(dataStat_lead, "Lead", "lpe");
  //leg->AddEntry(mcStat_lead, "Lead Tune v430", "l");

  TLegend* leg2 = new TLegend(0.6, 0.2525, 0.95, 0.32);
  //leg2->SetNColumns(2);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(dataStat_carbon, "Carbon", "lpe");
  //leg2->AddEntry(mcStat_carbon, "Carbon Tune v430", "l");
  leg2->AddEntry(dataStat_iron, "Iron", "lpe");
  //leg2->AddEntry(mcStat_iron, "Iron Tune v430", "l");
  leg2->AddEntry(dataStat_lead, "Lead", "lpe");
  //leg2->AddEntry(mcStat_lead, "Lead Tune v430", "l");

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
    gc->Print(doMultipliers ? Form("%s/CrossSectionRatio2D_Daisy_%s_pt_multiplier.png", outdir.c_str(), plist.c_str()) : Form("%s/CrossSectionRatio2D_Daisy_%s_pt.png", outdir.c_str(), plist.c_str()));
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
      gc->SetYTitle("d^{2}#sigma_{A}/dp_{t}dp_{||} / d^{2}#sigma_{CH}/dp_{t}dp_{||}");
  }

  // Adding title!
  TLegend* title2 = new TLegend(0.05, 0.95, 0.95, 1);

  title2->SetHeader("Cross-Section Ratios", "C"); // option "C" allows to center the header
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
    //gc2->Print(doMultipliers ? Form("%s/CrossSectionRatio2D_Daisy_%s_pz_multiplier.png", outdir.c_str(), plist.c_str()) : Form("%s/CrossSectionRatio2D_Daisy_%s_pz.png", outdir.c_str(),  plist.c_str()));
    //gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.C" : "nu-2d-xsec-comps-pz.C");
  }

}

int main(int argc, char* argv[])
{
  //makePlots(true,true,argv[1]); //multipliers on ratio don't make sense
  
  string indir = argv[1];
  string outdir = argv[2];
  const string playlist= argv[3];

  const std::string plist(playlist);

  makePlots(false,false, indir, outdir, plist);
  // do multipliers, do ratios

  return 0;
}
