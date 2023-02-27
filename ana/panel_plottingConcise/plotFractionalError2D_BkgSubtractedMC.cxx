#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvPlotter.h"

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

#include "myPlotStyle.h"

#include "plot.h"

#include <iostream>
#include <algorithm>

using namespace std;
using namespace PlotUtils;
// =====================================================================
// axis==1: x axis. axis==2: y axis

// Take a vector of 1D histograms which all have the same binning and
// stack them together into a 2D histogram. The axis argument says
// which axis to stack them on
TH2* concatenateHists(vector<TH1*>& hists1D, int axis)
{
  //cout << "Check 1" << endl;
  assert(hists1D.size());
  //cout << "Check 2" << endl;
  //cout << "concatenateHists with " << hists1D.size() << " hists, axis=" << axis << endl;
  //cout << "Check 3" << endl;
  const int nyBins=14;
  //cout << "Check 4" << endl;
  const double yBins[nyBins+1]={0, 0.07, 0.15, 0.25, 0.33, 0.4, 0.47, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5,4.5};
  //cout << "Check 5" << endl;
  const int nxBins=16;
  //cout << "Check 6" << endl;
  const double xBins[nxBins+1]={1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20,40,60};
  //cout << "Check 7" << endl;
  TH2* ret=0;
  //cout << "Check 8" << endl;
  if(axis==1){
  //cout << "Check 9" << endl;
    ret=new TH2D(uniq(), TString::Format(";%s", hists1D[0]->GetXaxis()->GetTitle()),
                 hists1D[0]->GetXaxis()->GetNbins(), hists1D[0]->GetXaxis()->GetXbins()->GetArray(),
                 nyBins, yBins);
  //cout << "Check 10" << endl;
  }
//  cout << "Check 11" << endl;
  else{
  //cout << "Check 12" << endl;  
    ret=new TH2D(uniq(), TString::Format(";;%s", hists1D[0]->GetXaxis()->GetTitle()),
                 nxBins, xBins,
                 hists1D[0]->GetXaxis()->GetNbins(), hists1D[0]->GetXaxis()->GetXbins()->GetArray());
  //cout << "Check 13" << endl;
  }
  //cout << "Check 14" << endl;
  ret->SetLineColor(hists1D[0]->GetLineColor());
  //cout << "Check 15" << endl;
  ret->SetLineStyle(hists1D[0]->GetLineStyle());
  //cout << "Check 16" << endl;
  for(unsigned int iHist=0; iHist<hists1D.size(); ++iHist){
  //cout << "Check 17" << endl;
    for(int j=0; j<hists1D[0]->GetXaxis()->GetNbins()+1; ++j){
  //cout << "Check 18" << endl;
      int ixBin=axis==1 ? j       : iHist+1;
  //cout << "Check 19" << endl;
      int iyBin=axis==1 ? iHist+1 : j;
  //cout << "Check 20" << endl;
      //cout << "iHist: " << iHist << endl;
      //cout << "j: " << j << endl;
      //cout << "hists1D[iHist]: " << hists1D[iHist] << endl;
      //cout << "hists1D.size(): " << hists1D.size() << endl;
      //cout << "hists1D[iHist]->GetXaxis()->GetNbins(): " << hists1D[iHist]->GetXaxis()->GetNbins() << endl;
      double content=hists1D[iHist]->GetBinContent(j);
  //cout << "Check 21" << endl;
      ret->SetBinContent(ixBin, iyBin, content);
  //cout << "Check 22" << endl;
    }
  //cout << "Check 23" << endl;
  }
  //cout << "Check 24" << endl;
  //cout << "This is my concat hist address " << ret << endl;
  //cout << "Check 25" << endl;
  return ret;
}

// =====================================================================
vector<std::pair<TH2*, const char*> > getSystHistsAndOpts(MnvH2D* data, bool pt, TLegend *&leg, string group = "")
{
  MnvPlotter plotter;
  plotter.ApplyStyle(kCCQENuInclusiveStyle);

  plotter.error_summary_group_map.clear();
  plotter.error_summary_group_map["Muon Energy"].push_back("Muon_Energy_Resolution");
  plotter.error_summary_group_map["Muon Energy"].push_back("Muon_Energy_MINOS");
  plotter.error_summary_group_map["Muon Energy"].push_back("Muon_Energy_MINERvA");

  plotter.error_summary_group_map["Muon Efficiency"].push_back("MINOS_Reconstruction_Efficiency");

  plotter.error_summary_group_map["Muon Angle"].push_back("BeamAngleX");
  plotter.error_summary_group_map["Muon Angle"].push_back("BeamAngleY");

  plotter.error_summary_group_map["Flux"].push_back("Flux");

  plotter.error_summary_group_map["Interaction Model"].push_back("RPA_HighQ2");
  plotter.error_summary_group_map["Interaction Model"].push_back("RPA_LowQ2");
  plotter.error_summary_group_map["Interaction Model"].push_back("Low_Recoil_2p2h_Tune");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_AhtBY");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_BhtBY");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_CCQEPauliSupViaKF");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_CV1uBY");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_CV2uBY");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_EtaNCEL");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_MaCCQE");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_MaNCEL");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_MaRES");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_MvRES");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_NormDISCC");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_NormNCRES");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvn1pi");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvn2pi");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvp1pi");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvp2pi");
  plotter.error_summary_group_map["Interaction Model"].push_back("GENIE_VecFFCCQEshape");

  plotter.error_summary_group_map["FSI"].push_back("GENIE_AGKYxF1pi");
  plotter.error_summary_group_map["FSI"].push_back("GENIE_RDecBR1gamma");
  plotter.error_summary_group_map["FSI"].push_back("GENIE_Theta_Delta2Npi");
  plotter.error_summary_group_map["FSI"].push_back("GENIE_FrAbs_N");
  plotter.error_summary_group_map["FSI"].push_back("GENIE_FrAbs_pi");
  plotter.error_summary_group_map["FSI"].push_back("GENIE_FrCEx_N");
  plotter.error_summary_group_map["FSI"].push_back("GENIE_FrCEx_pi");
  plotter.error_summary_group_map["FSI"].push_back("GENIE_FrElas_N");
  plotter.error_summary_group_map["FSI"].push_back("GENIE_FrElas_pi");
  plotter.error_summary_group_map["FSI"].push_back("GENIE_FrInel_N");
  plotter.error_summary_group_map["FSI"].push_back("GENIE_FrPiProd_N");
  plotter.error_summary_group_map["FSI"].push_back("GENIE_FrPiProd_pi");
  plotter.error_summary_group_map["FSI"].push_back("GENIE_MFP_N");
  plotter.error_summary_group_map["FSI"].push_back("GENIE_MFP_pi");

  plotter.error_summary_group_map["Hadrons"].push_back("response_em");
  plotter.error_summary_group_map["Hadrons"].push_back("response_low_proton");
  plotter.error_summary_group_map["Hadrons"].push_back("response_mid_proton");
  plotter.error_summary_group_map["Hadrons"].push_back("response_high_proton");
  plotter.error_summary_group_map["Hadrons"].push_back("response_meson");
  plotter.error_summary_group_map["Hadrons"].push_back("response_other");

  plotter.error_summary_group_map["Hadrons"].push_back("GEANT_Neutron");
  plotter.error_summary_group_map["Hadrons"].push_back("GEANT_Proton");
  plotter.error_summary_group_map["Hadrons"].push_back("GEANT_Pion");

  plotter.error_summary_group_map["Target Mass"].push_back("Target_Mass_C");
  plotter.error_summary_group_map["Target Mass"].push_back("Target_Mass_Fe");
  plotter.error_summary_group_map["Target Mass"].push_back("Target_Mass_Pb");
  plotter.error_summary_group_map["Target Mass"].push_back("Target_Mass_CH");
  plotter.error_summary_group_map["Target Mass"].push_back("Target_Mass_H2O");

  //Change colour
  plotter.error_color_map["Target Mass"] = 15;
  plotter.error_color_map["Muon Efficiency"] = 46;
  plotter.error_color_map["Particle Response"] = 65;
  
  //plotter.ApplyStyle(kCheck);
  vector<string> vertnames = data->GetVertErrorBandNames();
  vector<string> latnames = data->GetLatErrorBandNames();
  // For each bin in the other variable, make a vector of the
  // systematic histograms
  const int nBins=pt ? 16 : 15;
  vector<vector<TH1*> > histsPT;
  histsPT.resize(nBins);
   
  TCanvas c;
  // Get MnvPlotter to plot all the histograms, and slurp them into histsPT
  for(int i=0; i<nBins; ++i){
    // First plot the histograms in the dummy canvas...
    //TCanvas c;
    MnvH1D* proj=pt ? data->ProjectionY(uniq(), i+1, i+1) : data->ProjectionX(uniq(), i+1, i+1);
    if(group==""){
      // plotter.DrawErrorSummary(proj, "TR", true, true, -1, false,"New Particle Response", true);
      plotter.DrawErrorSummary(proj, "TR", true, true,-1, false, "", true); 
      //plotter.DrawErrorSummary(proj, "TR", true, true, -1, false,"", false); //setting asfrac to false gives absolute error
      if(i==0){
        //cout << "About to get the legend" << endl;
        leg = new TLegend(*getPadLegend(&c));
        leg->SetName("MyLegend");
        leg->SaveAs("Legend.root");
        //cout << "This is my legend " << leg << "\t" << leg->GetName() << endl;
      }
    }
    else{
      if(std::count(vertnames.begin(),vertnames.end(),group)){
        TH1D err = proj->GetVertErrorBand(group)->GetErrorBand(true);
        err.DrawClone();
      }
      else{
        TH1D err = proj->GetLatErrorBand(group)->GetErrorBand(true);
        err.DrawClone();
      }
    }
    
    std::vector<TH1*> padHists=getPadHists(&c);
    auto toRemove = std::remove_if(padHists.begin(), padHists.end(), [](const auto& hist) { return std::string(hist->GetName()).find("_copy") != std::string::npos; });
    padHists.erase(toRemove, padHists.end());
    histsPT[i]=padHists;
  }
  // concatenateHists wants a vector of hists for each of the bins of
  // a given systematic. But histsPT is the other way round (the inner
  // vector loops over systematics).  So we have this fiddly loop to
  // make a transposed version of the vector-of-vector

  // It would have been easier to just pass the original
  // vector-of-vector into concatenateHists, and tell it which
  // systematic we wanted, but I've written and debugged this now, so
  // not changing it

  //  First index is systematic, second
  // index is bin
  vector<vector<TH1*> > histsPT_transpose;
  int nSyst=histsPT[0].size();
  cout << "There are " << nSyst << " systematics" << endl;
  histsPT_transpose.resize(nSyst);

  for(int iSyst=0; iSyst<nSyst; ++iSyst){
    cout << histsPT[0][iSyst]->GetName() << endl;
    for(unsigned int iBin=0; iBin<histsPT.size(); ++iBin){  
      histsPT_transpose[iSyst].push_back(histsPT[iBin][iSyst]);
    }
  }
  vector<std::pair<TH2*, const char*> > histsPT2D;
  // TODO: Figure out why the last systematic is crashing
  for(int iSyst=0; iSyst<histsPT_transpose.size(); ++iSyst){
    cout << iSyst << endl;
    TH2* h2d=concatenateHists(histsPT_transpose[iSyst], pt ? 2 : 1);
    // We want to draw all of these histograms as graphs, and exclude
    // the zero bins, to get rid of ROOT artifacts. The "graph0" draw
    // option does that (and I made it safe to pass all graphs)
    histsPT2D.push_back(std::make_pair(h2d, "graph0 l"));
  }
  //cout << "Done getSystHistsAndOpts " << endl;
  return histsPT2D;
}

// =====================================================================
void makePlots(bool pt, bool drawGroups, int areanorm,  string indir, string outdir, int targetID, int targetZ, string plist)
{
  // This turns out to be complicated (well, I could have made it less
  // bad, but this way is complicated and generalizable):
  //
  // MnvPlotter knows how to make the histograms we need, including
  // fancy stuff like grouping systematics together and sticking to
  // colour schemes. But it doesn't know about GridCanvas, and it
  // goes off and draws the histograms in its own way
  //
  // So we're going to take projections, and ask MnvPlotter to plot
  // the projection in a dummy canvas. Then we'll grab all the
  // histograms from that canvas and hold onto them.
  //
  // We could just grab all those 1D histograms and plot them straight
  // into the appropriate panel of our GridCanvas, but I want to reuse
  // the plotpz1D and plotpT1D functions in plot.h, because they know
  // how to do fanciness like squashing the tail in pz. Those
  // functions take a vector of 2D histograms, so we have to take our
  // 1D histograms from MnvPlotter, and stack them back together into
  // 2D histograms. That's what the concatenateHists function does

  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0); //commented out Aug9
  gStyle->SetEndErrorSize(2);

  TString histFileName;
  histFileName = Form("%s/BkgSubtracted_EventSelection2D_%s_t%d_z%02d_sys.root", indir.c_str(), plist.c_str(), targetID, targetZ);
 
  TFile *f = new TFile( histFileName,"read" );

  MnvH2D* data=(MnvH2D*)f->Get("h_bkg_subtracted_mc_pZmu_pTmu");

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
  

  MnvH2D* zerobin = (MnvH2D*)data->Clone("zero");
  zerobin->Reset();
  zerobin->ClearAllErrorBands();
  for(int i=0;i<zerobin->GetNbinsX();i++){
    for(int j=0;j<zerobin->GetNbinsY();j++){
      if(i==5 && j==12){
	      zerobin->SetBinContent(i+1,j+1,0.0);
	      zerobin->SetBinError(i+1,j+1,0.0);
      }
      else{
	      zerobin->SetBinContent(i+1,j+1,1.0);
	      zerobin->SetBinError(i+1,j+1,0.0);
      }
    }
  }
  
  zerobin->AddMissingErrorBandsAndFillWithCV(*data);
  //zerobin->SaveAs("test.root");
  //  data->Multiply(data,zerobin);
  if(areanorm==1) data = new MnvH2D(data->GetAreaNormalizedCopy());
  TLegend *leg = NULL;
  TLegend *title = new TLegend(0.05, 0.95, 0.95, 1);
  vector<std::pair<TH2*, const char*> > histsPT2D=getSystHistsAndOpts(data,  pt, *&leg);
  //cout << "I AM HERE" << endl;
  //cout << "Now my legend is "<< leg << endl;
  vector<string> names = data->GetErrorBandNames();
  //cout << "Start Sysnames " << endl;
  if(!drawGroups){
    myPlotStyle();
    if(pt){
      GridCanvas* gcPT=plotXAxis1DRebinPz(histsPT2D, "Muon Longitudinal Momentum (GeV/c)", "P_{t}", 4,4,800,500,NULL);
      gcPT->cd(0);
      gcPT->SetYTitle("Fractional uncertainty");
      //gcPT->SetYLimits(0, 0.29);
      gcPT->SetYLimits(0, 0.5);
      gcPT->Modified();
      leg->SetNColumns(2);
      //0.78, 0.05, 0.95, 0.32
      leg->SetX1(0.58);
      leg->SetY1(0.10);
      leg->SetX2(1.0);
      leg->SetY2(0.32);
      leg->SetTextFont(42);
      leg->Draw("SAME");

      // Adding title!
      if(targetID==99){
        title->SetHeader(Form("%s: Bkg Subtracted (MC)", trueZ.c_str()),"C"); // option "C" allows to center the header
      }
      else {
        title->SetHeader(Form("Target %d %s: Bkg Subtracted (MC)", targetID, trueZ.c_str()),"C"); 
      };
      title->SetBorderSize(0);
      title->SetFillStyle(0);
      title->SetTextSize(0.04);
      title->Draw("SAME");
      
      //gcPT->Print(Form("errors-pt_areanorm_%d.eps",areanorm));
      gcPT->Print(Form("%s/BkgSubtracted_EventSelection2D_FracErr_t%d_z%02d_%s_pt_%d_MC.png", outdir.c_str(), targetID, targetZ, plist.c_str(),areanorm));
      //gcPT->Print(Form("errors-pt_areanorm_%d.C",areanorm));
      cout << "Done" << endl;

    }
    else{
      GridCanvas* gcPT=plotYAxis1DRebinPt(histsPT2D, "Muon Transverse Momentum (GeV/c)", "P_{||}",4,4,800,500,NULL);
      // Adding title!
      TLegend* title2 = new TLegend(0.05, 0.95, 0.95, 1);
      if(targetID==99){
        title2->SetHeader(Form("%s: Bkg Subtracted (MC)", trueZ.c_str()),"C"); // option "C" allows to center the header
      }
      else {
        title2->SetHeader(Form("Target %d %s: Bkg Subtracted (MC)", targetID, trueZ.c_str()),"C"); 
      };
      title2->SetBorderSize(0);
      title2->SetFillStyle(0);
      title2->SetTextSize(0.04);
      title2->Draw();

      gcPT->SetYTitle("Fractional uncertainty");
      gcPT->SetYLimits(0, 0.29);
      gcPT->Modified();
      //gcPT->Print(Form("errors-pz_areanorm_%d.eps",areanorm));
      gcPT->Print(Form("%s/BkgSubtracted_EventSelection2D_FracErr_t%d_z%02d_%s_pz_%d_MC.png", outdir.c_str(), targetID, targetZ, plist.c_str(),areanorm));
      //gcPT->Print(Form("errors-pz_areanorm_%d.C",areanorm));
    }
    //cout << "DONE with !drawGroups" << endl;
  }
  else{
    //cout << "Other method" << endl;
    for(int n=0;n<names.size();n++){
      string group = names[n];
      vector<std::pair<TH2*, const char*> > histsPT2D=getSystHistsAndOpts(data,  pt, *&leg, group);


      TLegend *leg2 = new TLegend(0.98, 0.05, 0.95, 0.32);
      leg2->SetLineColor(kWhite);
      leg2->SetFillColor(kWhite);
      leg2->SetTextFont(42);
      for(int j=0;j<histsPT2D.size();j++)leg2->AddEntry(histsPT2D[j].first,group.c_str(),"l");
      if(pt){
        GridCanvas* gcPT=plotXAxis1DRebinPz(histsPT2D, "P_{||} (GeV/c)", "P_{t}", 4,4,800,500,NULL);
        gcPT->SetYTitle("Fractional uncertainty");
        //gcPT->SetYLimits(0, 0.25);
        gcPT->SetYLimits(0, 0.5);
        leg2->Draw("SAME");

        // Adding title!
        TLegend* title3 = new TLegend(0.05, 0.95, 0.95, 1);
        if(targetID==99){
          title3->SetHeader(Form("%s: Bkg Subtracted (MC)", trueZ.c_str()),"C"); // option "C" allows to center the header
        }
        else {
          title3->SetHeader(Form("Target %d %s: Bkg Subtracted (MC)", targetID, trueZ.c_str()),"C"); 
        };
        title3->SetBorderSize(0);
        title3->SetFillStyle(0);
        title3->SetTextSize(0.04);
        title3->Draw("SAME");

        gcPT->Modified();
        //gcPT->Print(Form("errors-%s-pt_areanorm_%d.eps",group.c_str(),areanorm));
        gcPT->Print(Form("%s/BkgSubtracted_EventSelection2D_FracErr_t%d_z%02d_%s_pt_%d_MC.png", outdir.c_str(), targetID, targetZ, plist.c_str(),areanorm));
        //cPT->Print(Form("errors-%s-pt_areanorm_%d.C",group.c_str(),areanorm));
      }
      else{
        GridCanvas* gcPT=plotYAxis1D(histsPT2D, "P_{t} (GeV/c)", "P_{||}",4,4,800,500,NULL);
        // Adding title!
        TLegend* title4 = new TLegend(0.05, 0.95, 0.95, 1);
        if(targetID==99){
          title4->SetHeader(Form("%s: Bkg Subtracted (MC)", trueZ.c_str()),"C"); // option "C" allows to center the header
        }
        else {
          title4->SetHeader(Form("Target %d %s: Bkg Subtracted (MC)", targetID, trueZ.c_str()),"C"); 
        };
        title4->SetBorderSize(0);
        title4->SetFillStyle(0);
        title4->SetTextSize(0.04);
        title4->Draw();

        gcPT->SetYTitle("Fractional uncertainty");
        gcPT->SetYLimits(0, 0.25);
        gcPT->Modified();
        //gcPT->Print(Form("errors-%s-pz_areanorm_%d.eps",group.c_str(),areanorm));
        gcPT->Print(Form("%s/BkgSubtracted_EventSelection2D_FracErr_t%d_z%02d_%s_pz_%d_MC.png", outdir.c_str(), targetID, targetZ, plist.c_str(),areanorm));
       //gcPT->Print(Form("errors-%s-pz_areanorm_%d.C",group.c_str(),areanorm));
      }
    }
  }
}

int main(int argc, char* argv[])
{
  TH1::AddDirectory(false);
  
  string indir = argv[1];
  string outdir = argv[2];
  int targetID = atoi(argv[3]);
  int targetZ = atoi(argv[4]);
  const string playlist= argv[5];

  const std::string plist(playlist);


  //for(int i=0;i<2;i++){
  makePlots(true,false, 0, indir, outdir, targetID, targetZ, plist); // not area normalized
    //makePlots(true,false,1, indir, outdir, targetID, targetZ, plist); //area normalized
    //makePlots(false,argv[1],false,i);
    //makePlots(true,argv[1],true,i);
    //makePlots(false,argv[1],true,i);
  //}
  return 0;
}
