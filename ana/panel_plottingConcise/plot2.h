#ifndef PLOT_H
#define PLOT_H

#include "TClass.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TGraphErrors.h"

#include <iostream>
#include <vector>
#include <math.h>
#include <iostream>

using namespace std;

//======================================================================

// Return a unique string on every call, so we can name root hists
// without it clobbering them
TString uniq()
{
  static int i=0;
  return TString::Format("uniq%d", i++);
}


//======================================================================

// Convert a histogram to a graph so we can draw it without ROOT
// drawing extra lines at the end and other such nonsense
TGraphErrors* histToGraph(TH1* h, double multiplier=1, bool includeZeros=true)
{
  TGraphErrors* grE=new TGraphErrors;
  grE->SetLineColor(h->GetLineColor());
  grE->SetLineStyle(h->GetLineStyle());
  grE->SetLineWidth(h->GetLineWidth());

  for(int i=0; i<h->GetNbinsX(); ++i){
    int bin=i+1;
    double binVal=h->GetBinContent(bin);
    if(includeZeros || binVal!=0){
      int binNo=grE->GetN();
      grE->SetPoint(binNo, h->GetBinCenter(bin), multiplier*binVal);
      grE->SetPointError(binNo, 0, multiplier*h->GetBinError(bin));
    }
  }
  return grE;
}

//======================================================================

TH1* getPadAxisHist(TPad* pad)
{
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  while (( obj=next() )) {
    // cout << obj->GetName() << endl;
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      // cout << "getPadAxisHist returning hist with name " << obj->GetName() << endl;
      return (TH1*)obj;
    }
  }
  return NULL;
}

//======================================================================
std::vector<TH1*> getPadHists(TPad* pad)
{
  std::vector<TH1*> ret;
  std::vector<string> names;
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  while (( obj=next() )) {
    //    cout << obj->GetName() << endl;
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      string name = obj->GetName();
      //      cout << "getPadAxisHist returning hist with name " << obj->GetName() <<"t"<< obj->GetTitle() << endl;
      bool isFound = false;
      for(unsigned int i=0;i<names.size();i++){
	if(names[i]==name){
	  isFound = true;
	  break;
	}
      }
      if(!isFound){
	ret.push_back((TH1*)obj);
	names.push_back(name);
      }

    }
  }
  return ret;
}

//======================================================================
double getPadMax(TPad* pad)
{
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  Double_t runningMax=-9e99;//Hparam.ymax;
  while (( obj=next() )) {
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      TH1* curHist=(TH1*)obj;
      const double thisMax=curHist->GetBinContent(curHist->GetMaximumBin());
      if (thisMax > runningMax) {
        runningMax=thisMax;
      }
    }
  }
  return runningMax;
}

//======================================================================
void reMax(TPad* pad, double ymin=0, double headroom=1.1)
{
  TH1* firstHist=getPadAxisHist(pad);
  if(!firstHist){
    std::cerr << "reMax: No histogram in pad " << pad->GetName() << ". Can't reMax" << std::endl;
    return;
  }
  firstHist->GetYaxis()->SetRangeUser(ymin, headroom*getPadMax(pad));
  pad->Update();
}

//======================================================================
void drawBinRange(TH2* h, int axis, int bin, const char* varName, const char* numFormatStr=".2f", bool left=false)
{
  double varmin=axis==1 ? h->GetXaxis()->GetBinLowEdge(bin) : h->GetYaxis()->GetBinLowEdge(bin);
  double varmax=axis==1 ? h->GetXaxis()->GetBinUpEdge(bin) :  h->GetYaxis()->GetBinUpEdge(bin);
  
  TString formatStr(TString::Format("%%%s < %%s < %%%s", numFormatStr, numFormatStr));

  TLatex* la=0;
  TString text(TString::Format(formatStr.Data(), varmin, varName, varmax));
  if(left){
    la=new TLatex(gPad->GetLeftMargin()+0.02,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(13); // top left
  }
  else{
    la=new TLatex(1-gPad->GetRightMargin()-0.01,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(33); // top right
  }

  la->SetNDC();
  la->SetTextFont(42);
  la->SetTextSize(0.03);
  la->Draw();
}
//======================================================================

void drawBinRangeSpecial(TH2* h, int axis, int bin, const char* varName, const char* numFormatStr=".2f", bool left=false, bool smallrange=false)
{
  double varmin=axis==1 ? h->GetXaxis()->GetBinLowEdge(bin) : h->GetYaxis()->GetBinLowEdge(bin);
  double varmax=axis==1 ? h->GetXaxis()->GetBinUpEdge(bin) :  h->GetYaxis()->GetBinUpEdge(bin);
  
  TString formatStr;
  if(bin==3 && axis==1) formatStr.Form("%%%s  < %%s < %%%s", numFormatStr, numFormatStr);
  else if(bin==10 && axis==2) formatStr.Form("   %%%s< %%s< %%%s", numFormatStr, numFormatStr);
  else if(bin==1 && axis==2 && smallrange)  formatStr.Form("%%%s< %%s<  %%%s", numFormatStr, numFormatStr);
  else if(bin==11 && axis==2) formatStr.Form("     %%%s< %%s<%%%s", numFormatStr, numFormatStr);
  else if(bin==12 && axis==2 && smallrange) formatStr.Form("%%%s    <%%s<%%%s", numFormatStr, numFormatStr);
  else if(bin==12 && axis==2) formatStr.Form("%%%s  < %%s< %%%s", numFormatStr, numFormatStr);
  else if(bin==13 && axis==2 && smallrange) formatStr.Form("%%%s<    %%s< %%%s", numFormatStr, numFormatStr);
  else if(bin==13 && axis==2) formatStr.Form("%%%s<  %%s < %%%s", numFormatStr, numFormatStr);
  else formatStr.Form("%%%s < %%s < %%%s", numFormatStr, numFormatStr);

  TLatex* la=0;
  TString text(TString::Format(formatStr.Data(), varmin, varName, varmax));
  double leftmargin= axis==1? 0.023:0.01;
  if(left){
    la=new TLatex(gPad->GetLeftMargin()+leftmargin,
		  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(13); // top left

  }
  else{
    la=new TLatex(1-gPad->GetRightMargin()-0.01,
		  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(33); // top right

  }


  la->SetNDC();
  la->SetTextFont(42);
  la->SetTextSize(0.03);
  la->Draw();
}




//======================================================================
void drawVarBinRange(TH2* h, int axis, int minbin, int maxbin, const char* varName, const char* numFormatStr=".2f", bool left=false)
{
  double varmin=axis==1 ? h->GetXaxis()->GetBinLowEdge(minbin) : h->GetYaxis()->GetBinLowEdge(minbin);
  double varmax=axis==1 ? h->GetXaxis()->GetBinUpEdge(maxbin) :  h->GetYaxis()->GetBinUpEdge(maxbin);
  
  TString formatStr(TString::Format("%%%s < %%s < %%%s", numFormatStr, numFormatStr));

  TLatex* la=0;
  TString text(TString::Format(formatStr.Data(), varmin, varName, varmax));
  if(left){
    la=new TLatex(gPad->GetLeftMargin()+0.02,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(13); // top left
  }
  else{
    la=new TLatex(1-gPad->GetRightMargin()-0.01,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(33); // top right
  }

  la->SetNDC();
  la->SetTextFont(42);
  la->SetTextSize(0.03);
  la->Draw();
}

//======================================================================
GridCanvas* plotpT1D(std::vector<std::pair<TH2*, const char*> > histAndOpts,
                     double* multipliers=NULL)
{
  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();


  int grid_x = sqrt(nbins_pz)+1;
  int grid_y = nbins_pz/(grid_x-1);
  //int grid_y = nbins_pt/(grid_x-1)-1;
  cout << "Plotting plotpT1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pz; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionY(uniq(), i+1, i+1);
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }

        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 1, i+1, "Emu", ".1f", false);
    //drawBinRange(histAndOpts[0].first, 1, i+1, "y", ".1f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Hadronic Energy (GeV)");
  //gc->SetXTitle("x");
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotYAxis1D(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle,
			double* multipliers=NULL)
{
  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();


  int grid_x = sqrt(nbins_pz)+1;
  int grid_y = nbins_pz/(grid_x-1);
  if(grid_x*grid_y-nbins_pz==grid_x) grid_y-=1;
  cout << nbins_pz - grid_x*grid_y << endl;
  cout << "Plotting plotYAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pz; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionY(uniq(), i+1, i+1);
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }

        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 1, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle(xaxistitle.c_str());
  //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}
///////////////For ratios/////////////////////////////////////
//======================================================================
GridCanvas* plotYfrac1D(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle,
			double* multipliers=NULL)
{
  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();


  int grid_x = sqrt(nbins_pz)+1;
  int grid_y = nbins_pz/(grid_x-1);
  if(grid_x*grid_y-nbins_pz==grid_x) grid_y-=1;
  cout << nbins_pz - grid_x*grid_y << endl;
  cout << "Plotting plotYAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pz; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    TH2* hmc=histAndOpts[1/*histAndOpts.size() -2*/].first;
    TH1* mcproj=hmc->ProjectionY(uniq(), i+1, i+1);
    


    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionY(uniq(), i+1, i+1);
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);

      proj->Divide(mcproj);
      

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }

        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 1, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle(xaxistitle.c_str());
  //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
//======================================================================
GridCanvas* plotYAxis1D(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, int ncolumns, int nrows, int xwidth, int ywidth,
			double* multipliers=NULL)
{
  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();


  GridCanvas* gc=new GridCanvas(uniq(), ncolumns, nrows, xwidth, ywidth);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pz; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionY(uniq(), i+1, i+1);
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }

        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 1, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle(xaxistitle.c_str());
  //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
GridCanvas* plotpT1DAntiNu(std::vector<std::pair<TH2*, const char*> > histAndOpts,
                     double* multipliers=NULL)
{
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), 4, 3, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<11; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionY(uniq(), i+1, i+1);
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }

        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 1, i+1, "p_{||}/GeV", ".1f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Muon transverse momentum (GeV)");
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/nucleon)");
  gc->SetTitleSize(22.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
TH1* rebinpz(TH1* h)
{
  const int nBins=16;
  // Less aggressive version
  // double newBins[nBins+1]={1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 8, 10, 12, 14};
  //1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20,40,60
  double newBins[nBins+1]={1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8,9,10,13,16,20,24};

  TH1* ret=new TH1D(uniq(), TString::Format("%s;%s;%s", h->GetTitle(), h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle()),
                    nBins, newBins);

  ret->SetLineColor(h->GetLineColor());
  ret->SetLineStyle(h->GetLineStyle());
  ret->SetLineWidth(h->GetLineWidth());
  ret->SetMarkerStyle(h->GetMarkerStyle());
  ret->SetMarkerSize(h->GetMarkerSize());
  
  // Do under/overflow too
  for(int i=0; i<nBins+2; ++i){
    ret->SetBinContent(i, h->GetBinContent(i));
    ret->SetBinError(i, h->GetBinError(i));
  }
  return ret;
}

//======================================================================
GridCanvas* plotpz1D(std::vector<std::pair<TH2*, const char*> > histAndOpts,
                     double* multipliers=NULL, bool is3D=false)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x-1)-1;
  cout << "Plotting plotpZ1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* projtmp=h2d->ProjectionX(uniq(), i+1, i+1);
      TH1* proj=is3D?(TH1*)projtmp->Clone("clone"):rebinpz(projtmp);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, "Ehad", ".2f", false);
    //drawBinRange(histAndOpts[0].first, 2, i+1, "x", ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  //const int nLabels=3;
  //const double positions[nLabels]={5, 8, 10};
  //const char* valueStrings[nLabels]={"5", "10", "20"};
  
  if(!is3D){
    //gc->SetManualXLabels(nLabels, positions, valueStrings);
    gc->SetXTitle("Muon Energy (GeV)");
    gc->SetXTitle("y");
  }
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}
//======================================================================

GridCanvas* plotfracpz(std::vector<std::pair<TH2*, const char*> > histAndOpts,
		       double* multipliers=NULL, bool doGENIEs=false, bool doPrelimLabel=false, bool smallrange=false)
{
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), 5, 3, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<13; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    
    TH2* hmc=histAndOpts[1/*histAndOpts.size() -2*/].first;
    TH1* mctmp=hmc->ProjectionX(uniq(), i+1, i+1);
    TH1* mcproj=rebinpz(mctmp);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      h2d->SetBinContent(4,12,0);
      h2d->SetBinError(4,12,0);
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* projtmp=h2d->ProjectionX(uniq(), i+1, i+1);
      TH1* proj=rebinpz(projtmp);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);
      
      proj->Divide(mcproj);
      


      if(multipliers) proj->Scale(multipliers[i]);
      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
	gStyle->SetLabelSize(0.045);
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }
    if (i==0 && !smallrange) drawBinRangeSpecial(histAndOpts[0].first, 2, i+1, " p_{T}/GeV ", ".2f", true);//false);
    else if (i==0 && smallrange)  drawBinRangeSpecial(histAndOpts[0].first, 2, i+1, "p_{T}/GeV", ".2f", true, true);//false);
    //else if (i==9)     drawBinRangeSpecial(histAndOpts[0].first, 2, i+1, "p_{T}/GeV", ".2f", true);
    else if (i==10 && smallrange)     drawBinRangeSpecial(histAndOpts[0].first, 2, i+1, "p_{T}/GeV", ".2f", true);
    else if (i==11 && smallrange)     drawBinRangeSpecial(histAndOpts[0].first, 2, i+1, "p_{T}/GeV", ".2f", true, smallrange);
    else if (i==12)     drawBinRangeSpecial(histAndOpts[0].first, 2, i+1, " p_{T}/GeV", ".2f", true, smallrange);
    else drawBinRange(histAndOpts[0].first, 2, i+1, "p_{T}/GeV", ".2f", false);

   // drawBinRange(histAndOpts[0].first, 1, i+1, "Ehad (GeV)", ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  const int nLabels=3;
  const double positions[nLabels]={5, 8, 10};
  const char* valueStrings[nLabels]={"5", "10", "20"};
  double yoffset=0.015;
  gc->SetManualXLabels(nLabels, positions, valueStrings,yoffset);
  gc->SetXTitle("Muon longitudinal momentum (GeV)");
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(22.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  if (doPrelimLabel){
    TLatex* labl=new TLatex(0.15, 0.96,"MINER#nuA Preliminary      POT: 3.34 x10^{20}      #LTE_{#nu}#GT~3.5GeV" );
    labl->SetTextSize(0.04);
    labl->SetTextFont(112);//22);                                                                            
    labl->SetTextColor(4);
    labl->Draw();
  }


  return gc;
}
//======================================================================
GridCanvas* plotfracpt(std::vector<std::pair<TH2*, const char*> > histAndOpts,
		       double* multipliers=NULL, bool doGENIEs=false, bool doPrelimLabel=false, bool smallrange=false)
{
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), 4, 3, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<12; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    //TH2* hdatastat=histAndOpts[0].first;
    //TH1* dsproj=hdatastat->ProjectionY(uniq(), i+1, i+1);
    //double area=dsproj->Integral();

    TH2* hmc=histAndOpts[1].first;//[histAndOpts.size() - 2].first;
    TH1* mcproj=hmc->ProjectionY(uniq(), i+1, i+1);
    //if (doAreaNorm) mcproj->Scale(area/mcproj->Integral());

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      h2d->SetBinContent(4,12,0);
      h2d->SetBinError(4,12,0);
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionY(uniq(), i+1, i+1);
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);
      //if (doAreaNorm) proj->Scale(area/proj->Integral());
      proj->Divide(mcproj);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }

        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    //drawBinRange(histAndOpts[0].first, 2, i+1, "Emu (GeV)", ".2f", false);
    //bool abovepanel=false;
    //if (i<4) abovepanel=true;
    if (i==0 && smallrange) drawBinRange(histAndOpts[0].first, 1, i+1, "  p_{||}/GeV", ".1f", true);
    //else if (i==1) drawBinRange(histAndOpts[0].first, 1, i+1, "p_{||}  /GeV", ".1f", true);
    else if (i==2 && smallrange) drawBinRangeSpecial(histAndOpts[0].first, 1, i+1, " p_{||}/   GeV", ".1f", true);
    //else if (i==3 && smallrange) drawBinRange(histAndOpts[0].first, 1, i+1, " p_{||} /  GeV", ".1f", true);//false);
    else if (i==4 && smallrange) drawBinRange(histAndOpts[0].first, 1, i+1, "p_{||}/GeV  ", ".1f", true);//false);
    else if (i==5 && smallrange) drawBinRange(histAndOpts[0].first, 1, i+1, "p_{||}/GeV ", ".1f", true);//false);
    else drawBinRange(histAndOpts[0].first, 1, i+1, "p_{||}/GeV", ".1f", true);//false);
    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Hadronic Energy (GeV)");
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(22.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw(); 

  /*  TLegend* leg_bin4=new TLegend(0.885, 0.684, 0.895, 0.695);//0.871, 0.683, 0.98, 0.88);                      
  leg_bin4->SetFillColor(kWhite);
  leg_bin4->SetBorderSize(0);
  leg_bin4->Draw("SAME");

  TLegend* leg_bin4b=new TLegend(0.885, 0.67, 0.895, 0.677);//0.871, 0.683, 0.98, 0.88);                      
  leg_bin4b->SetFillColor(kWhite);
  leg_bin4b->SetBorderSize(0);
  //leg_bin4b->Draw("SAME");
  */
  if (doPrelimLabel){
    TLatex* labl=new TLatex(0.15, 0.96,"MINER#nuA Preliminary      POT: 3.34 x10^{20}      #LTE_{#nu}#GT~3.5GeV" );//(0.4, 0.96,"MINER#nuA Preliminary         POT: 3.34 x10^{20}" );
    labl->SetTextSize(0.04);
    labl->SetTextFont(112);//22);                                                                            
    labl->SetTextColor(4);
    labl->Draw();
  }


  return gc;
}
////////////////////////////For Ratio Plots//////////////////////////////////

//======================================================================
GridCanvas* plotXfrac1D(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle,
                     double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x-1);
  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  cout << "Plotting plotXfrac1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    TH2* hmc=histAndOpts[1/*histAndOpts.size() -2*/].first;
    TH1* mcproj=hmc->ProjectionX(uniq(), i+1, i+1);
   // TH1* mcproj=rebinpz(mctmp);
    


    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->SetFillColor(h2d->GetFillColor());
      proj->SetFillStyle(h2d->GetFillStyle());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);
      

      proj->Divide(mcproj);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw("E2");

  return gc;
}


//
//


//======================================================================
GridCanvas* plotXAxis1D(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle,
                     double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x-1);
  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  cout << "Plotting plotXAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->SetFillColor(h2d->GetFillColor());
      proj->SetFillStyle(h2d->GetFillStyle());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw("E2");

  return gc;
}


//======================================================================
GridCanvas* plotXAxis1D_ReducedXRange(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, double xminval, double xmaxval,
                     double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x-1);
  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  cout << "Plotting plotXAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800,500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);
      proj->SetFillStyle(h2d->GetFillStyle());

      proj->GetXaxis()->SetRangeUser(xminval,xmaxval);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
    GridCanvas* plotXAxis1D_ReducedXRange(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, int grid_x, int grid_y, double xminval, double xmaxval, double xpix, double ypix,
                     double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  //  int grid_x = sqrt(nbins_pt)+1;
  //  int grid_y = nbins_pt/(grid_x-1);
  //  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  cout << "Plotting plotXAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, xpix,ypix);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);
      proj->SetFillStyle(h2d->GetFillStyle());

      proj->GetXaxis()->SetRangeUser(xminval,xmaxval);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
GridCanvas* plotXAxis1D(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, int ncolumns, int nrows, int xwidth, int ywidth,
                     double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  GridCanvas* gc=new GridCanvas(uniq(), ncolumns, nrows, xwidth, ywidth);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);
      proj->SetFillColor(h2d->GetFillColor());
      proj->SetFillStyle(h2d->GetFillStyle());


      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
GridCanvas* plotXAxis1DRebinPz(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, int ncolumns, int nrows, int xwidth, int ywidth,
                     double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  GridCanvas* gc=new GridCanvas(uniq(), ncolumns, nrows, xwidth, ywidth);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();


      TH1* projtmp=h2d->ProjectionX(uniq(), i+1, i+1);
      TH1* proj=rebinpz(projtmp);


      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }


    

    //BINdrawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);
    drawBinRange(histAndOpts[0].first, 1, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
        //gcPT->SetYLimits(0, 0.29);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  const int nLabels=5;
  const double positions[nLabels]={4,10,16,20,24};
  const char* valueStrings[nLabels]={"4","10","20", "40", "60"};
  
  //gc->SetManualXLabels(nLabels, positions, valueStrings,0.03);
  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotXAxis1D_Nue(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle,
                     double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x-1);
  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  cout << "Plotting plotXAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), 2, 3, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<5; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);
      proj->GetXaxis()->SetRangeUser(0,100);
      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotpz1DAntiNu(std::vector<std::pair<TH2*, const char*> > histAndOpts,
                     double* multipliers=NULL)
{
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), 3, 2, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<6; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);
      //      TH1* proj=rebinpz(projtmp);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, "p_{T}/GeV", ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Muon longitudinal momentum (GeV)");
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/nucleon)");
  gc->SetTitleSize(22.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotvtx1D(std::vector<std::pair<TH2*, const char*> > histAndOpts,
		      double* multipliers=NULL, TH2* denominator=NULL, bool track2=false)
{
  cout << "Doing ratios" << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=track2? new GridCanvas(uniq(), 3, 2, 4000, 2000): new GridCanvas(uniq(), 3, 3, 4000, 2000);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();
  
  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them


  const int nbins = track2 ? 4:8;  
  int minb2[4]={1,3,6,9};
  int maxb2[4]={2,5,8,13};

  int minb[8]={1,3,4,5,6,7,8,9};
  int maxb[8]={2,3,4,5,6,7,8,13};

  
  for(int i=0; i<nbins; ++i){
    // Change to the appropriate pad in the canvas
    
    int val = i>1?i+2:i+1;
    if(!track2)val=i+1;
    TPad* pad=(TPad*)gc->cd(val);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();
      TH1* projtmp = track2? h2d->ProjectionX(uniq(), minb2[i],maxb2[i]):h2d->ProjectionX(uniq(), minb[i],maxb[i]);
      TH1* demproj = track2? denominator->ProjectionX(uniq(), minb2[i],maxb2[i]):denominator->ProjectionX(uniq(), minb[i],maxb[i]);
      projtmp->Divide(demproj);

      TH1* proj=(TH1*)projtmp->Clone();
      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
	if(opt.Contains("histe")){
	  proj->SetFillColor(kWhite);
	  proj->SetFillStyle(0);
	  proj->DrawClone("histc");
	  proj->SetLineWidth(2);
	  proj->SetFillColor(kRed);
	  proj->SetFillStyle(3003);
	  proj->Draw("E5SAME");
	}
	else proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
	  
      }

    }

    if(track2) drawVarBinRange(histAndOpts[0].first, 2, minb2[i],maxb2[i], "p_{T}/GeV", ".2f", false);
    else  drawVarBinRange(histAndOpts[0].first, 2, minb[i],maxb[i], "p_{T}/GeV", ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Vertex Energy");
  gc->SetYTitle("Events");
  gc->SetTitleSize(100.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotvtxrate1D(std::vector<std::pair<TH2*, const char*> > histAndOpts,
			  double* multipliers=NULL,bool track2=false, bool binbybin =false)
{
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=track2? new GridCanvas(uniq(), 3, 2, 4000, 2000): new GridCanvas(uniq(), 3, 3, 4000, 2000);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();
  
  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them


  const int nbins = track2 ? 4:8;  
  int minb2[4]={1,3,6,9};
  int maxb2[4]={2,5,8,13};

  int minb[8]={1,3,4,5,6,7,8,9};
  int maxb[8]={2,3,4,5,6,7,8,13};


  
  for(int i=0; i<nbins; ++i){
    // Change to the appropriate pad in the canvas
    
    int val = i>1?i+2:i+1;
    if(!track2)val=i+1;
    TPad* pad=(TPad*)gc->cd(val);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      //      h2d->SaveAs(Form("%s-%s-%d.root",uniq().Data(),h2d->GetName(),track2));
      opt.ToLower();
      TH1* proj= NULL;
      if(!binbybin) proj = track2? h2d->ProjectionX(uniq(), minb2[i],maxb2[i]):h2d->ProjectionX(uniq(), minb[i],maxb[i]);
      else proj = h2d->ProjectionX(uniq(),i+1,i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);
      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
	if(opt.Contains("histe")){
	  proj->SetFillColor(kWhite);
	  proj->SetFillStyle(0);
	  proj->DrawClone("histc");
	  proj->SetLineWidth(2);
	  proj->SetFillColor(kRed);
	  proj->SetFillStyle(3003);
	  proj->Draw("E5SAME");
	}
	else proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
	  
      }

    }
    if(!binbybin){
      if(track2) drawVarBinRange(histAndOpts[0].first, 2, minb2[i],maxb2[i], "p_{T}/GeV", ".2f", false);
      else  drawVarBinRange(histAndOpts[0].first, 2, minb[i],maxb[i], "p_{T}/GeV", ".2f", false);
    }
    else{
      drawBinRange(histAndOpts[0].first, 2, i+1, "p_{T}/GeV", ".2f", false);
    }

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Vertex Energy");
  gc->SetYTitle("Events");
  gc->SetTitleSize(100.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

vector<double> GetScales(std::vector<std::pair<TH2*, const char*> >histopts, bool pxProj, double plotMax, double fracMax, bool limitRange = false, double rmin=0, double rmax=100){
  vector<double> tmpvect;
  int nbins = histopts[0].first->GetNbinsX()+2;
  if(pxProj) nbins = histopts[0].first->GetNbinsY()+2;
  for(int i=1;i<nbins-1;i++){
    double maxval = 0;
    for(uint j=0;j<histopts.size();j++){
      TH1D *tmp = pxProj? histopts[j].first->ProjectionX("tmp",i,i): histopts[j].first->ProjectionY("tmp",i,i);
      if(limitRange) tmp->GetXaxis()->SetRangeUser(rmin,rmax);
      int maxbin = tmp->GetMaximumBin();
      double content = tmp->GetBinContent(maxbin);
      if(content>maxval) maxval=content;
    }
    //we want about fracMax of the plotMax. 
    double scale = (plotMax*fracMax)/maxval;
    if(scale>1){
      int tmpscale = floor(scale*10);
      scale = tmpscale/10.0;
    }
    else{
      int tmpscale = ceil(scale*10);
      scale = tmpscale/10.0;
    }
    cout << scale << endl;
    tmpvect.push_back(scale);
  }
    return tmpvect;
}




#endif
