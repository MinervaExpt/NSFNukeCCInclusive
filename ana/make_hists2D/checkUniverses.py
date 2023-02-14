import ROOT
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis, gPad, TLine
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
import os,sys
from ROOT import PlotUtils

rootfile = sys.argv[1] #sys.argv is the list of command line arguments passed to a Python script. argv[0] is the script name

infile = ROOT.TFile(str(rootfile))
canvas1 = ROOT.TCanvas() #need to declare canvas before calling MnvPlotter
mnv = PlotUtils.MnvPlotter()

ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)

response2D_hist = infile.Get("response2d_daisy_0_pZmu_pTmu_migration")
vertbands = response2D_hist.GetVertErrorBandNames()
for band in vertbands:
  print(band)
  universes = response2D_hist.GetVertErrorBand(band).GetNHists()
  for univ in range(universes):
    hist = response2D_hist.GetVertErrorBand(band).GetHist(univ)
    hist.SetEntries(1)
    hist.Draw()
    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("%s_universe%s.png"%(band, univ))