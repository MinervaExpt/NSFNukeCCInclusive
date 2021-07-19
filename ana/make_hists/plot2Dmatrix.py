import ROOT
import os,sys
#import PlotUtils
from PlotUtils import MnvH1D,MnvH2D,MnvPlotter
from ROOT import gStyle

infile= ROOT.TFile("Hists_Migration_t1235_z26_AntiNu_v1_nosys.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = MnvPlotter()

#mcPOT = infile.Get("MCPOT").GetVal()
#dataPOT = infile.Get("DataPOT").GetVal()

#mcScale =  dataPOT/mcPOT
gStyle. SetPalette(1)
gStyle.SetOptTitle(0)

mc_hist = infile.Get("response1d_Enu_migration")
mc_hist.Draw("COLZ TEXTnn")
mc_hist.GetXaxis().SetTitle("Reconstructed Neutrino Energy [GeV]")
mc_hist.GetXaxis().CenterTitle()
mc_hist.GetYaxis().SetTitle("True Neutrino Energy [GeV]")
mc_hist.GetYaxis().CenterTitle()
mc_hist.SetTitle("All Iron Targets")
mnv.AddHistoTitle("All Iron Targets", 0.05, 1)
canvas1.Modified()

canvas1.Print("ME6A_T1235Fe_Migration_nosys.png")

#canvas1.Print("NukeCC_MEJ_Emu_T1.png")

raw_input("Done")
