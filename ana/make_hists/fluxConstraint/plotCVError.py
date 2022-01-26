import ROOT
import os,sys
from ROOT import PlotUtils

infile= ROOT.TFile("/minerva/app/users/anezkak/MAT_GitHub/NSFNukeCCInclusive/ana/make_hists/fluxConstraint/Hists_EventSelection_minervame1A_FluxConstraint_sys_t99_z99_Nu.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

#mcPOT = infile.Get("MCPOT").GetVal()

#myfile = ROOT.TFile(sys.argv[1])
#myvariable = sys.argv[1]

mc_hist = infile.Get("selected_mc_truth_Enu")

#if myvariable == "Enu":
#    mc_hist.GetXaxis().SetTitle("Neutrino Energy [GeV]")
#    mc_hist.GetYaxis().SetTitle("Events/GeV")

#if myvariable == "vtxz":
#    mc_hist.GetXaxis().SetTitle("Vertex z [cm]")
#    mc_hist.GetYaxis().SetTitle("Events/cm")

mc_hist.GetXaxis().CenterTitle()
mc_hist.GetYaxis().CenterTitle()
print(mc_hist.GetNbinsX())
#mc_hist.GetXaxis().SetRangeUser(0,0.001)

mnv.DrawMCWithErrorBand(mc_hist)

#ndf = 0
#chi2 = mnv.Chi2DataMC(data_hist, mc_hist, ndf, mcScale, False, True)
#ndf = mc_hist.GetNbinsX()
#print(chi2)
#print(ndf)

#mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TR", 0.04, 0.0, True, False)

'''
Chi2DataMC( dataHist, mcHist, ndf, mcScale, useDataErrorMatrix, useOnlyShapeErrors)
If your data distribution has ONLY STATISTICAL uncertainties, useDataErrorMatrix should be FALSE
If you are using POT-normalization, then useOnlyShapeErrors should be FALSE
'''

mnv.AddHistoTitle("Tracker", 0.05, 1)
#mnv.WritePreliminary(0.66, 0.80, 0.035, 0.11);
#mnv.AddPOTNormBox(dataPOT,mcPOT, 0.7, 0.8)

canvas1.Modified()
canvas1.Print("ME1A_Enu_sys.png")

mc_hist.SetDirectory(0)

raw_input("Done")

