import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle

ROOT.gROOT.SetBatch(True)

ROOT.TH1.AddDirectory(False)
dirpwd = sys.argv[1]
plist = sys.argv[2]

infile= ROOT.TFile(str(dirpwd)+"/Migration2D_Daisy_%s_t99_z99_sys.root"%(plist))
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()
mnv.SetRedHeatPalette()
#mnv.SetCorrelationPalette()

trueZ = "Tracker"
mat = "CH"

vars = ["pZmu_pTmu"] #"ThetamuDeg"]

for var in vars:
        for petal in range(0,12):
                print("---------------------------------------")
                print(var)
                mc_hist = infile.Get("response2d_daisy_%d_%s_migration"%(petal,var))

                gStyle.SetOptTitle(0)
                gStyle.SetPaintTextFormat("2.0f")

                # ROW NORMALIZED
                mnv.DrawNormalizedMigrationHistogram(mc_hist, True, False, True, True)
                # Each bigger box is a pII bin while the smaller boxes in the big one pT
                
                keys = canvas1.GetListOfPrimitives();

                for k in keys:
                        if(k.ClassName().find("TH2D")!=-1): # to change axis titles
                                if var == "pZmu_pTmu":
                                        k.GetXaxis().SetTitle("Reconstructed Muon P_{||} per Muon P_{t} bins")
                                        k.GetYaxis().SetTitle("True Muon P_{||} per Muon P_{t} bins")

                                k.GetZaxis().SetTitle("Row Normalized Event Rate (%)")
                                
                                k.GetXaxis().CenterTitle()
                                k.GetXaxis().SetTitleFont(42)
                                k.GetXaxis().SetTitleSize(0.05)

                                k.GetYaxis().SetTitleOffset(1)
                                k.GetYaxis().CenterTitle()
                                k.GetYaxis().SetTitleFont(42)
                                k.GetYaxis().SetTitleSize(0.05)

                                k.GetZaxis().CenterTitle()
                                k.GetZaxis().SetRangeUser(0,100)
                                k.GetZaxis().SetTitleFont(42)
                                k.GetZaxis().SetTitleOffset(1.1)
                                k.GetZaxis().SetTitleSize(0.045)

                                k.GetXaxis().SetTitleSize(0.05)
                                k.GetYaxis().SetTitleSize(0.05)
                                k.GetYaxis().SetTitleOffset(1.)
                                k.GetZaxis().SetTitleSize(0.04)
                                k.GetZaxis().SetTitleOffset(1.2)
                                k.GetZaxis().SetLabelSize(0.04)
                '''
                void MnvPlotter::DrawNormalizedMigrationHistogram(
                        const TH2D* h_migration,
                        const bool drawAsMatrix,
                        const bool coarseContours, /* = false */
                        const bool includeFlows, /* = true */
                        const bool noText /* = false */
                        )
                {
                '''
           
                mnv.AddHistoTitle("Daisy %s: Petal %d"%(trueZ, petal), 0.05, 1)
                canvas1.Print("Migration2D_daisy_%ss_t99_z99_%s_p%d_RowNorm.png"%(plist, var, petal))

print("DONE %s 99 99 Daisy"%(plist))
