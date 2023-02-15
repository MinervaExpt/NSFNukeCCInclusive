import ROOT
import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle

ROOT.gROOT.SetBatch(True)

def colNormalize(hist):
        norm_hist = hist.Clone()
        norm_hist.Reset()

        nXBins = hist.GetNbinsX()+2 # to add overflow + underflow
        nYBins = hist.GetNbinsY()+2 # to add overflow + underflow
        
        for i in range(0, nXBins): # column
                col_integral = 0
                for j in range(0,nYBins): # row
                        col_integral += hist.GetBinContent(i,j)
                for j in range(0,nYBins):
                        if col_integral==0: norm_hist.SetBinContent(i,j,0)
                        else: norm_hist.SetBinContent(i,j, hist.GetBinContent(i,j)*100/col_integral)
        return norm_hist


def rowNormalize(hist):
        norm_hist = hist.Clone()
        norm_hist.Reset()

        nXBins = hist.GetNbinsX()+2 # to add overflow + underflow
        nYBins = hist.GetNbinsY()+2 # to add overflow + underflow
        
        for i in range(0, nYBins): # column
                row_integral = 0
                for j in range(0,nXBins): # row
                        row_integral += hist.GetBinContent(j,i)
                for j in range(0,nXBins):
                        if row_integral==0: norm_hist.SetBinContent(j,i,0)
                        else: norm_hist.SetBinContent(j,i, hist.GetBinContent(j,i)*100/row_integral)
        return norm_hist


ROOT.TH1.AddDirectory(False)
dirpwd = sys.argv[1]
plist = sys.argv[2]

infile= ROOT.TFile(str(dirpwd)+"/Migration_Daisy_%s_t99_z99_sys.root"%(plist))
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

vars = ["Enu", "x", "pTmu", "pZmu"] #"ThetamuDeg"]

for var in vars:
        for petal in range(0,12):
                print("---------------------------------------")
                print(var)
                new = infile.Get("selected_Migration_daisy_%s_%s"%(petal, var))

                # to add overflow + underflow

                gStyle.SetPalette(1)
                gStyle.SetOptTitle(0)
                gStyle.SetPaintTextFormat("2.0f")

                if var == "Enu":
                        new.GetXaxis().SetTitle("Reconstructed E_{#bar{#nu}}")
                        new.GetYaxis().SetTitle("True E_{#bar{#nu}}")
                        new.SetMarkerSize(1.2)


                if var == "x":
                        new.GetXaxis().SetTitle("Reconstructed Bjorken x")
                        new.GetYaxis().SetTitle("True Bjorken x")
                        new.SetMarkerSize(1.2)

                if var == "pTmu":
                        new.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
                        new.GetYaxis().SetTitle("True Muon p_{T} (GeV/c)")
                        new.SetMarkerSize(1.2)
                
                if var == "pZmu":
                        new.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)")
                        new.GetYaxis().SetTitle("True Muon p_{Z} (GeV/c)")
                        new.SetMarkerSize(1.2)

                '''
                if myvariable == "x09":
                        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x Bins")
                        mc_hist.GetYaxis().SetTitle("True Bjorken x Bins")

                if myvariable == "xfine":
                        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x Bins")
                        mc_hist.GetYaxis().SetTitle("True Bjorken x Bins")

                if myvariable == "xBrian":
                        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x Bins")
                        mc_hist.GetYaxis().SetTitle("True Bjorken x Bins")
                '''


                new.GetXaxis().CenterTitle()
                new.GetXaxis().SetTitleFont(42)
                new.GetXaxis().SetTitleSize(0.05)

                new.GetYaxis().SetTitleOffset(1)
                new.GetYaxis().CenterTitle()
                new.GetYaxis().SetTitleFont(42)
                new.GetYaxis().SetTitleSize(0.05)

                new.GetZaxis().CenterTitle()
                new.GetZaxis().SetRangeUser(0,1000)
                new.GetZaxis().SetTitleFont(42)
                new.GetZaxis().SetTitleOffset(1.1)
                new.GetZaxis().SetTitleSize(0.045)


                '''
                mc_hist.GetXaxis().SetBinLabel(0,'underflow')
                mc_hist.GetXaxis().SetBinLabel(1,'0.0 - 0.05')
                mc_hist.GetXaxis().SetBinLabel(2,'0.05 - 0.1')
                mc_hist.GetXaxis().SetBinLabel(3,'0.1 - 0.2')
                mc_hist.GetXaxis().SetBinLabel(4,'0.2 - 0.4')
                mc_hist.GetXaxis().SetBinLabel(5,'0.4 - 1.0')
                mc_hist.GetXaxis().SetBinLabel(6,'1.0 - 2.2')
                mc_hist.GetXaxis().SetBinLabel(7,' > 2.2')
                #{0.0, 0.05, 0.1, 0.2, 0.4, 1.0, 2.2}

                mc_hist.GetYaxis().SetBinLabel(0,'underflow')
                mc_hist.GetYaxis().SetBinLabel(1,'0.0 - 0.05')
                mc_hist.GetYaxis().SetBinLabel(2,'0.05 - 0.1')
                mc_hist.GetYaxis().SetBinLabel(3,'0.1 - 0.2')
                mc_hist.GetYaxis().SetBinLabel(4,'0.2 - 0.4')
                mc_hist.GetYaxis().SetBinLabel(5,'0.4 - 1.0')
                mc_hist.GetYaxis().SetBinLabel(6,'1.0 - 2.2')
                mc_hist.GetYaxis().SetBinLabel(7,' > 2.2')
                '''


                new.GetXaxis().SetTitleSize(0.05)
                new.GetYaxis().SetTitleSize(0.05)
                new.GetYaxis().SetTitleOffset(1.)
                new.GetZaxis().SetTitleSize(0.04)
                new.GetZaxis().SetTitleOffset(1.2)
                new.GetZaxis().SetLabelSize(0.04)

                # SIMPLE OCCUPANCY
                new.Draw("COLZ")
                mnv.SetRedHeatPalette()
                mnv.AddHistoTitle("Daisy Tracker Petal %s"%(petal), 0.05, 1)
                new.GetZaxis().SetTitle("Event Rate")
                canvas1.Print("Migration_daisy_%s_t99_z99_%s_p%s_Occupancy_Binned.png"%(plist, var, petal))

                # COLUMN NORMALIZED
                new.GetZaxis().SetRangeUser(0,100)
                new.GetZaxis().SetTitle("Column Normalized Event Rate (%)")
                colnorm_hist_new = colNormalize(new)
                colnorm_hist_new.Draw("COLZ")
                mnv.AddHistoTitle("Daisy Tracker Petal %s"%(petal), 0.05, 1)
                canvas1.Print("Migration_daisy_%s_t99_z99_%s_p%s_ColumnNorm_Binned.png"%(plist, var, petal))

                # ROW NORMALIZED
                new.GetZaxis().SetTitle("Row Normalized Event Rate (%)")
                rownorm_hist_new = rowNormalize(new)
                rownorm_hist_new.Draw("COLZ")
                mnv.AddHistoTitle("Daisy Tracker Petal %s"%(petal), 0.05, 1)
                canvas1.Print("Migration_daisy_%s_t99_z99_%s_p%s_RowNorm_Binned.png"%(plist, var, petal))

#mnv.DrawNormalizedMigrationHistogram(mc_hist, True, False, True, False)
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

print("DONE %s 99 99 Daisy"%(plist))
