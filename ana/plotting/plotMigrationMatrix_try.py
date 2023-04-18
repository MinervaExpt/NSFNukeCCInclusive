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
targetID = sys.argv[2] 
targetZ = sys.argv[3]
plist = sys.argv[4]

infile= ROOT.TFile(str(dirpwd)+"/Migration_%s_t%s_z%02s_sys_p4_NukeCC_AntiNu_Tgt5_Fe_MCTuning.root"%(plist, targetID, targetZ))
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mat = None
trueZ = None

if targetZ == "26":
    trueZ = "Iron"
    mat = "Fe"

if targetZ == "82":
    trueZ = "Lead"
    mat = "Pb"

if targetZ == "06":
    trueZ = "Carbon"
    mat = "C"

if targetZ == "99":
    trueZ = "Tracker"
    mat = "CH"

vars = [ "Ehad"] #"ThetamuDeg"]

for var in vars:
        print("---------------------------------------")
        print(var)
        mc_hist = infile.Get("response1d_%s_migration"%var)
        new = infile.Get("selected_Migration_%s"%var)

        # to add overflow + underflow
        mc_hist.GetXaxis().SetRangeUser(-1,mc_hist.GetNbinsX()+1)
        mc_hist.GetYaxis().SetRangeUser(-1,mc_hist.GetNbinsY()+1)

        gStyle. SetPalette(1)
        gStyle.SetOptTitle(0)
        gStyle.SetPaintTextFormat("2.0f")

        if var == "Enu":
                mc_hist.GetXaxis().SetTitle("Reconstructed E_{#bar{#nu}} bin number")
                mc_hist.GetYaxis().SetTitle("True E_{#bar{#nu}} bin number")
                mc_hist.SetMarkerSize(1.2)
                new.GetXaxis().SetTitle("Reconstructed E_{#bar{#nu}}")
                new.GetYaxis().SetTitle("True E_{#bar{#nu}}")
                new.SetMarkerSize(1.2)


        if var == "x":
                mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x bin number")
                mc_hist.GetYaxis().SetTitle("True Bjorken x bin number")
                mc_hist.SetMarkerSize(1.5)
                new.GetXaxis().SetTitle("Reconstructed Bjorken x")
                new.GetYaxis().SetTitle("True Bjorken x")
                new.SetMarkerSize(1.2)
        
        if var == "pTmu1D":
                mc_hist.GetXaxis().SetTitle("Reconstructed Muon p_{T} bin number")
                mc_hist.GetYaxis().SetTitle("True Muon p_{T} bin number")
                mc_hist.SetMarkerSize(1.5)
                new.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
                new.GetYaxis().SetTitle("True Muon p_{T} (GeV/c)")
                new.SetMarkerSize(1.2)
                
        if var == "pZmu1D":
                mc_hist.GetXaxis().SetTitle("Reconstructed Muon p_{Z} bin number")
                mc_hist.GetYaxis().SetTitle("True Muon p_{Z} bin number")
                mc_hist.SetMarkerSize(1.5)
                new.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)")
                new.GetYaxis().SetTitle("True Muon p_{Z} (GeV/c)")
                new.SetMarkerSize(1.2)

        if var == "ThetamuDeg":
                mc_hist.GetXaxis().SetTitle("Reconstructed Muon #theta_{#mu} bin number")
                mc_hist.GetYaxis().SetTitle("True Muon #theta_{#mu} bin number")
                mc_hist.SetMarkerSize(1.5)
                new.GetXaxis().SetTitle("Reconstructed Muon #theta_{#mu} (Deg)")
                new.GetYaxis().SetTitle("True Muon #theta_{#mu} (Deg)")
                new.SetMarkerSize(1.2)

        if var == "Ehad":
                mc_hist.GetXaxis().SetTitle("Reconstructed Recoil bin number")
                mc_hist.GetYaxis().SetTitle("True Recoil bin number")
                mc_hist.SetMarkerSize(1.5)
                new.GetXaxis().SetTitle("Reconstructed Recoil (GeV)")
                new.GetYaxis().SetTitle("True Recoil (GeV)")
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
        mc_hist.GetXaxis().CenterTitle()
        mc_hist.GetXaxis().SetTitleFont(42)
        mc_hist.GetXaxis().SetTitleSize(0.05)

        mc_hist.GetYaxis().SetTitleOffset(1)
        mc_hist.GetYaxis().CenterTitle()
        mc_hist.GetYaxis().SetTitleFont(42)
        mc_hist.GetYaxis().SetTitleSize(0.05)

        mc_hist.GetZaxis().CenterTitle()
        mc_hist.GetZaxis().SetRangeUser(0,1000)
        mc_hist.GetZaxis().SetTitleFont(42)
        mc_hist.GetZaxis().SetTitleOffset(1.1)
        mc_hist.GetZaxis().SetTitleSize(0.045)

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

        mc_hist.GetXaxis().SetTitleSize(0.05)
        mc_hist.GetYaxis().SetTitleSize(0.05)
        mc_hist.GetYaxis().SetTitleOffset(1.)
        mc_hist.GetZaxis().SetTitleSize(0.04)
        mc_hist.GetZaxis().SetTitleOffset(1.2)
        mc_hist.GetZaxis().SetLabelSize(0.04)

        new.GetXaxis().SetTitleSize(0.05)
        new.GetYaxis().SetTitleSize(0.05)
        new.GetYaxis().SetTitleOffset(1.)
        new.GetZaxis().SetTitleSize(0.04)
        new.GetZaxis().SetTitleOffset(1.2)
        new.GetZaxis().SetLabelSize(0.04)

        # SIMPLE OCCUPANCY
        mc_hist.Draw("COLZTEXT")
        if targetZ == "99":
                mnv.AddHistoTitle("%s"%(trueZ), 0.05, 1)
        else:
                mnv.AddHistoTitle("Target %s %s"%(targetID, trueZ), 0.05, 1)
        mnv.SetRedHeatPalette()
        #mnv.SetROOT6Palette(109)
        mc_hist.GetZaxis().SetTitle("Event Rate")
        canvas1.Print("Migration_%s_t%s_z%02s_%s_Occupancy_p4_AntiNu_Tgt5_Fe_MCTuning.png"%(plist, targetID, targetZ, var))

        new.Draw("COLZ")
        if targetZ == "99":
                mnv.AddHistoTitle("%s"%(trueZ), 0.05, 1)
        else:
                mnv.AddHistoTitle("Target %s %s"%(targetID, trueZ), 0.05, 1)
        new.GetZaxis().SetTitle("Event Rate")
        canvas1.Print("Migration_%s_t%s_z%02s_%s_Occupancy_Binned_p4_AntiNu_Tgt5_Fe_MCTuning.png"%(plist, targetID, targetZ, var))

        # COLUMN NORMALIZED
        mc_hist.GetZaxis().SetRangeUser(0,100)
        mc_hist.GetZaxis().SetTitle("Column Normalized Event Rate (%)")
        colnorm_hist = colNormalize(mc_hist)
        colnorm_hist.Draw("COLZTEXT")
        if targetZ == "99":
                mnv.AddHistoTitle("%s"%(trueZ), 0.05, 1)
        else:
                mnv.AddHistoTitle("Target %s %s"%(targetID, trueZ), 0.05, 1)
        canvas1.Print("Migration_%s_t%s_z%02s_%s_ColumnNorm_p4_AntiNu_Tgt5_Fe_MCTuning.png"%(plist, targetID, targetZ, var))


        new.GetZaxis().SetRangeUser(0,100)
        new.GetZaxis().SetTitle("Column Normalized Event Rate (%)")
        colnorm_hist_new = colNormalize(new)
        colnorm_hist_new.Draw("COLZ")
        if targetZ == "99":
                mnv.AddHistoTitle("%s"%(trueZ), 0.05, 1)
        else:
                mnv.AddHistoTitle("Target %s %s"%(targetID, trueZ), 0.05, 1)
        canvas1.Print("Migration_%s_t%s_z%02s_%s_ColumnNorm_Binned_p4_AntiNu_Tgt5_Fe_MCTuning.png"%(plist, targetID, targetZ, var))

        # ROW NORMALIZED
        mc_hist.GetZaxis().SetTitle("Row Normalized Event Rate (%)")
        rownorm_hist = rowNormalize(mc_hist)
        rownorm_hist.Draw("COLZTEXT")
        if targetZ == "99":
                mnv.AddHistoTitle("%s"%(trueZ), 0.05, 1)
        else:
                mnv.AddHistoTitle("Target %s %s"%(targetID, trueZ), 0.05, 1)
        canvas1.Print("Migration_%s_t%s_z%02s_%s_RowNorm_p4_AntiNu_Tgt5_Fe_MCTuning.png"%(plist, targetID, targetZ, var))

        new.GetZaxis().SetTitle("Row Normalized Event Rate (%)")
        rownorm_hist_new = rowNormalize(new)
        rownorm_hist_new.Draw("COLZ")
        if targetZ == "99":
                mnv.AddHistoTitle("%s"%(trueZ), 0.05, 1)
        else:
                mnv.AddHistoTitle("Target %s %s"%(targetID, trueZ), 0.05, 1)
        canvas1.Print("Migration_%s_t%s_z%02s_%s_RowNorm_Binned_p4_AntiNu_Tgt5_Fe_MCTuning.png"%(plist, targetID, targetZ, var))

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

print("DONE %s %s %02s"%(plist, targetID, targetZ))
