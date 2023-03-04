import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLegend
from ROOT import TLine

ROOT.gROOT.SetBatch(True)

dirpwd = sys.argv[1]
targetZ = sys.argv[2]
plist = sys.argv[3]
scale = sys.argv[4]

mat = None
trueZ = None

if targetZ == "26":
    trueZ = "Iron"
    mat = "Fe"

if targetZ == "82":
    trueZ = "Lead"
    mat = "Pb"


if targetZ == "26":
    target2 = ROOT.TFile(str(dirpwd)+"/CrossSection_Daisy_t2_z26_%s.root"%(plist))
    target3 = ROOT.TFile(str(dirpwd)+"/CrossSection_Daisy_t3_z26_%s.root"%(plist))
    target5 = ROOT.TFile(str(dirpwd)+"/CrossSection_Daisy_t5_z26_%s.root"%(plist))
if targetZ == "82":
    target2 = ROOT.TFile(str(dirpwd)+"/CrossSection_Daisy_t2_z82_%s.root"%(plist))
    target3 = ROOT.TFile(str(dirpwd)+"/CrossSection_Daisy_t3_z82_%s.root"%(plist))
    target4 = ROOT.TFile(str(dirpwd)+"/CrossSection_Daisy_t4_z82_%s.root"%(plist))
    target5 = ROOT.TFile(str(dirpwd)+"/CrossSection_Daisy_t5_z82_%s.root"%(plist))

canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = target2.Get("MCPOT").GetVal()
dataPOT = target2.Get("DataPOT").GetVal()

mcScale = dataPOT/mcPOT
if scale == "1":
    mcScale = 1

vars = ["x"] #, "pZmu1D", "pTmu1D", "ThetamuDeg"]

#steps = ['unfolded','unfolded_effCorrected', 'crossSection', 'crossSection_total']
steps = ['crossSection']

for step in steps:

    for var in vars:

        data_hist2 = target2.Get("%s_data_%s"%(step, var))
        data_hist3 = target3.Get("%s_data_%s"%(step, var))
        data_hist5 = target5.Get("%s_data_%s"%(step, var))
        if targetZ == "82":
            data_hist4 = target4.Get("%s_data_%s"%(step, var))

        mc_hist2 = target2.Get("simEventRate_%s_mc_%s"%(step, var))
        mc_hist3 = target3.Get("simEventRate_%s_mc_%s"%(step, var))
        mc_hist5 = target5.Get("simEventRate_%s_mc_%s"%(step, var))
        if targetZ == "82":
            mc_hist4 = target4.Get("simEventRate_%s_mc_%s"%(step, var))


        if var == "Enu":
            mc_hist2.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
            if step == "crossSection":
                mc_hist2.GetYaxis().SetTitle("d#sigma/dE_{#bar{#nu}} (10^{-39} cm^{2}/GeV/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist2.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist2.GetYaxis().SetTitle("Events #times 10^{3}/GeV")

        if var == "x":
            mc_hist2.GetXaxis().SetTitle("Bjorken x")
            if step == "crossSection":
                mc_hist2.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist2.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist2.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")

        if var == "pTmu1D":
            mc_hist.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
            if step == "crossSection":
                mc_hist2.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist2.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist2.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")
       
        if var == "pZmu1D":
            mc_hist2.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)")
            if step == "crossSection":
                mc_hist2.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist2.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist2.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")

        if var == "ThetamuDeg":
            mc_hist2.GetXaxis().SetTitle("Reconstructed Muon Angle (Deg)")
            if step == "crossSection":
                mc_hist2.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist2.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist2.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")


        mc_hist2.GetXaxis().CenterTitle()
        mc_hist2.GetYaxis().CenterTitle()
        #print(data_hist.GetBinContent(data_hist.GetNbinsX())) # last bin
        #print(data_hist.GetBinContent(data_hist.GetNbinsX()+1)) # overflow
        #mc_hist.GetXaxis().SetRangeUser(66,173)

        #mnv.ApplyStyle(1)
        #if var == "Enu":
        #    mc_hist.GetXaxis().SetRangeUser(2, 20)

        data_hist_stat2 =  data_hist2.GetCVHistoWithStatError() # stat error
        data_hist_total2 = data_hist2.GetCVHistoWithError() # total error
        data_hist_sys2 = data_hist2.GetCVHistoWithError(False) # sys error (bool is include stat)
        mc_hist_stat2 = mc_hist2.GetCVHistoWithStatError() 

        data_hist_stat3 =  data_hist3.GetCVHistoWithStatError() # stat error
        data_hist_total3 = data_hist3.GetCVHistoWithError() # total error
        data_hist_sys3 = data_hist3.GetCVHistoWithError(False) # sys error (bool is include stat)
        mc_hist_stat3 = mc_hist3.GetCVHistoWithStatError() 

        data_hist_stat5 =  data_hist5.GetCVHistoWithStatError() # stat error
        data_hist_total5 = data_hist5.GetCVHistoWithError() # total error
        data_hist_sys5 = data_hist5.GetCVHistoWithError(False) # sys error (bool is include stat)
        mc_hist_stat5 = mc_hist5.GetCVHistoWithStatError() 

        if targetZ == "82":
            data_hist_stat4 =  data_hist4.GetCVHistoWithStatError() # stat error
            data_hist_total4 = data_hist4.GetCVHistoWithError() # total error
            data_hist_sys4 = data_hist4.GetCVHistoWithError(False) # sys error (bool is include stat)
            mc_hist_stat4 = mc_hist4.GetCVHistoWithStatError() 

        
        # MC
        mc_hist2.SetLineWidth(3)
        mc_hist2.SetLineColor(2)
        # MC error
        mc_hist_stat2.SetFillColor(ROOT.kRed-10)
        mc_hist_stat2.SetFillStyle(1001)
        mc_hist_stat2.SetMarkerStyle(0)
        mc_hist_stat2.SetLineWidth(3)
        mc_hist_stat2.SetLineColor(2)

        mc_hist2.SetLineColor(ROOT.kRed)
        mc_hist2.SetLineWidth(2)

        mc_hist2.SetMaximum(data_hist2.GetMaximum()*1.5)

        
        ratio2 = data_hist_stat2.Clone()
        ratio2.Divide(ratio2,mc_hist_stat2) # stat
        if var == "x":
            ratio2.GetXaxis().SetTitle("Bjorken x")
        if var == "pTmu1D":
            ratio2.GetXaxis().SetTitle("Muon p_{T} (GeV/c)")
        if var == "pZmu1D":
            ratio2.GetXaxis().SetTitle("Muon p_{Z} (GeV/c)")
        if var == "ThetamuDeg":
            ratio2.GetXaxis().SetTitle("Muon #theta_{#mu} (Deg)")
        ratio2.GetYaxis().SetTitle("Data/MC")
        ratio2.GetYaxis().CenterTitle()
        ratio2.GetXaxis().CenterTitle()
        ratio2.Draw("X0")
  
        ratio2.SetMaximum(1.5)
        ratio2.SetMinimum(0.5)
        ratio2.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions
        # Copied from MINERvA-101-Cross-Section/backgroundStack.py
        # same as in void MnvPlotter::DrawDataMCRatio() in MnvPlotter
        # systematic error centered at y = 1
        ratio_tot2 = data_hist_total2.Clone()
        ratio_tot2.Divide(ratio_tot2, mc_hist_stat2)

        # line at 1
        line = TLine(ratio2.GetXaxis().GetXmin(), 1, ratio2.GetXaxis().GetXmax(), 1)
        line.SetLineColor(46)
        line.SetLineWidth(2)
        line.SetLineStyle(9)
        line.Draw()
        
        ratio2.SetLineColor(ROOT.kMagenta+2)
        ratio2.SetMarkerColor(ROOT.kMagenta+2)
        ratio2.Draw("X0 SAME E1")
        #ratio_tot2.Draw("E1 SAME")

        ratio3 = data_hist_stat3.Clone()
        ratio3.Divide(ratio3,mc_hist_stat3) # stat
        ratio3.SetLineColor(ROOT.kYellow+3)
        ratio3.SetMarkerColor(ROOT.kYellow+3)
        ratio3.Draw("X0 SAME E1")


        if targetZ == "82":
            ratio4 = data_hist_stat4.Clone()
            ratio4.Divide(ratio4,mc_hist_stat4) # stat
            ratio4.SetLineColor(ROOT.kRed-7)
            ratio4.SetMarkerColor(ROOT.kRed-7)
            ratio4.Draw("X0 SAME E1")


        ratio5 = data_hist_stat5.Clone()
        ratio5.Divide(ratio5,mc_hist_stat5) # stat
        ratio5.SetLineColor(ROOT.kCyan+1)
        ratio5.SetMarkerColor(ROOT.kCyan+1)
        ratio5.Draw("X0 SAME E1")


        #chi2 = mnv.Chi2DataMC(data_hist, mc_hist,mcScale,False,True )
        #ndf = mc_hist.GetNbinsX()
        #print(chi2)
        if step == "unfolded":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.7, 0.32)
        elif step == "total_unfolded_effCorrected":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.7, 0.32)
        elif step == "crossSection_total":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.7, 0.32)
        else:
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.7, 0.32)


        mnv.AddHistoTitle("Target %s: %s"%(trueZ, step), 0.04, 1)


        mc_hist2.GetYaxis().SetTitleOffset(0.96)


        legend = TLegend(0.55,0.5,0.80,0.89)
        legend = TLegend(0.20,0.55,0.55,0.89)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.035)
        legend.SetTextFont(42)
        legend.AddEntry(ratio2, " Target 2", "lepl")
        legend.AddEntry(ratio3, " Target 3", "lep")
        if targetZ == "82":
            legend.AddEntry(ratio4, " Target 4 ", "lep")
        legend.AddEntry(ratio5, " Target 5", "lep")

        legend.Draw()
        gStyle.SetErrorX(0)
        
        canvas1.SetLogx(False)
        if var == "x":
            canvas1.SetLogx()
        canvas1.Modified()
        canvas1.Print("DataMCratio_%s_Daisy_z%02s_%s_%s_comparison.png"%(step, targetZ, var, plist))

    
data_hist2.SetDirectory(0)
mc_hist2.SetDirectory(0)

print("DONE %s %02s"%(plist, targetZ))
