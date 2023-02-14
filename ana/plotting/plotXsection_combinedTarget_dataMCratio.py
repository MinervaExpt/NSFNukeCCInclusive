import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLegend
from ROOT import TLine

ROOT.gROOT.SetBatch(True)

dirpwd = sys.argv[1]
targetID = sys.argv[2] 
targetZ = sys.argv[3]
plist = sys.argv[4]
scale = sys.argv[5]

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

if targetZ == "99":
    infile= ROOT.TFile(str(dirpwd)+"/CrossSection_t%s_z%02s_%s.root"%(targetID, targetZ, plist))
else:
    infile= ROOT.TFile(str(dirpwd)+"/CrossSection_Daisy_t%s_z%02s_%s.root"%(targetID, targetZ, plist))
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale = dataPOT/mcPOT
if scale == "1":
    mcScale = 1

vars = ["Enu", "x"]

#steps = ['unfolded','unfolded_effCorrected', 'crossSection', 'crossSection_total']
steps = ['total_unfolded_effCorrected','crossSection', 'crossSection_total']

for step in steps:

    for var in vars:

        data_hist = infile.Get("%s_data_%s"%(step, var))

        if step == "unfolded":
            mc_hist = infile.Get("total_efficiency_numerator_%s"%(var))
        elif step == "total_unfolded_effCorrected":
            mc_hist = infile.Get("total_simEventRate_%s"%(var))
        else:
            mc_hist = infile.Get("simEventRate_%s_mc_%s"%(step, var))
    

        if var == "Enu":
            mc_hist.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dE_{#bar{#nu}} (10^{-39} cm^{2}/GeV/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3}/GeV")

        if var == "x":
            mc_hist.GetXaxis().SetTitle("Bjorken x")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")


        mc_hist.GetXaxis().CenterTitle()
        mc_hist.GetYaxis().CenterTitle()
        #print(data_hist.GetBinContent(data_hist.GetNbinsX())) # last bin
        #print(data_hist.GetBinContent(data_hist.GetNbinsX()+1)) # overflow
        #mc_hist.GetXaxis().SetRangeUser(66,173)

        #mnv.ApplyStyle(1)
        if step == "unfolded":
            mc_hist.Scale(mcScale)
        elif step == "total_unfolded_effCorrected":
            mc_hist.Scale(mcScale)
        else:
            mc_hist.Scale(1)

        if var == "Enu":
            mc_hist.GetXaxis().SetRangeUser(2, 20)

        data_hist_stat =  data_hist.GetCVHistoWithStatError() # stat error
        data_hist_total = data_hist.GetCVHistoWithError() # total error
        data_hist_sys = data_hist.GetCVHistoWithError(False) # sys error (bool is include stat)
        mc_hist_stat = mc_hist.GetCVHistoWithStatError() 
        
        # MC
        mc_hist.SetLineWidth(3)
        mc_hist.SetLineColor(2)
        # MC error
        mc_hist_stat.SetFillColor(ROOT.kRed-10)
        mc_hist_stat.SetFillStyle(1001)
        mc_hist_stat.SetMarkerStyle(0)
        mc_hist_stat.SetLineWidth(3)
        mc_hist_stat.SetLineColor(2)

        mc_hist.SetLineColor(ROOT.kRed)
        mc_hist.SetLineWidth(2)
        if step == "crossSection_total":
            mc_hist.GetYaxis().SetRangeUser(0, 10)
        else:
            mc_hist.SetMaximum(data_hist.GetMaximum()*1.5)

        
        ratio = data_hist_stat.Clone()
        ratio.Divide(ratio,mc_hist_stat) # stat
        if var == "Enu":
            ratio.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
        if var == "x":
            ratio.GetXaxis().SetTitle("Bjorken x")
        ratio.GetYaxis().SetTitle("Data/MC")
        ratio.GetYaxis().CenterTitle()
        ratio.GetXaxis().CenterTitle()
        ratio.Draw("X0")
        if var == "Enu":
            ratio.SetMaximum(2.0)
            ratio.SetMinimum(0.0)
        else:
            ratio.SetMaximum(1.5)
            ratio.SetMinimum(0.5)
        ratio.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions
        # Copied from MINERvA-101-Cross-Section/backgroundStack.py
        # same as in void MnvPlotter::DrawDataMCRatio() in MnvPlotter
        # systematic error centered at y = 1
        ratio_tot = data_hist_total.Clone()
        ratio_tot.Divide(ratio_tot, mc_hist_stat)

        # line at 1
        line = TLine(ratio.GetXaxis().GetXmin(), 1, ratio.GetXaxis().GetXmax(), 1)
        line.SetLineColor(46)
        line.SetLineWidth(2)
        line.SetLineStyle(9)
        line.Draw()
        
        ratio.Draw("X0 SAME E1")
        ratio_tot.Draw("E1 SAME")

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

        if targetZ == "99":
            mnv.AddHistoTitle("%s"%(trueZ), 0.05, 1)
        else:
            mnv.AddHistoTitle("Target %s %s: %s"%(targetID, trueZ, step), 0.04, 1)


        mc_hist.GetYaxis().SetTitleOffset(0.96)

        canvas1.Modified()
        if targetZ == "99":
            canvas1.Print("DataMCratio_%s_t%s_z%02s_%s_%s.png"%(step, targetID, targetZ, var, plist))
        else:
            canvas1.Print("DataMCratio_%s_Daisy_t%s_z%02s_%s_%s.png"%(step, targetID, targetZ, var, plist))

    
data_hist.SetDirectory(0)
mc_hist.SetDirectory(0)

print("DONE %s %s %02s"%(plist, targetID, targetZ))
