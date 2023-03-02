import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import TLegend
from ROOT import THStack
from ROOT import gStyle
from ROOT import gPad

ROOT.gROOT.SetBatch(True)

dirpwd = sys.argv[1]
targetID = sys.argv[2] 
targetZ = sys.argv[3]
plist = sys.argv[4]
scale = sys.argv[5]

infile= ROOT.TFile(str(dirpwd)+"EventSelection_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ))
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale = dataPOT/mcPOT
if scale == "1":
    mcScale = 1

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

vars = ["Enu", "x", "pZmu1D", "pTmu1D", "ThetamuDeg"]

for var in vars:
    total = infile.Get("selected_mc_reco_%s"%var)
    signal = infile.Get("selected_mc_reco_signal_%s"%var)
    # All background categories
    #MaterialOrTarget = infile.Get("selected_mc_sb_%s_WrongMaterialOrTarget"%var)
    DS_hist = infile.Get("selected_mc_DSplastic_%s"%var) # only stat err
    US_hist = infile.Get("selected_mc_USplastic_%s"%var)
    other_hist = infile.Get("selected_mc_other_%s"%var)
    WrongSign = infile.Get("selected_mc_WrongSign_%s"%var)
    NC = infile.Get("selected_mc_NC_%s"%var)
    NotEmu = infile.Get("selected_mc_NotEmu_%s"%var)

    data_hist = infile.Get("selected_data_reco_%s"%var)

    # Background percentage
    signal_per = round(signal.Integral()*100/total.Integral(),1)
    US_hist_per = round(US_hist.Integral()*100/total.Integral(),1)
    DS_hist_per = round(DS_hist.Integral()*100/total.Integral(),1)
    other_hist_per = round(other_hist.Integral()*100/total.Integral(),1)
    WrongSign_per = round(WrongSign.Integral()*100/total.Integral(),1)
    NC_per = round(NC.Integral()*100/total.Integral(),1)
    NotEmu_per = round(NotEmu.Integral()*100/total.Integral(),1) #100 - signal_per - US_hist_per - DS_hist_per - other_hist_per - WrongSign_per - NC_per
    print(NotEmu_per)

    total.Scale(total.GetBinWidth(1), "width")
    signal.Scale(signal.GetBinWidth(1), "width")
    DS_hist.Scale(DS_hist.GetBinWidth(1), "width")
    US_hist.Scale(US_hist.GetBinWidth(1), "width")
    other_hist.Scale(other_hist.GetBinWidth(1), "width")
    WrongSign.Scale(WrongSign.GetBinWidth(1), "width")
    NC.Scale(NC.GetBinWidth(1), "width")
    NotEmu.Scale(NotEmu.GetBinWidth(1), "width")

    total.Scale(mcScale)
    signal.Scale(mcScale)
    DS_hist.Scale(mcScale)
    US_hist.Scale(mcScale)
    other_hist.Scale(mcScale)
    WrongSign.Scale(mcScale)
    NC.Scale(mcScale)
    NotEmu.Scale(mcScale)

    signal.SetLineColor(ROOT.kGray+2)
    #signal.SetLineColor(46)
    #bkg.SetLineColor(12)
    DS_hist.SetLineColor(ROOT.kOrange-2)
    US_hist.SetLineColor(30)
    other_hist.SetLineColor(38)
    #MaterialOrTarget.SetLineColor(39)
    WrongSign.SetLineColor(ROOT.kRed-6)
    NC.SetLineColor(ROOT.kBlue-5)
    NotEmu.SetLineColor(ROOT.kGray)

    signal.SetFillColor(ROOT.kGray+2)
    #signal.SetFillColor(46)
    signal.SetFillStyle(3002)
    DS_hist.SetFillColor(ROOT.kOrange-2)
    US_hist.SetFillColor(30)
    other_hist.SetFillColor(38)
    #MaterialOrTarget.SetFillColor(39)
    WrongSign.SetFillColor(ROOT.kRed-6)
    NC.SetFillColor(ROOT.kBlue-5)
    NotEmu.SetFillColor(ROOT.kGray)

    signal.SetLineWidth(1)
    DS_hist.SetLineWidth(2)
    US_hist.SetLineWidth(2)
    other_hist.SetLineWidth(2)
    #MaterialOrTarget.SetLineWidth(2)
    WrongSign.SetLineWidth(2)
    NC.SetLineWidth(2)
    NotEmu.SetLineWidth(2)

    
    data_hist.Scale(data_hist.GetBinWidth(1), "width")
    data_hist.SetMarkerStyle(20)
    data_hist.SetMarkerSize(1)
    data_hist.SetMarkerColor(1)
    data_hist.SetLineWidth(1)
    data_hist.SetLineStyle(1)
    data_hist.SetLineColor(1)
    data_hist.SetMaximum(data_hist.GetMaximum()*1.3)
    data_hist.Draw("HIST p E1 X0") # for error bars, suppressed error bars along X
    
    
    stack = ROOT.THStack("stack","stack")
    if NotEmu.GetEntries() != 0:
        stack.Add(NotEmu)
    if NC.GetEntries() != 0:
        stack.Add(NC)
    stack.Add(WrongSign)
    stack.Add(other_hist)
    stack.Add(DS_hist)
    stack.Add(US_hist)
    stack.Add(signal)
    stack.Draw("SAME HIST")

    gPad.RedrawAxis()
    gPad.Update()
    

    '''
    signal.Draw("SAME HIST")
    if NotEmu.GetEntries() != 0:
        NotEmu.Draw("HIST SAME")

    NC.Draw("HIST SAME")
    WrongSign.Draw("HIST SAME")
    other_hist.Draw("HIST SAME")
    DS_hist.Draw("HIST SAME")
    US_hist.Draw("HIST SAME")
    #stack.Draw("SAME HIST")
    '''

    if var == "Enu":
        data_hist.GetXaxis().SetTitle("Reconstructed Antineutrino Energy (GeV)")
        data_hist.GetYaxis().SetTitle("Events/GeV")

    if var == "x":
        data_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        data_hist.GetYaxis().SetTitle("Events (norm.)")
        data_hist.SetMaximum(data_hist.GetMaximum()*1.2)
    
    if var == "pTmu1D":
        data_hist.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
        data_hist.GetYaxis().SetTitle("Events/(GeV/c)")

    if var == "pZmu1D":
        data_hist.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)")
        data_hist.GetYaxis().SetTitle("Events/(GeV/c)")   
    
    if var == "ThetamuDeg":
        data_hist.GetXaxis().SetTitle("Reconstructed Muon Angle (Deg)")
        data_hist.GetYaxis().SetTitle("Events/Deg")


    '''
    if var == "pTmu1D":
        stack.GetXaxis().SetTitle("Reconstructed Muon p_{T}[GeV/c] ")
        stack.GetYaxis().SetTitle("Events/(GeV/c)")
        stack.SetMaximum(600)

    if var == "pZmu":
        stack.GetXaxis().SetTitle("Reconstructed Muon p_{||}[GeV/c] ")
        stack.GetYaxis().SetTitle("Events/(GeV/c)")
        stack.SetMaximum(600)
    '''

    data_hist.GetXaxis().CenterTitle()
    data_hist.GetYaxis().CenterTitle()
    data_hist.GetYaxis().SetTitleOffset(0.8)
    data_hist.GetXaxis().SetTitleOffset(1.2)

    if var == "ThetamuDeg":
        data_hist.SetMaximum(data_hist.GetMaximum()*1.2)
    if var == "x":
        data_hist.SetMaximum(data_hist.GetMaximum()*1.1)


    data_hist.Draw("SAME HIST p E1 X0") # for error bars, suppressed error bars along X

    mnv.AddHistoTitle("Target %s %s"%(targetID, trueZ), 0.05, 1)
    
    legend = TLegend(0.40,0.60,0.80,0.89)
    legend.SetTextSize(0.035)
    if var == "ThetamuDeg":
        legend = TLegend(0.45,0.62,0.68,0.89)
    if var == "x":
        legend = TLegend(0.17,0.7,0.82,0.89)
        legend.SetNColumns(2)
        legend.SetTextSize(0.028)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.AddEntry(data_hist, " Data", "ep")
    legend.AddEntry(signal, " Signal ("+ str(signal_per) +"%)", "f")
    legend.AddEntry(US_hist, " Upstream plastic ("+ str(US_hist_per) +"%)", "fl")
    legend.AddEntry(DS_hist, " Downstream plastic ("+ str(DS_hist_per) +"%)", "fl")
    legend.AddEntry(other_hist, " Other ("+ str(other_hist_per) +"%)", "fl")
    if WrongSign.GetEntries() != 0:
        legend.AddEntry(WrongSign, " Wrong sign ("+ str(WrongSign_per) +"%)", "fl")
    if NC.GetEntries() != 0:
        legend.AddEntry(NC, " Neutral current ("+ str(NC_per) +"%)", "fl")
    if NotEmu.GetEntries() != 0:
        legend.AddEntry(NotEmu, " Outside muon energy ("+ str(NotEmu_per) +"%)", "fl")
    legend.SetTextFont(42)
    legend.Draw()

    canvas1.SetLogx(False)
    if var == "x":
        canvas1.SetLogx()
    canvas1.Modified()
    canvas1.Print("EventSelection_Bkg_t%s_z%02s_%s_%s.png"%(targetID, targetZ, var, plist))

    #data_hist.SetMaximum(data_hist.GetMaximum()*10)
    #data_hist.SetMinimum(0.1)
    #canvas1.SetLogy()
    #canvas1.Print("EventSelection_BkgBreadown_t%s_z%02s_%s_%s_log.png"%(targetID, targetZ, var, plist))


print("DONE %s %s %02s"%(plist, targetID, targetZ))