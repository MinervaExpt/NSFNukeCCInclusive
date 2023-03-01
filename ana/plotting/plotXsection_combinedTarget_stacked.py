import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLegend

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

vars = ["Enu", "x", "pZmu1D", "pTmu", "ThetamuDeg"]

#steps = ['unfolded','unfolded_effCorrected', 'crossSection', 'crossSection_total']
steps = ['total_unfolded_effCorrected','crossSection', 'crossSection_total']

for step in steps:

    for var in vars:

        data_hist = infile.Get("%s_data_%s"%(step, var))

        if step == "unfolded":
            mc_hist = infile.Get("total_efficiency_numerator_%s"%(var))
        elif step == "total_unfolded_effCorrected":
            mc_hist = infile.Get("total_simEventRate_%s"%(var))
            mc_hist_QE = infile.Get("total_simEventRate_QE_%s"%(var))
            mc_hist_RES = infile.Get("total_simEventRate_RES_%s"%(var))
            mc_hist_DIS = infile.Get("total_simEventRate_DIS_%s"%(var))
            mc_hist_Other = infile.Get("total_simEventRate_Other_%s"%(var))
            mc_hist_2p2h = infile.Get("total_simEventRate_2p2h_%s"%(var))
        else:
            mc_hist = infile.Get("simEventRate_%s_mc_%s"%(step, var))
            mc_hist_QE = infile.Get("simEventRate_QE_%s_mc_%s"%(step, var))
            mc_hist_RES = infile.Get("simEventRate_RES_%s_mc_%s"%(step, var))
            mc_hist_DIS = infile.Get("simEventRate_DIS_%s_mc_%s"%(step, var))
            mc_hist_Other = infile.Get("simEventRate_Other_%s_mc_%s"%(step, var))
            mc_hist_2p2h = infile.Get("simEventRate_2p2h_%s_mc_%s"%(step, var))
    

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

        if var == "pTmu":
            mc_hist.GetXaxis().SetTitle("Muon p_{T} (GeV/c)")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")
       
        if var == "pZmu1D":
            mc_hist.GetXaxis().SetTitle("Muon p_{Z} (GeV/c)")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")

        if var == "ThetamuDeg":
            mc_hist.GetXaxis().SetTitle("Muon Angle (Deg)")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")

        if step == "total_unfolded_effCorrected":
            mc_hist.Scale(1/1E3)
            mc_hist_QE.Scale(1/1E3)
            mc_hist_RES.Scale(1/1E3)
            mc_hist_DIS.Scale(1/1E3)
            mc_hist_Other.Scale(1/1E3)
            mc_hist_2p2h.Scale(1/1E3)
            data_hist.Scale(1/1E3)

        if step == "crossSection":
            mc_hist.Scale(1E39)
            mc_hist_QE.Scale(1E39)
            mc_hist_RES.Scale(1E39)
            mc_hist_DIS.Scale(1E39)
            mc_hist_Other.Scale(1E39)
            mc_hist_2p2h.Scale(1E39)
            data_hist.Scale(1E39)
        
        if step == "crossSection_total":
            mc_hist.Scale(1E38)
            mc_hist_QE.Scale(1E38)
            mc_hist_RES.Scale(1E38)
            mc_hist_DIS.Scale(1E38)
            mc_hist_Other.Scale(1E38)
            mc_hist_2p2h.Scale(1E38)
            data_hist.Scale(1E38)

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
            mc_hist_QE.Scale(mcScale)
            mc_hist_RES.Scale(mcScale)
            mc_hist_DIS.Scale(mcScale)
            mc_hist_Other.Scale(mcScale)
            mc_hist_2p2h.Scale(mcScale)
        else:
            mc_hist.Scale(1)
            mc_hist_QE.Scale(1)
            mc_hist_RES.Scale(1)
            mc_hist_DIS.Scale(1)
            mc_hist_Other.Scale(1)
            mc_hist_2p2h.Scale(1)

        #if var == "Enu":
        #    mc_hist.GetXaxis().SetRangeUser(2, 20)
        
        mc_hist.Scale(mc_hist.GetNormBinWidth(), "width")
        mc_hist_QE.Scale(mc_hist_QE.GetNormBinWidth(), "width")
        mc_hist_RES.Scale(mc_hist_RES.GetNormBinWidth(), "width")
        mc_hist_DIS.Scale(mc_hist_DIS.GetNormBinWidth(), "width")
        mc_hist_Other.Scale(mc_hist_Other.GetNormBinWidth(), "width")
        mc_hist_2p2h.Scale(mc_hist_2p2h.GetNormBinWidth(), "width")
        data_hist.Scale(mc_hist.GetNormBinWidth(), "width")
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
            mc_hist.GetYaxis().SetRangeUser(0, data_hist.GetMaximum()*1.3)
        elif step == "total_unfolded_effCorrected":
            mc_hist.SetMaximum(data_hist.GetMaximum()*1.6)
            if var == "x":
                mc_hist.SetMaximum(data_hist.GetMaximum()*2.0)
                if targetZ == "82":
                    mc_hist.SetMaximum(data_hist.GetMaximum()*2.2)    
        else:
        elif var == "x":
            mc_hist.SetMaximum(data_hist.GetMaximum()*1.3)
            if step == "total_unfolded_effCorrected":
                mc_hist.SetMaximum(data_hist.GetMaximum()*2.0)
                if targetZ == "82":
                    mc_hist.SetMaximum(data_hist.GetMaximum()*2.2)
        else:
            mc_hist.SetMaximum(data_hist.GetMaximum()*1.4)

        # Int channels
        mc_hist_QE.SetLineWidth(3)
        mc_hist_QE.SetLineColor(38)
        mc_hist_QE.SetFillColor(38)
        mc_hist_RES.SetLineColor(30)
        mc_hist_RES.SetFillColor(30)
        mc_hist_DIS.SetLineColor(40)
        mc_hist_DIS.SetFillColor(40)
        mc_hist_Other.SetLineColor(14)
        mc_hist_Other.SetFillColor(14)
        mc_hist_2p2h.SetLineColor(41)
        mc_hist_2p2h.SetFillColor(41)

        #gStyle.SetErrorX()
        mc_hist.Draw("HIST")
        #gStyle.SetErrorX(0.5)
        #mc_hist_stat.Draw("E2 SAME")
        #gStyle.SetErrorX(0)
        #mc_hist.Draw("HIST SAME")

        stack = ROOT.THStack("stack","stack")
        stack.Add(mc_hist_QE)
        stack.Add(mc_hist_RES)
        stack.Add(mc_hist_2p2h)
        stack.Add(mc_hist_DIS)
        stack.Add(mc_hist_Other)
        stack.Draw("SAME HIST")

        #mc_hist_QE.Draw("HIST SAME")
        #mc_hist_RES.Draw("HIST SAME")
        #mc_hist_DIS.Draw("HIST SAME")
        #mc_hist_Other.Draw("HIST SAME")
        #mc_hist_2p2h.Draw("HIST SAME")
        data_hist.Draw("SAME E1 X0")
        #data_hist_stat.Draw("E2 SAME")
        data_hist_total.Draw("E1 SAME X0")


        ''' void MnvPlotter::DrawDataMCWithErrorBand(
        const MnvH1D* dataHist,
        const MnvH1D* mcHist,
        const Double_t mcScale /*= 1.0*/,
        const std::string& legPos /*= "L"*/,
        const bool useHistTitles /*=false*/,
        const MnvH1D* bkgdHist /*= NULL*/,
        const MnvH1D* dataBkgdHist /*= NULL*/,
        const bool covAreaNormalize/*= false, Area Normalize considerations for covariance matrix*/,
        const bool statPlusSys /* = false */,
        const bool isSmooth /* = false */)
        '''

        #chi2 = mnv.Chi2DataMC(data_hist, mc_hist,mcScale,False,True )
        #ndf = mc_hist.GetNbinsX()
        #print(chi2)
        if step == "unfolded":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
        elif step == "total_unfolded_effCorrected":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
        elif step == "crossSection_total":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.7, 0.82)
        else:
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)

        '''
        Chi2DataMC( dataHist, mcHist, ndf, mcScale, useDataErrorMatrix, useOnlyShapeErrors)
        If your data distribution has ONLY STATISTICAL uncertainties, useDataErrorMatrix should be FALSE
        If you are using POT-normalization, then useOnlyShapeErrors should be FALSE
        '''
        if step == 'crossSection':
            mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)
        
        elif step == 'crossSection_total':
            mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TR", 0.035, 0.0, True, False)
        
        if targetZ == "99":
            mnv.AddHistoTitle("%s"%(trueZ), 0.05, 1)
        else:
            mnv.AddHistoTitle("Target %s %s: %s"%(targetID, trueZ, step), 0.04, 1)

        mc_hist.GetYaxis().SetTitleOffset(0.96)

        if step == "crossSection_total":
            legend = TLegend(0.20,0.45,0.50,0.89)
        elif var == "x" or "ThetamuDeg":
            legend = TLegend(0.70,0.55,0.85,0.89)
        else:
            legend = TLegend(0.55,0.55,0.80,0.89)
        
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.035)
        #legend.AddEntry(mc_hist_stat, " Simulation", "fl")
        legend.AddEntry(mc_hist_QE, " QE", "fl")
        legend.AddEntry(mc_hist_RES, " RES", "fl")
        legend.AddEntry(mc_hist_2p2h, " 2p2h", "fl")
        legend.AddEntry(mc_hist_DIS, " DIS", "fl")
        legend.AddEntry(mc_hist_Other, " Other", "fl")
        legend.AddEntry(data_hist, " Data", "lep")
        #legend.SetTextFont(42)
        legend.Draw()
        canvas1.RedrawAxis()

        canvas1.SetLogx(False)
        if var == "x":
            canvas1.SetLogx()
        canvas1.Modified()
        if targetZ == "99":
            canvas1.Print("Stacked_%s_t%s_z%02s_%s_%s.png"%(step, targetID, targetZ, var, plist))
        else:
            canvas1.Print("Stacked_%s_Daisy_t%s_z%02s_%s_%s.png"%(step, targetID, targetZ, var, plist))

data_hist.SetDirectory(0)
mc_hist.SetDirectory(0)

print("DONE %s %s %02s"%(plist, targetID, targetZ))
