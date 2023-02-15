import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLegend
from ROOT import gPad

ROOT.gROOT.SetBatch(True)

dirpwd = sys.argv[1]
plist = sys.argv[2]
scale = sys.argv[3]

targetID = 99 
targetZ = 99


infile= ROOT.TFile(str(dirpwd)+"/CrossSection_Daisy_t%s_z%02s_%s.root"%(targetID, targetZ, plist))
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()
material = 'tracker_daisy'

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale = dataPOT/mcPOT
if scale == "1":
    mcScale = 1

mat = None
trueZ = None

if targetZ == "99":
    trueZ = "Tracker"
    mat = "CH"

vars = ["Enu", "x", "pTmu", "pZmu"]

steps = [ 'crossSection', 'crossSection_total']
#steps = ['crossSection', 'crossSection_total']

for step in steps:

    for var in vars:

        data_hist_C = infile.Get("%s_carbon_data_%s"%(step, var))
        data_hist_Fe = infile.Get("%s_iron_data_%s"%(step, var))
        data_hist_Pb = infile.Get("%s_lead_data_%s"%(step, var))

        mc_hist_C = infile.Get("simEventRate_%s_carbon_mc_%s"%(step, var))
        mc_hist_Fe = infile.Get("simEventRate_%s_iron_mc_%s"%(step, var))
        mc_hist_Pb = infile.Get("simEventRate_%s_lead_mc_%s"%(step, var))

        if var == "Enu":
            mc_hist_C.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
            if step == "crossSection":
                mc_hist_C.GetYaxis().SetTitle("d#sigma/dE_{#bar{#nu}} (10^{-39} cm^{2}/GeV/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist_C.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events/GeV")

        if var == "x":
            mc_hist_C.GetXaxis().SetTitle("Bjorken x")
            if step == "crossSection":
                mc_hist_C.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist_C.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist_C.GetYaxis().SetTitle("Events (norm.)")

        if var == "pTmu":
            mc_hist.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")
       
        if var == "pZmu":
            mc_hist.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")

        if step == "crossSection":
            mc_hist_C.Scale(1E39)
            mc_hist_Fe.Scale(1E39)
            mc_hist_Pb.Scale(1E39)
            data_hist_C.Scale(1E39)
            data_hist_Fe.Scale(1E39)
            data_hist_Pb.Scale(1E39)
        
        if step == "crossSection_total":
            mc_hist_C.Scale(1E38)
            mc_hist_Fe.Scale(1E38)
            mc_hist_Pb.Scale(1E38)
            data_hist_C.Scale(1E38)
            data_hist_Fe.Scale(1E38)
            data_hist_Pb.Scale(1E38)

        if var == "Enu":
            mc_hist_C.GetXaxis().SetRangeUser(2, 20)
            mc_hist_Fe.GetXaxis().SetRangeUser(2, 20)
            mc_hist_Pb.GetXaxis().SetRangeUser(2, 20)
            data_hist_C.GetXaxis().SetRangeUser(2, 20)
            data_hist_Fe.GetXaxis().SetRangeUser(2, 20)
            data_hist_Pb.GetXaxis().SetRangeUser(2, 20)

        mc_hist_C.GetXaxis().CenterTitle()
        mc_hist_C.GetYaxis().CenterTitle()
        #print(data_hist.GetBinContent(data_hist.GetNbinsX())) # last bin
        #print(data_hist.GetBinContent(data_hist.GetNbinsX()+1)) # overflow
        #mc_hist.GetXaxis().SetRangeUser(66,173)

        #mnv.ApplyStyle(1)
  
        mc_hist_C.Scale(mc_hist_C.GetNormBinWidth(), "width")
        mc_hist_Fe.Scale(mc_hist_Fe.GetNormBinWidth(), "width")
        mc_hist_Pb.Scale(mc_hist_Pb.GetNormBinWidth(), "width")

        data_hist_C.Scale(mc_hist_C.GetNormBinWidth(), "width")
        data_hist_Fe.Scale(mc_hist_Fe.GetNormBinWidth(), "width")
        data_hist_Pb.Scale(mc_hist_Pb.GetNormBinWidth(), "width")

        data_hist_stat_C =  data_hist_C.GetCVHistoWithStatError() # stat error
        data_hist_total_C = data_hist_C.GetCVHistoWithError() # total error
        data_hist_sys_C = data_hist_C.GetCVHistoWithError(False) # sys error (bool is include stat)
        mc_hist_stat_C = mc_hist_C.GetCVHistoWithStatError() 

        data_hist_stat_Fe =  data_hist_Fe.GetCVHistoWithStatError() # stat error
        data_hist_total_Fe = data_hist_Fe.GetCVHistoWithError() # total error
        data_hist_sys_Fe = data_hist_Fe.GetCVHistoWithError(False) # sys error (bool is include stat)
        mc_hist_stat_Fe = mc_hist_Fe.GetCVHistoWithStatError() 

        data_hist_stat_Pb =  data_hist_Pb.GetCVHistoWithStatError() # stat error
        data_hist_total_Pb = data_hist_Pb.GetCVHistoWithError() # total error
        data_hist_sys_Pb = data_hist_Pb.GetCVHistoWithError(False) # sys error (bool is include stat)
        mc_hist_stat_Pb = mc_hist_Pb.GetCVHistoWithStatError() 
        
        # MC
        mc_hist_C.SetLineWidth(3)
        mc_hist_C.SetLineColor(46)
        mc_hist_Fe.SetLineWidth(3)
        mc_hist_Fe.SetLineColor(38)
        mc_hist_Pb.SetLineWidth(3)
        mc_hist_Pb.SetLineColor(30)

        mc_hist_stat_C.SetFillColor(ROOT.kRed-10)
        mc_hist_stat_C.SetFillStyle(1001)
        mc_hist_stat_C.SetMarkerStyle(0)
        mc_hist_stat_C.SetLineWidth(3)
        mc_hist_stat_C.SetLineColor(46)

        if step == "crossSection_total":
            mc_hist_C.GetYaxis().SetRangeUser(0,10)
        else:
            mc_hist_C.SetMaximum(data_hist_C.GetMaximum()*1.2)

        mc_hist_stat_Fe.SetFillColor(ROOT.kBlue-10)
        mc_hist_stat_Fe.SetFillStyle(1001)
        mc_hist_stat_Fe.SetMarkerStyle(0)
        mc_hist_stat_Fe.SetLineWidth(3)
        mc_hist_stat_Fe.SetLineColor(38)

        mc_hist_stat_Pb.SetFillColor(ROOT.kGreen-10)
        mc_hist_stat_Pb.SetFillStyle(1001)
        mc_hist_stat_Pb.SetMarkerStyle(0)
        mc_hist_stat_Pb.SetLineWidth(3)
        mc_hist_stat_Pb.SetLineColor(30)

        #gStyle.SetErrorX()
        data_hist_C.SetMarkerColor(46)
        data_hist_C.SetLineColor(46)
        data_hist_total_C.SetLineColor(46)
        data_hist_total_C.SetMarkerColor(46)
        data_hist_Fe.SetMarkerColor(38)
        data_hist_Fe.SetLineColor(38)
        data_hist_total_Fe.SetLineColor(38)
        data_hist_total_Fe.SetMarkerColor(38)
        data_hist_Pb.SetMarkerColor(30)
        data_hist_Pb.SetLineColor(30)
        data_hist_total_Pb.SetLineColor(30)
        data_hist_total_Pb.SetMarkerColor(30)

        mc_hist_C.Draw("HIST")
        mc_hist_stat_C.Draw("E2 SAME")
        mc_hist_C.Draw("HIST SAME")
        data_hist_C.Draw("SAME E1 X0")
        data_hist_total_C.Draw("SAME E1 X0")
        #gStyle.SetErrorX(0.5)
        #mc_hist_stat_C.Draw("E2 SAME")
        #gStyle.SetErrorX(0)
        #mc_hist_C.Draw("HIST SAME")
        #data_hist_stat_C.Draw("SAME E1 X0")
        #data_hist_stat.Draw("E2 SAME")
        #data_hist_total_C.Draw("E1 SAME X0")

        mc_hist_Fe.Draw("HIST SAME")
        mc_hist_stat_Fe.Draw("E2 SAME")
        mc_hist_Fe.Draw("HIST SAME")
        data_hist_Fe.Draw("SAME E1 X0")
        data_hist_total_Fe.Draw("SAME E1 X0")

        mc_hist_Pb.Draw("HIST SAME")
        mc_hist_stat_Pb.Draw("E2 SAME")
        mc_hist_Pb.Draw("HIST SAME")
        data_hist_Pb.Draw("SAME E1 X0")
        data_hist_total_Pb.Draw("SAME E1 X0")

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

        '''
        Chi2DataMC( dataHist, mcHist, ndf, mcScale, useDataErrorMatrix, useOnlyShapeErrors)
        If your data distribution has ONLY STATISTICAL uncertainties, useDataErrorMatrix should be FALSE
        If you are using POT-normalization, then useOnlyShapeErrors should be FALSE
        '''
        #if step == 'crossSection':
        #    mnv.AddChi2Label(data_hist_C, mc_hist_C, mcScale, "TL", 0.035, 0.0, True, False)
        
        #if step == 'crossSection_total':
        #    mnv.AddChi2Label(data_hist_C, mc_hist_C, mcScale, "TL", 0.035, 0.0, True, False)
        
        mnv.AddHistoTitle("%s %s"%(material,step), 0.05, 1)
        if step == "crossSection_total":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.7, 0.82)
        else:
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.7, 0.45)

        mc_hist_C.GetYaxis().SetTitleOffset(0.96)

        if step == "crossSection_total":
            legend = TLegend(0.20,0.45,0.50,0.89)
        else:
            legend = TLegend(0.55,0.5,0.80,0.89)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.035)
        legend.AddEntry(mc_hist_stat_C, " Simulation C", "fl")
        legend.AddEntry(mc_hist_stat_Fe, " Simulation Fe", "fl")
        legend.AddEntry(mc_hist_stat_Pb, " Simulation Pb", "fl")

        legend.AddEntry(data_hist_C, " Data C", "lep")
        legend.AddEntry(data_hist_Fe, " Data Fe", "lep")
        legend.AddEntry(data_hist_Pb, " Data Pb", "lep")
        legend.SetTextFont(42)
        legend.Draw()

        #gPad.SetLogx(20)

        gStyle.SetErrorX(0)

        canvas1.Modified()
        canvas1.Print("%s_Daisy_t%s_z%02s_%s_%s.png"%(step, targetID, targetZ, var, plist))

        

raw_input("Done")
