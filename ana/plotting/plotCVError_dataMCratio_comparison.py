import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLine
from ROOT import TLegend

ROOT.gROOT.SetBatch(True)

targetID = sys.argv[1] 
targetZ = sys.argv[2]
scale = sys.argv[3]

infile = ROOT.TFile("/minerva/data/users/anezkak/05-11-2023_v4/1D/combined/EventSelection/EventSelection_minervame5A6A6B6C6D6E6F6G6H6I6J_t%s_z%02s_sys.root"%(targetID, targetZ))
infile_v1 = ROOT.TFile("/minerva/data/users/anezkak/05-11-2023_v1/1D/combined/EventSelection/EventSelection_minervame5A6A6B6C6D6E6F6G6H6I6J_t%s_z%02s_sys.root"%(targetID, targetZ))
infile_lowq2 = ROOT.TFile("/minerva/data/users/anezkak/04-26-2023_lowQ2sup/1D/combined/EventSelection/EventSelection_minervame5A6A6B6C6D6E6F6G6H6I6J_t%s_z%02s_sys.root"%(targetID, targetZ))
infile_v431 = ROOT.TFile("/minerva/data/users/anezkak/05-11-2023_v431/1D/combined/EventSelection/EventSelection_minervame5A6A6B6C6D6E6F6G6H6I6J_t%s_z%02s_sys.root"%(targetID, targetZ))
infile_v430 = ROOT.TFile("/minerva/data/users/anezkak/05-12-2023_v430/1D/combined/EventSelection/EventSelection_minervame5A6A6B6C6D6E6F6G6H6I6J_t%s_z%02s_sys.root"%(targetID, targetZ))
infile_nonres = ROOT.TFile("/minerva/data/users/anezkak/05-15-2023/nonrespiwarp_v1/combined/EventSelection/EventSelection_minervame5A6A6B6C6D6E6F6G6H6I6J_t%s_z%02s_sys.root"%(targetID, targetZ))
infile_v430_nonres = ROOT.TFile("/minerva/data/users/anezkak/05-15-2023/nonrespiwarp_v430/combined/EventSelection/EventSelection_minervame5A6A6B6C6D6E6F6G6H6I6J_t%s_z%02s_sys.root"%(targetID, targetZ))
#if syss=="False":
#    infile= ROOT.TFile(str(dirpwd)+"EventSelection_%s_t%s_z%02s_nosys.root"%(plist, targetID, targetZ))

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

if targetZ == "99":
    trueZ = "Tracker"
    mat = "CH"

vars = ["Enu", "x", "pZmu1D", "pTmu1D", "ThetamuDeg"]#, "ThetamuDeg"]

for var in vars:

    mc_hist = infile.Get("selected_mc_reco_%s"%var)
    data_hist = infile.Get("selected_data_reco_%s"%var)

    mc_hist_v1 = infile_v1.Get("selected_mc_reco_%s"%var)
    data_hist_v1 = infile_v1.Get("selected_data_reco_%s"%var)
    
    mc_hist_lowq2 = infile_lowq2.Get("selected_mc_reco_%s"%var)
    data_hist_lowq2 = infile_lowq2.Get("selected_data_reco_%s"%var)

    mc_hist_v431 = infile_v431.Get("selected_mc_reco_%s"%var)
    data_hist_v431 = infile_v431.Get("selected_data_reco_%s"%var)

    mc_hist_v430 = infile_v430.Get("selected_mc_reco_%s"%var)
    data_hist_v430 = infile_v430.Get("selected_data_reco_%s"%var)

    mc_hist_nonres = infile_nonres.Get("selected_mc_reco_%s"%var)
    data_hist_nonres = infile_nonres.Get("selected_data_reco_%s"%var)

    mc_hist_v430_nonres = infile_v430_nonres.Get("selected_mc_reco_%s"%var)
    data_hist_v430_nonres = infile_v430_nonres.Get("selected_data_reco_%s"%var)

    #############################################################################
    ###################### Pop vertical error band ##############################
    ###################### Pop response error band! #############################
    #############################################################################

    mc_hist.PopVertErrorBand("response_em")
    mc_hist.PopVertErrorBand("response_high_proton")
    mc_hist.PopVertErrorBand("response_low_proton")
    mc_hist.PopVertErrorBand("response_meson")
    mc_hist.PopVertErrorBand("response_mid_proton")
    mc_hist.PopVertErrorBand("response_other")

    mc_hist_v1.PopVertErrorBand("response_em")
    mc_hist_v1.PopVertErrorBand("response_high_proton")
    mc_hist_v1.PopVertErrorBand("response_low_proton")
    mc_hist_v1.PopVertErrorBand("response_meson")
    mc_hist_v1.PopVertErrorBand("response_mid_proton")
    mc_hist_v1.PopVertErrorBand("response_other")

    mc_hist_lowq2.PopVertErrorBand("response_em")
    mc_hist_lowq2.PopVertErrorBand("response_high_proton")
    mc_hist_lowq2.PopVertErrorBand("response_low_proton")
    mc_hist_lowq2.PopVertErrorBand("response_meson")
    mc_hist_lowq2.PopVertErrorBand("response_mid_proton")
    mc_hist_lowq2.PopVertErrorBand("response_other")

    mc_hist_v431.PopVertErrorBand("response_em")
    mc_hist_v431.PopVertErrorBand("response_high_proton")
    mc_hist_v431.PopVertErrorBand("response_low_proton")
    mc_hist_v431.PopVertErrorBand("response_meson")
    mc_hist_v431.PopVertErrorBand("response_mid_proton")
    mc_hist_v431.PopVertErrorBand("response_other")

    mc_hist_v430.PopVertErrorBand("response_em")
    mc_hist_v430.PopVertErrorBand("response_high_proton")
    mc_hist_v430.PopVertErrorBand("response_low_proton")
    mc_hist_v430.PopVertErrorBand("response_meson")
    mc_hist_v430.PopVertErrorBand("response_mid_proton")
    mc_hist_v430.PopVertErrorBand("response_other")

    mc_hist_nonres.PopVertErrorBand("response_em")
    mc_hist_nonres.PopVertErrorBand("response_high_proton")
    mc_hist_nonres.PopVertErrorBand("response_low_proton")
    mc_hist_nonres.PopVertErrorBand("response_meson")
    mc_hist_nonres.PopVertErrorBand("response_mid_proton")
    mc_hist_nonres.PopVertErrorBand("response_other")

    mc_hist_v430_nonres.PopVertErrorBand("response_em")
    mc_hist_v430_nonres.PopVertErrorBand("response_high_proton")
    mc_hist_v430_nonres.PopVertErrorBand("response_low_proton")
    mc_hist_v430_nonres.PopVertErrorBand("response_meson")
    mc_hist_v430_nonres.PopVertErrorBand("response_mid_proton")
    mc_hist_v430_nonres.PopVertErrorBand("response_other")





    if var == "Enu":
        mc_hist.GetXaxis().SetTitle("Reconstructed Antineutrino Energy (GeV)")
        mc_hist.GetYaxis().SetTitle("Events/GeV")

    if var == "x":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")

    if var == "ThetamuDeg":
        mc_hist.GetXaxis().SetTitle("Reconstructed Muon Angle (Deg)")
        mc_hist.GetYaxis().SetTitle("Events/Deg")
    
    if var == "pTmu1D":
        mc_hist.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
        mc_hist.GetYaxis().SetTitle("Events/(GeV/c)")

    if var == "pZmu1D":
        mc_hist.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)")
        mc_hist.GetYaxis().SetTitle("Events/(GeV/c)")

    '''
    if myvariable == "xfine":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")


    if myvariable == "xBrian":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")

    if myvariable == "x09":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")
    ''' 
    mc_hist.GetXaxis().CenterTitle()
    mc_hist.GetYaxis().CenterTitle()
    mc_hist.Scale(mcScale)

    data_hist_stat =  data_hist.GetCVHistoWithStatError() # stat error
    mc_hist_stat = mc_hist.GetCVHistoWithStatError() 
    mc_hist_total = mc_hist.GetCVHistoWithError()

    data_hist_stat_v1 =  data_hist_v1.GetCVHistoWithStatError() # stat error
    mc_hist_stat_v1 = mc_hist_v1.GetCVHistoWithStatError() 
    mc_hist_total_v1 = mc_hist_v1.GetCVHistoWithError()


    data_hist_stat_lowq2 =  data_hist_lowq2.GetCVHistoWithStatError() # stat error
    mc_hist_stat_lowq2 = mc_hist_lowq2.GetCVHistoWithStatError() 
    mc_hist_total_lowq2 = mc_hist_lowq2.GetCVHistoWithError()


    data_hist_stat_v431 =  data_hist_v431.GetCVHistoWithStatError() # stat error
    mc_hist_stat_v431 = mc_hist_v431.GetCVHistoWithStatError() 
    mc_hist_total_v431 = mc_hist_v431.GetCVHistoWithError()


    data_hist_stat_v430 =  data_hist_v430.GetCVHistoWithStatError() # stat error
    mc_hist_stat_v430 = mc_hist_v430.GetCVHistoWithStatError() 
    mc_hist_total_v430 = mc_hist_v430.GetCVHistoWithError()


    data_hist_stat_nonres =  data_hist_nonres.GetCVHistoWithStatError() # stat error
    mc_hist_stat_nonres = mc_hist_nonres.GetCVHistoWithStatError() 
    mc_hist_total_nonres = mc_hist_nonres.GetCVHistoWithError()


    data_hist_stat_v430_nonres =  data_hist_v430_nonres.GetCVHistoWithStatError() # stat error
    mc_hist_stat_v430_nonres = mc_hist_v430_nonres.GetCVHistoWithStatError() 
    mc_hist_total_v430_nonres = mc_hist_v430_nonres.GetCVHistoWithError()

    
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

    
    ratio = data_hist_stat.Clone()
    ratio.Divide(ratio,mc_hist_stat) # stat
    if var == "Enu":
        ratio.GetXaxis().SetTitle("Reconstructed Antineutrino Energy (GeV)")
    if var == "x":
        ratio.GetXaxis().SetTitle("Reconstructed Bjorken x")
    if var == "pTmu1D":
        ratio.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
    if var == "pZmu1D":
        ratio.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)")
    if var == "ThetamuDeg":
        ratio.GetXaxis().SetTitle("Reconstructed Muon Angle (Deg)")
        
    ratio.GetYaxis().SetTitle("Data/MC")
    ratio.GetYaxis().CenterTitle()
    ratio.GetXaxis().CenterTitle()
    ratio.Draw("X0")
    if var == "Enu":
        ratio.SetMaximum(2.0)
        ratio.SetMinimum(0.0)
    else:
        ratio.SetMaximum(2)
        ratio.SetMinimum(0.0)
    ratio.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions
    # Copied from MINERvA-101-Cross-Section/backgroundStack.py
    # same as in void MnvPlotter::DrawDataMCRatio() in MnvPlotter
    # systematic error centered at y = 1

    # line at 1
    line = TLine(ratio.GetXaxis().GetXmin(), 1, ratio.GetXaxis().GetXmax(), 1)
    line.SetLineColor(46)
    line.SetLineWidth(2)
    line.SetLineStyle(9)
    line.Draw()
    
    ratio.Draw("X0 SAME E1")

    ratio_v1 = data_hist_stat_v1.Clone()
    ratio_v1.Divide(ratio_v1,mc_hist_total_v1) # stat

    v1_colour = 221
    ratio_v1.SetLineColor(v1_colour)
    ratio_v1.SetMarkerColor(v1_colour)
    ratio_v1.SetMarkerStyle(29)
    ratio_v1.SetMarkerSize(2.2)

    ratio_lowq2 = data_hist_stat_lowq2.Clone()
    ratio_lowq2.Divide(ratio_lowq2,mc_hist_total_lowq2) # stat

    lowq2_colour = 65
    ratio_lowq2.SetLineColor(lowq2_colour)
    ratio_lowq2.SetMarkerColor(lowq2_colour)
    ratio_lowq2.SetMarkerStyle(22)
    ratio_lowq2.SetMarkerSize(1.5)

    ratio_v431 = data_hist_stat_v431.Clone()
    ratio_v431.Divide(ratio_v431,mc_hist_total_v431) # stat

    v431_colour = 96
    ratio_v431.SetLineColor(v431_colour)
    ratio_v431.SetMarkerColor(v431_colour)
    ratio_v431.SetMarkerStyle(21)
    ratio_v431.SetMarkerSize(1.2)

    ratio_v430 = data_hist_stat_v430.Clone()
    ratio_v430.Divide(ratio_v430,mc_hist_total_v430) # stat

    v430_colour = 210
    ratio_v430.SetLineColor(v430_colour)
    ratio_v430.SetMarkerColor(v430_colour)
    ratio_v430.SetMarkerStyle(33)
    ratio_v430.SetMarkerSize(1.2)

    ratio_nonres = data_hist_stat_nonres.Clone()
    ratio_nonres.Divide(ratio_nonres,mc_hist_total_nonres) # stat

    nonres_colour = 60
    ratio_nonres.SetLineColor(nonres_colour)
    ratio_nonres.SetMarkerColor(nonres_colour)
    ratio_nonres.SetMarkerStyle(24)
    ratio_nonres.SetMarkerSize(1.2)

    ratio_v430_nonres = data_hist_stat_v430_nonres.Clone()
    ratio_v430_nonres.Divide(ratio_v430_nonres,mc_hist_total_v430_nonres) # stat

    v430_nonres_colour = 14
    ratio_v430_nonres.SetLineColor(v430_nonres_colour)
    ratio_v430_nonres.SetMarkerColor(v430_nonres_colour)
    ratio_v430_nonres.SetMarkerStyle(25)
    ratio_v430_nonres.SetMarkerSize(1.2)

    gStyle.SetErrorX(0.0001)

    ratio_v1.Draw("X0 SAME E1")
    
    ratio_lowq2.Draw("X0 SAME E1")

    ratio_v431.Draw("X0 SAME E1")

    ratio_v430.Draw("X0 SAME E1")

    ratio_nonres.Draw("X0 SAME E1")

    legend = TLegend(0.18,0.16,0.9,0.36)
    legend.SetNColumns(3)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(ratio_v1, " v1", "lep")
    legend.AddEntry(ratio, " v4", "lep")
    legend.AddEntry(ratio_v430, " v4.3.0", "lep")
    legend.AddEntry(ratio_v431, " v4.3.1", "lep")
    legend.AddEntry(ratio_lowq2, " v4+lowQ2 sup.", "lep")
    legend.AddEntry(ratio_nonres, " v1 - nonrespi", "lep")
    #legend.AddEntry(ratio_v430_nonres, " v4.3.0 - nonrespi", "lep")

    legend.SetTextFont(42)
    legend.Draw()
    
    if targetZ == "99":
        mnv.AddHistoTitle("%s"%(trueZ), 0.05, 1)
    else:
        mnv.AddHistoTitle("Target %s %s"%(targetID, trueZ), 0.05, 1)
    #mnv.AddPOTNormBox(dataPOT,mcPOT, 0.7, 0.32)

     #################################################################################
    ########################## Chi2 calculation #####################################
    #################################################################################
    print("#######################################  Target " + str(targetID) + " " + str(targetZ) + " ###############################################" )
    print("############################### Chi2/ndf calculation " + str(var) +  " #####################################" )

    #chi2_v1 = mnv.Chi2DataMC(data_hist_v1, mc_hist_v1, mcScale, True, False, False)
    # data error, no statistical error from model
    chi2_v1 = mnv.Chi2DataMC(data_hist_v1, mc_hist_v1, mcScale, False, False, True)
    # mc stat error + data err
    '''
    Chi2DataMC( dataHist, mcHist, mcScale, useDataErrorMatrix, useOnlyShapeErrors, useModelStat)
    If your data distribution has ONLY STATISTICAL uncertainties, useDataErrorMatrix should be FALSE
    If you are using POT-normalization, then useOnlyShapeErrors should be FALSE
    '''
    chi2_v1_perndf = chi2_v1/data_hist_v1.GetNbinsX()+1
    print("v1 histo chi2/ndf =                                  " + str(chi2_v1_perndf))

    chi2_v4 = mnv.Chi2DataMC(data_hist, mc_hist, mcScale, False, False, True)
    chi2_v4_perndf = chi2_v4/data_hist.GetNbinsX()+1
    print("v4 histo chi2/ndf =                                  " + str(chi2_v4_perndf))

    chi2_lowq2 = mnv.Chi2DataMC(data_hist_lowq2, mc_hist_lowq2, mcScale, False, False, True)
    chi2_lowq2_perndf = chi2_lowq2/data_hist_lowq2.GetNbinsX()+1
    print("v4 + lowq2 sup. histo chi2/ndf =                     " + str(chi2_lowq2_perndf))

    chi2_v430 = mnv.Chi2DataMC(data_hist_v430, mc_hist_v430, mcScale, False, False, True)
    chi2_v430_perndf = chi2_v430/data_hist_v430.GetNbinsX()+1
    print("v4.3.0 histo chi2/ndf =                              " + str(chi2_v430_perndf))

    chi2_v431 = mnv.Chi2DataMC(data_hist_v431, mc_hist_v431, mcScale, False, False, True)
    chi2_v431_perndf = chi2_v431/data_hist_v431.GetNbinsX()+1
    print("v4.3.1 histo chi2/ndf =                              " + str(chi2_v431_perndf))

    chi2_nonres = mnv.Chi2DataMC(data_hist_nonres, mc_hist_nonres, mcScale, False, False, True)
    chi2_nonres_perndf = chi2_nonres/data_hist_nonres.GetNbinsX()+1
    print("v1 - nonrespi histo chi2/ndf =                              " + str(chi2_nonres_perndf))
    
    canvas1.SetLogx(False)
    if var == "x":
        canvas1.SetLogx()
    canvas1.Modified()
    canvas1.Print("DataMCratio_EventSelection_t%s_z%02s_%s.png"%(targetID, targetZ, var))



data_hist.SetDirectory(0)
mc_hist.SetDirectory(0)
print("DONE %s %02s"%(targetID, targetZ))