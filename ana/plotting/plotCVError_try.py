import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle

ROOT.gROOT.SetBatch(True)

dirpwd = sys.argv[1]
targetID = sys.argv[2] 
targetZ = sys.argv[3]
plist = sys.argv[4]
scale = sys.argv[5]

infile= ROOT.TFile(str(dirpwd)+"/EventSelection_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ))
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

#vars = ["Enu", "x", "pZmu1D", "pTmu1D", "ThetamuDeg"]#, "ThetamuDeg"]
vars = ["Enu", "Ehad"]#, "ThetamuDeg"]


for var in vars:

    mc_hist = infile.Get("selected_mc_reco_%s"%var)
    #data_hist = infile.Get("selected_data_reco_%s"%var)
    #bkg_mc = infile.Get("selected_mc_reco_bkg_%s"%var)

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

    if var == "Ehad":
        mc_hist.GetXaxis().SetTitle("Reconstructed Recoil (GeV)")
        mc_hist.GetYaxis().SetTitle("Events/(GeV)")

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
    #print(data_hist.GetBinContent(data_hist.GetNbinsX())) # last bin
    #print(data_hist.GetBinContent(data_hist.GetNbinsX()+1)) # overflow
    #mc_hist.GetXaxis().SetRangeUser(0,0.001)

    #mnv.ApplyStyle(1)

    mnv.DrawMCWithErrorBand(mc_hist, mcScale)
    #mnv.DrawDataMCWithErrorBand(data_hist, mc_hist, mcScale, "TR", False, None, None, False, True, False)

    #mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)

    '''
    Chi2DataMC( dataHist, mcHist, ndf, mcScale, useDataErrorMatrix, useOnlyShapeErrors)
    If your data distribution has ONLY STATISTICAL uncertainties, useDataErrorMatrix should be FALSE
    If you are using POT-normalization, then useOnlyShapeErrors should be FALSE
    '''
    if targetZ == "99":
        mnv.AddHistoTitle("%s"%(trueZ), 0.05, 1)
    else:
        mnv.AddHistoTitle("Target %s %s"%(targetID, trueZ), 0.05, 1)
    #mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)

    canvas1.SetLogx(False)
    if var == "x":
        canvas1.SetLogx()
    canvas1.Modified()
    canvas1.Print("EventSelection_t%s_z%02s_%s_%s.png"%(targetID, targetZ, var, plist))


    # ----------------------------------------------------------------------------
    # Event selection with simulated background

    #mnv.DrawDataMCWithErrorBand(data_hist, mc_hist, mcScale, "TR", False, bkg_mc, None, False, True, False)
    #mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)
    #if targetZ == "99":
        #mnv.AddHistoTitle("%s"%(trueZ), 0.05, 1)
    #else:
        #mnv.AddHistoTitle("Target %s %s"%(targetID, trueZ), 0.05, 1)
    #mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
    #canvas1.Print("EventSelection_t%s_z%02s_%s_%s_wsimbkg.png"%(targetID, targetZ, var, plist))


    # ----------------------------------------------------------------------------

    # Fractional Error Groups
    mnv2 = PlotUtils.MnvPlotter()
    mnv2.ApplyStyle(7)

    mnv2.error_summary_group_map.clear()
    '''
    mnv2.error_summary_group_map["Hadrons"].push_back("response_em")
    mnv2.error_summary_group_map["Hadrons"].push_back("response_low_proton")
    mnv2.error_summary_group_map["Hadrons"].push_back("response_mid_proton")
    mnv2.error_summary_group_map["Hadrons"].push_back("response_high_proton")
    mnv2.error_summary_group_map["Hadrons"].push_back("response_meson")
    mnv2.error_summary_group_map["Hadrons"].push_back("response_other")
    '''
    # Change colour
    mnv2.error_color_map["Particle Response"] = 65
    
    #mnv2.ApplyStyle(7)
    mnv2.DrawErrorSummary(mc_hist, "TL", True, True, 0.0, False, "",True);
    # last boolean decides whether frac or not

    keys = canvas1.GetListOfPrimitives();
    for k in keys:
        if(k.ClassName().find("TH1")!=-1):
            if var == "Enu":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum()*1.02) # to 0.35 for Enu
            if var == "x":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum()*1.05)
            if var == "pTmu1D":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum()*1.05)
            if var == "pZmu1D":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum()* 1.05)
            if var == "ThetamuDeg":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum()*1.05)

        if(k.ClassName().find("Legend")!=-1):
            if var == "Enu":
                k.SetNColumns(2)
                k.SetX1(2) #Enu
                k.SetX2(24) #Enu
                k.SetY1(0.05) #Enu
                if targetZ == "99":
                    k.SetY1(0.14) #Enu
            
            if var == "Ehad":
                k.SetNColumns(2)
                k.SetX1(0) #Enu
                k.SetX2(29) #Enu
                k.SetY1(0.15) #Enu
                if targetZ == "99":
                    k.SetY1(0.14) #Enu


            if var == "x":
                k.SetNColumns(2)
                k.SetX2(0.5) #x
                k.SetY1(0.16) #x
                if targetZ == "99":
                    k.SetY1(0.18) #Enu

            if var == "pTmu1D":
                k.SetNColumns(2)
                k.SetX2(3) 
                k.SetY1(0.17) #Enu
                if targetZ == "99":
                    k.SetY1(0.16) #Enu
                if targetZ == "82":
                    k.SetY1(0.18) #Enu
                    if targetID == "3":
                       k.SetY1(0.24) #Enu 
                    if targetID == "5":
                        k.SetY1(0.22) #Enu 
                if targetZ == "06":
                    k.SetY1(0.18) #Enu

            if var == "pZmu1D":
                k.SetNColumns(2)
                k.SetX2(20)
                if targetZ == "26":
                    k.SetY1(0.16) #Enu
                if targetZ == "26":
                    k.SetY1(0.16) #Enu
                if targetZ == "99":
                    k.SetY1(0.15) #Enu
        
            if var == "ThetamuDeg":
                k.SetNColumns(2)
                k.SetX2(17) 
                k.SetY1(0.16)
                if targetZ == "99":
                    k.SetY1(0.14) #Enu
                if targetZ == "06":
                    k.SetY1(0.15) #Enu
                if targetZ == "82":
                    k.SetY1(0.15)


    if targetZ == "99":
        mnv2.AddHistoTitle("%s"%(trueZ), 0.05, 1)
    else:
        mnv2.AddHistoTitle("Target %s %s"%(targetID, trueZ), 0.05, 1)

    canvas1.Modified()
    canvas1.Print("EventSelection_t%s_z%02s_%s_%s_FracErr.png"%(targetID, targetZ, var, plist))


#data_hist.SetDirectory(0)
mc_hist.SetDirectory(0)
print("DONE %s %s %02s"%(plist, targetID, targetZ))