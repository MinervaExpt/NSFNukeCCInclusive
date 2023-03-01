import ROOT
from ROOT import *
import os,sys
from ROOT import PlotUtils
import ctypes
import math as m

ROOT.gROOT.SetBatch(True)

def plotStackRatio(mnvPlotter, mcArray, data, mcScale, tuning, sideband, var):#, #tuning, sideband,):
    #Apply a different color for each MC category
    mcColors = MnvColors.GetColors(MnvColors.kOkabeItoDarkPalette)
    nextColor = 0

    totalMC = mcArray[0].Clone()
    totalMC.Add(mcArray[1])
    totalMC.Add(mcArray[2])
    totalMC.Add(mcArray[3])
    totalMC.Scale(mcScale)

    for hist in mcArray:
        hist.SetLineColor(mcColors[nextColor])
        hist.SetFillColor(mcColors[nextColor])
        nextColor = nextColor + 1

    dataWithStatErrors = data.GetCVHistoWithError().Clone()

    #Create a TCanvas on which to draw plots and split it into 2 panels
    overall = TCanvas("Data/MC for " + var, "plot", 800, 800)
    top = TPad("Overlay", "Overlay", 0, bottomFraction, 1, 1)
    bottom = TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction + margin)
    #Thou shalt Draw() new TPads lest they be blank!
    top.Draw()
    bottom.Draw()

    top.cd()
    stack = THStack()
    for hist in mcArray:
        hist.Scale(hist.GetNormBinWidth(), "width")
        hist.Scale(mcScale)

    stack.Add(mcArray[0])
    stack.Add(mcArray[1])
    stack.Add(mcArray[2])
    stack.Add(mcArray[3])
   
    stack.Draw("HIST")
    stack.SetMinimum(0.001)

    if var == "Enu":
        stack.GetYaxis().SetTitle("Events/GeV")
    if var == "x":
        stack.GetYaxis().SetTitle("Events (norm.)")
    if var == "planeDNN":
        stack.GetYaxis().SetTitle("Events (norm.)")
    if var == "pTmu":
        stack.GetXaxis().SetTitle("Events/(GeV/c)")
    if var == "pZmu1D":
        stack.GetXaxis().SetTitle("Events/(GeV/c)")
    if var == "ThetamuDeg":
        stack.GetXaxis().SetTitle("Events/Deg")

    stack.GetYaxis().CenterTitle()
    stack.GetYaxis().SetTitleOffset(0.9)
    dataNorm = dataWithStatErrors.Clone()
    dataNorm.Scale(data.GetNormBinWidth(), "width")
    dataNorm.Draw("X0 E1 SAME")
    dataNorm.SetMarkerSize(1)
    dataNorm.SetMarkerStyle(20)
    stack.SetMaximum(1.2*dataNorm.GetMaximum())
    #mnvPlotter.DrawDataStackedMC(data, mcArray, 1, "TR", "Data", -2, -2, 3001)
    mnvPlotter.AddHistoTitle("%s %s: Target %s %s"%(tuning, sideband, targetID, trueZ), 0.04)

    # Chi2 
    """
    ndfStat = ctypes.c_int(1)
    print(ndfStat)
    chi2Sys = mnv.Chi2DataMC(data, totalMC, ndfStat, 1 )
    ndfStat = ndfStat.value -1
    """
    ndfStat = ctypes.c_int(1)
    chi2Sys = mnv.Chi2DataMC(data, totalMC, ndfStat, 1 )
    ndfStat = 0
    for bin in range(totalMC.GetNbinsX() + 1):
        if data.GetBinContent(bin) > 0.:
            ndfStat = ndfStat+1

    ndfStat = ndfStat - 1
    print(ndfStat)


    labelStat = "#chi^{2}/ndf = %3.2f/%d = %3.2f"%(chi2Sys, ndfStat, chi2Sys/ndfStat)
    mnv.AddPlotLabel( labelStat, .35, .875, .035 )

    legend = TLegend(0.55,0.65,0.8,0.89)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(data, " Data", "ep")
    if sideband == "Upstream":
        legend.AddEntry(mcArray[2], " Upstream", "f")
        legend.AddEntry(mcArray[3], " Downstream", "f")
    if sideband == "Downstream":
        legend.AddEntry(mcArray[3], " Downstream", "f")
        legend.AddEntry(mcArray[2], " Upstream", "f")
        
    legend.AddEntry(mcArray[0], " Signal ", "f")
    legend.AddEntry(mcArray[1], " Other", "f")
    legend.SetTextFont(62)
    legend.Draw()
    
    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.4)
    
    ratio = dataWithStatErrors.Clone() # stat

    ratio.Divide(ratio,totalMC.GetCVHistoWithError()) # stat
    ratio.GetYaxis().SetTitle("Data/MC")
    ratio.GetYaxis().CenterTitle()
    ratio.GetYaxis().SetTitleSize(0.12)
    ratio.GetYaxis().SetTitleOffset(0.3)
    ratio.GetXaxis().CenterTitle()
    ratio.GetYaxis().CenterTitle()
    ratio.Draw("X0")

    ratio.SetMaximum(1.5)
    ratio.SetMinimum(0.5)
    ratio.GetXaxis().SetTitleSize(titleSize)
    ratio.GetXaxis().SetLabelSize(labelSize)
    ratio.GetYaxis().SetLabelSize(0.11)
    ratio.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions
    ratio.GetYaxis().SetTitleSize(titleSize)

    # Copied from MINERvA-101-Cross-Section/backgroundStack.py
    # same as in void MnvPlotter::DrawDataMCRatio() in MnvPlotter
    # systematic error centered at y = 1
    sysError = totalMC.GetTotalError(False, True, False) # False for stat error, True for as frac, False for covAreaNorm
    for whichBin in range(1, sysError.GetXaxis().GetNbins()+1):
        sysError.SetBinError(whichBin, max(sysError.GetBinContent(whichBin), 1e-9))
        sysError.SetBinContent(whichBin, 1)

    sysError.SetFillColor(ROOT.kRed-10)
    sysError.SetFillStyle(1001)
    sysError.SetMarkerStyle(0)
    sysError.Draw("E2 SAME")

    # line at 1
    line = TLine(ratio.GetXaxis().GetXmin(), 1, ratio.GetXaxis().GetXmax(), 1)
    line.SetLineColor(46)
    line.SetLineWidth(lineSize)
    line.SetLineStyle(9)
    line.Draw()

    ratio.Draw("X0 SAME")

    if var == "Enu":
        ratio.GetXaxis().SetTitle("Reconstructed Antineutrino Energy (GeV)")

    if var == "x":
        ratio.GetXaxis().SetTitle("Reconstructed Bjorken x")
    if var == "planeDNN":
        ratio.GetXaxis().SetTitle("Plane Number")
    if var == "pTmu":
        ratio.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
    if var == "pZmu1D":
        ratio.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)")
    if var == "ThetamuDeg":
        ratio.GetXaxis().SetTitle("Reconstructed Muon Angle (Deg)")

    
    # ratio equivalent to this:
    #mnvPlotter.DrawDataMCRatio(data, totalMC, 1, True, True, 0.7, 1.5, "Data/Total MC" )

    overall.Print("T" + targetID + "_z" + targetZ + "_" + sideband+"_"+tuning+ "_" + var+ "_" + plist + ".png")


# ---------------------------------------------------------------------------------------------------

gROOT.SetBatch() #Don't render histograms to a window.  Also gets filled areas correct.

bottomFraction = 0.15
margin = 0.116 #Tuned by hand
labelSize = 0.15
lineSize = 2
titleSize = 0.16

TH1.AddDirectory(False)

# Run like: python bkgStack+ratio.py 3 26
# To plot plots for target 3 iron
dirpwd = sys.argv[1]
outdir = sys.argv[2]
targetID = sys.argv[3] 
targetZ = sys.argv[4]
plist = sys.argv[5]

infileUntuned = TFile.Open(str(dirpwd)+"PlasticBkg_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ))
scaleFactors = TFile.Open(str(outdir)+"/Plastic_ScaleFactors_t%s_z%02s_%s.root"%(targetID,targetZ, plist))
mnv = PlotUtils.MnvPlotter()

mcPOT = infileUntuned.Get("MCPOT").GetVal()
dataPOT = infileUntuned.Get("DataPOT").GetVal()
mcScale = dataPOT/mcPOT

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

vars = ["Enu", "x", "planeDNN", "pZmu1D", "pTmu", "ThetamuDeg"]

for var in vars:
    # ----------------------------------------------------------------------------
    #                             UNTUNED HISTOGRAMS
    # ----------------------------------------------------------------------------
    # Read in untuned histograms
    # Upstream
    hists_US_data = infileUntuned.Get("selected_data_reco_US_%s"%var)
    hists_US_mat = infileUntuned.Get("US_%s_%s_%s"%(str(mat),str(trueZ), str(var)))
    hists_US_other = infileUntuned.Get("US_other_%s_%s"%(str(trueZ), str(var)))
    hists_US_regUS = infileUntuned.Get("US_regUS_%s_%s"%(str(trueZ), str(var)))
    hists_US_regDS = infileUntuned.Get("US_regDS_%s_%s"%(str(trueZ), str(var)))
    # Downstream
    hists_DS_data = infileUntuned.Get("selected_data_reco_DS_%s"%var)
    hists_DS_mat = infileUntuned.Get("DS_%s_%s_%s"%(str(mat),str(trueZ), str(var)))
    hists_DS_other = infileUntuned.Get("DS_other_%s_%s"%(str(trueZ), str(var)))
    hists_DS_regUS = infileUntuned.Get("DS_regUS_%s_%s"%(str(trueZ), str(var)))
    hists_DS_regDS = infileUntuned.Get("DS_regDS_%s_%s"%(str(trueZ), str(var)))

    print("Upstream\n")
    print("Total")
    tot = hists_US_regUS.GetEntries()+hists_US_mat.GetEntries()+hists_US_other.GetEntries()+hists_US_regDS.GetEntries()
    nofit=hists_US_mat.GetEntries()+hists_US_other.GetEntries()+hists_US_regDS.GetEntries()
    print(tot)
    print("Not fitted")
    print(nofit*100/tot)
    print(hists_US_regUS.GetEntries()*100/tot)
    print(hists_US_mat.GetEntries()*100/tot)
    print(hists_US_other.GetEntries()*100/tot)
    print(hists_US_regDS.GetEntries()*100/tot)

    print("Downstream\n")
    print("Total")
    tot = hists_DS_regDS.GetEntries()+hists_DS_mat.GetEntries()+hists_DS_other.GetEntries()+hists_DS_regUS.GetEntries()
    nofit=hists_DS_mat.GetEntries()+hists_DS_other.GetEntries()+hists_DS_regUS.GetEntries()
    print(tot)
    print("Not fitted")
    print(nofit*100/tot)
    print(hists_DS_regDS.GetEntries()*100/tot)
    print(hists_DS_mat.GetEntries()*100/tot)
    print(hists_DS_other.GetEntries()*100/tot)
    print(hists_DS_regUS.GetEntries()*100/tot)



 
    # Clone histograms to use for plotting of untuned histograms + ratios
    # Upstream
    hists_US_data_Untuned = hists_US_data.Clone()
    hists_US_mat_Untuned = hists_US_mat.Clone()
    hists_US_other_Untuned = hists_US_other.Clone()
    hists_US_regUS_Untuned = hists_US_regUS.Clone()
    hists_US_regDS_Untuned = hists_US_regDS.Clone()
    # Downstream
    hists_DS_data_Untuned = hists_DS_data.Clone()
    hists_DS_mat_Untuned = hists_DS_mat.Clone()
    hists_DS_other_Untuned = hists_DS_other.Clone()
    hists_DS_regUS_Untuned = hists_DS_regUS.Clone()
    hists_DS_regDS_Untuned = hists_DS_regDS.Clone()

    # Clone histograms to use for plotting of tuned histograms + ratios
    # Upstream
    hists_US_data_Tuned = hists_US_data.Clone()
    hists_US_mat_Tuned = hists_US_mat.Clone()
    hists_US_other_Tuned = hists_US_other.Clone()
    hists_US_regUS_Tuned = hists_US_regUS.Clone()
    hists_US_regDS_Tuned = hists_US_regDS.Clone()
    # Downstream
    hists_DS_data_Tuned = hists_DS_data.Clone()
    hists_DS_mat_Tuned = hists_DS_mat.Clone()
    hists_DS_other_Tuned = hists_DS_other.Clone()
    hists_DS_regUS_Tuned = hists_DS_regUS.Clone()
    hists_DS_regDS_Tuned = hists_DS_regDS.Clone()

    # Create object arrrays to pass into the plotting function
    # Upstream
    mcUS = TObjArray()
    mcUS.Add(hists_US_mat_Untuned)
    mcUS.Add(hists_US_other_Untuned)
    mcUS.Add(hists_US_regUS_Untuned)
    mcUS.Add(hists_US_regDS_Untuned)
    # Downstream
    mcDS = TObjArray()
    mcDS.Add(hists_DS_mat_Untuned)
    mcDS.Add(hists_DS_other_Untuned)
    mcDS.Add(hists_DS_regUS_Untuned)
    mcDS.Add(hists_DS_regDS_Untuned)
 
    # Plot Untuned
    plotStackRatio(mnv, mcUS, hists_US_data_Untuned, mcScale, "Untuned","Upstream", var) # Upstream
    plotStackRatio(mnv, mcDS, hists_DS_data_Untuned, mcScale, "Untuned","Downstream", var) # Downstream

    # ----------------------------------------------------------------------------
    #                             TUNED HISTOGRAMS
    # ----------------------------------------------------------------------------
    # Different Chi2 than from 'get_scale_factors" script producing TunedPlasticSidebands_%s_t%s_z%s_%s_minervame6A.root files
    # Old version was not propagating the statistical error
    # Multiplication procedure, i.e. propagation of the stat and sys error consistent with my background subtraction

    # Read in scale factors
    scaleUS = scaleFactors.Get("scaleFactor_US_%s"%var)
    scaleDS = scaleFactors.Get("scaleFactor_DS_%s"%var)

    # Systematic error is constant, cv value is constant.
    # Statistical error saved in the MnvH1D is bin-by-bin statistical error which corresponds
    # to the stat error in the sideband itself. The uncertainty you propagate is the uncertainty
    # on the scale factor due to the statistical power of the sideband sample you fit, i.e should be consistent across the bin range.
    # Typical inclusive fits: stat error would correspond approximately to the sqrt(integral of the sideband) 
    # given the scale factor integrates over the sideband.

    # Statistical error: fractional 1/sqrt(N) where N number of events in the whole sideband (just given)
    # Upstream: Remove original stat error and update it with 1/sqrt(N) in each bin
    integralUS = hists_US_regUS_Tuned.Integral()
    statErrUS = 1/ m.sqrt(integralUS)
    for bin in range(scaleUS.GetNbinsX() + 1):
	    scaleUS.SetBinError(bin,0+statErrUS)

    # Downstream: Remove original stat error and update it with 1/sqrt(N) in each bin
    integralDS = hists_DS_regDS_Tuned.Integral()
    statErrDS = 1/ m.sqrt(integralDS)
    for bin in range(scaleDS.GetNbinsX() + 1):
	    scaleDS.SetBinError(bin,0+statErrDS)	

    # Multiply untuned histograms by the correspoding scale factors (MnvH1Ds)
    hists_US_regUS_Tuned.Multiply(scaleUS, hists_US_regUS_Tuned)
    hists_DS_regDS_Tuned.Multiply(scaleDS, hists_DS_regDS_Tuned)

    # Create object arrrays to pass into the plotting function
    # Upstream
    mcUSTuned = TObjArray()
    mcUSTuned.Add(hists_US_mat_Tuned)
    mcUSTuned.Add(hists_US_other_Tuned)
    mcUSTuned.Add(hists_US_regUS_Tuned)
    mcUSTuned.Add(hists_US_regDS_Tuned)
    # Downstream
    mcDSTuned = TObjArray()
    mcDSTuned.Add(hists_DS_mat_Tuned)
    mcDSTuned.Add(hists_DS_other_Tuned)
    mcDSTuned.Add(hists_DS_regUS_Tuned)
    mcDSTuned.Add(hists_DS_regDS_Tuned)
    
    # Plot Tuned
    plotStackRatio(mnv, mcUSTuned, hists_US_data_Tuned, mcScale, "Tuned","Upstream", var) # Upstream
    plotStackRatio(mnv, mcDSTuned, hists_DS_data_Tuned, mcScale, "Tuned","Downstream", var) # Downstream

    out = ROOT.TFile("Tuned_Untuned_Plastic_Sidebands.root","RECREATE")
    out.cd()
    hists_DS_regDS_Untuned.Write()
    hists_US_regUS_Untuned.Write()
    hists_DS_regDS_Tuned.Write()
    hists_US_regUS_Tuned.Write()
    out.Close()


# Left here to read in the produced root files from the 'get_scale_factors" script if needed
'''
bands = ["US", "DS"]
for band in bands:
    for var in vars:
        tunedFile = TFile.Open("TunedPlasticSidebands_%s_t%s_z%s_%s_minervame6A.root"%(band, targetID, targetZ, var))
        data = tunedFile.Get("tuned_data")
        signal = tunedFile.Get("tuned_mc_signal")
        other = tunedFile.Get("tuned_mc_other")
        us = tunedFile.Get("tuned_mc_us")
        ds = tunedFile.Get("tuned_mc_ds")

        mcTuned = TObjArray()
        mcTuned.Add(signal)
        mcTuned.Add(other)
        mcTuned.Add(us)
        mcTuned.Add(ds)

        if band == "US":
            plotStackRatio(mnv, mcTuned, data, mcScale, "TunedOrig","Upstream", var)
        if band == "DS":
            plotStackRatio(mnv, mcTuned, data, mcScale, "TunedOrig","Downstream", var)
'''

print("DONE %s %s %02s"%(plist, targetID, targetZ))
