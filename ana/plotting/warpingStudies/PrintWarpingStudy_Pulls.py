import ROOT
from ROOT import PlotUtils
import sys
from ROOT import gStyle
from ROOT import *
import math as m

tDirName = "Pull_Histograms"
histoName = "Stat_0_Iter_1_pull"
histoName2 = "Stat_0_Iter_2_pull"
histoName3 = "Stat_0_Iter_3_pull"
histoName4 = "Stat_0_Iter_4_pull"
histoName5 = "Stat_0_Iter_5_pull"
histoName6 = "Stat_0_Iter_6_pull"
histoName7 = "Stat_0_Iter_7_pull"
histoName8 = "Stat_0_Iter_8_pull"
histoName9 = "Stat_0_Iter_9_pull"
histoName10 = "Stat_0_Iter_10_pull"
histoName20 = "Stat_0_Iter_20_pull"
histoName50 = "Stat_0_Iter_50_pull"
histoName100 = "Stat_0_Iter_100_pull"

lineWidth = 3

can = ROOT.TCanvas("pulls")
mnv = PlotUtils.MnvPlotter()
#mnv.SetRedHeatPalette()
mnv.SetROOT6Palette(57) #kBird

mcColors = MnvColors.GetColors(MnvColors.kOkabeItoDarkPalette)
nextColor = 0

for fileName in sys.argv[1:]:
    myFile = ROOT.TFile.Open(fileName)

    #Try to infer a useful universe name from the file name
    univName = "TODO"
    #univName = fileName[fileName.find("merged") + len("merged") + 1:fileName.find(".root")]
    #if fileName.find("SuSA") != -1: #The SuSA warp is a stand-alone CV, so it needs special treatment
    #  univName = "SuSA"

    spread = myFile.Get(tDirName).Get(histoName)
    spread2 = myFile.Get(tDirName).Get(histoName2)
    spread3 = myFile.Get(tDirName).Get(histoName3)
    spread4 = myFile.Get(tDirName).Get(histoName4)
    spread5 = myFile.Get(tDirName).Get(histoName5)
    spread6 = myFile.Get(tDirName).Get(histoName6)
    spread7 = myFile.Get(tDirName).Get(histoName7)
    spread8 = myFile.Get(tDirName).Get(histoName8)
    spread9 = myFile.Get(tDirName).Get(histoName9)
    spread10 = myFile.Get(tDirName).Get(histoName10)
    spread20 = myFile.Get(tDirName).Get(histoName20)
    spread50 = myFile.Get(tDirName).Get(histoName50)
    spread100 = myFile.Get(tDirName).Get(histoName100)

    #Infer number of degrees of freedom from y axis title
    axisTitle = spread.GetYaxis().GetTitle()
    #yNDF = int(axisTitle[axisTitle.find("ndf=") + 4:axisTitle.find(")")])

    # SET STATISTICAL BIN ERROR
    histos = [spread, spread2, spread3, spread4, spread5, spread6, spread7, spread8, spread9, spread10, spread20, spread50, spread100]
    for histo in histos:
        for bin in range(1, histo.GetNbinsX()+1):
            content = histo.GetBinContent(bin)
            err = m.sqrt(abs(content))
            histo.SetBinError(bin, err)
            print(histo.GetBinError(bin))


    

    #spread.SetTitle("Pull Distribution " + univName)
    spread.GetXaxis().CenterTitle()
    spread.SetTitleOffset(1.1, "X")
    spread.GetYaxis().CenterTitle()
    spread.SetTitleOffset(1.1, "Y")
    spread.SetMinimum(-2)
    spread.SetMaximum(2)
    if 'x.root' in str(myFile):
        spread.GetXaxis().SetTitle("Bjorken x")
    else: 
        spread.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
    spread.GetYaxis().SetTitle("(reco-true)/unc")
    #spread.GetXaxis().SetRangeUser(1, 30)
    spread.SetLineColor(40)
    spread.SetMarkerColor(40)
    spread.Draw(" E")
    spread.GetYaxis().SetNdivisions(505)
    #spread.GetYaxis().SetNDivisions(500)

    spread.SetLineColor(30)
    spread.SetMarkerColor(30)
    spread2.Draw("SAME  E")
    spread3.SetLineColor(42)
    spread3.SetMarkerColor(42)
    spread3.Draw("SAME  E")
    spread4.SetLineColor(46)
    spread4.SetMarkerColor(46)
    spread4.Draw("SAME  E")
    spread5.SetLineColor(48)
    spread5.SetMarkerColor(48)
    spread5.Draw("SAME  E")
    #spread6.SetLineColor(30)
    #spread6.Draw("SAME HIST")
    #spread7.SetLineColor(38)
    #spread7.Draw("SAME HIST")
    #spread8.SetLineColor(12)
    #spread8.Draw("SAME HIST")
    #spread9.SetLineColor(6)
    #spread9.Draw("SAME HIST")
    spread10.SetLineColor(12)
    spread10.SetMarkerColor(12)
    spread10.Draw("SAME  E")
    spread20.SetLineColor(36)
    spread20.SetMarkerColor(36)
    spread10.Draw("SAME  E")
    spread50.SetLineColor(38)
    spread50.SetMarkerColor(38)
    spread50.Draw("SAME  E")
    spread100.SetLineColor(40)
    spread100.SetMarkerColor(40)
    spread100.Draw("SAME  E")

    leg = ROOT.TLegend(0.35, 0.70, 0.8, 0.9)
    leg.SetNColumns(2)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.03)
    leg.SetTextFont(42)
    leg.AddEntry(spread, "Iteration 1", "l")
    leg.AddEntry(spread2, "Iteration 2", "l")
    leg.AddEntry(spread3, "Iteration 3", "l")
    leg.AddEntry(spread4, "Iteration 4", "l")
    leg.AddEntry(spread5, "Iteration 5", "l")
    #leg.AddEntry(spread6, "Iteration 6", "l")
    #leg.AddEntry(spread7, "Iteration 7", "l")
    #leg.AddEntry(spread8, "Iteration 8", "l")
    #leg.AddEntry(spread9, "Iteration 9", "l")
    leg.AddEntry(spread10, "Iteration 10", "l")
    leg.AddEntry(spread20, "Iteration 20", "l")
    leg.AddEntry(spread50, "Iteration 50", "l")
    leg.AddEntry(spread100, "Iteration 100", "l")
    #leg.AddEntry(doubleNDFLine, "2x Number of Bins", "l")
    #leg.AddEntry(iterLine, str(iterChosen) + " iterations", "l")
    leg.Draw()

    can.Modified()
    can.Update()
    mnv.AddHistoTitle("Pull Distribution: CV ", 0.04, 1)
    can.Print(fileName + "_pulls.png")