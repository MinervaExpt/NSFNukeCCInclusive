import ROOT
from ROOT import PlotUtils
import sys
from ROOT import gStyle
from ROOT import *

tDirName = "Ratio_Unfolded_True_Histograms"
histoName = "Ratio_modelData_trueData_Stat_0_Iter_1"
histoName2 = "Ratio_modelData_trueData_Stat_0_Iter_2"
histoName3 = "Ratio_modelData_trueData_Stat_0_Iter_3"
histoName4 = "Ratio_modelData_trueData_Stat_0_Iter_4"
histoName5 = "Ratio_modelData_trueData_Stat_0_Iter_5"
histoName6 = "Ratio_modelData_trueData_Stat_0_Iter_6"
histoName7 = "Ratio_modelData_trueData_Stat_0_Iter_7"
histoName8 = "Ratio_modelData_trueData_Stat_0_Iter_8"
histoName9 = "Ratio_modelData_trueData_Stat_0_Iter_9"
histoName10 = "Ratio_modelData_trueData_Stat_0_Iter_10"
histoName20 = "Ratio_modelData_trueData_Stat_0_Iter_20"
histoName50 = "Ratio_modelData_trueData_Stat_0_Iter_50"
histoName100 = "Ratio_modelData_trueData_Stat_0_Iter_100"

lineWidth = 3

can = ROOT.TCanvas("chi2")
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

    spread.SetTitle("Stat Universe: " + univName)
    spread.GetXaxis().CenterTitle()
    spread.SetTitleOffset(1.1, "X")
    spread.GetYaxis().CenterTitle()
    spread.SetTitleOffset(1.1, "Y")
    spread.SetMinimum(0.9)
    spread.SetMaximum(1.1)
    spread.GetXaxis().SetTitle("Bjorken x")
    spread.GetYaxis().SetTitle("Unfolded Data/True Data")
    #spread.GetXaxis().SetRangeUser(1, 30)
    spread.SetLineColor(40)
    spread.Draw("HIST")
    spread.GetYaxis().SetNdivisions(505)
    #spread.GetYaxis().SetNDivisions(500)

    spread.SetLineColor(30)
    spread2.Draw("SAME HIST")
    spread3.SetLineColor(42)
    spread3.Draw("SAME HIST")
    spread4.SetLineColor(46)
    spread4.Draw("SAME HIST")
    spread5.SetLineColor(48)
    spread5.Draw("SAME HIST")
    #spread6.SetLineColor(30)
    #spread6.Draw("SAME HIST")
    #spread7.SetLineColor(38)
    #spread7.Draw("SAME HIST")
    #spread8.SetLineColor(12)
    #spread8.Draw("SAME HIST")
    #spread9.SetLineColor(6)
    #spread9.Draw("SAME HIST")
    spread10.SetLineColor(12)
    spread10.Draw("SAME HIST")
    spread20.SetLineColor(36)
    spread10.Draw("SAME HIST")
    spread50.SetLineColor(38)
    spread50.Draw("SAME HIST")
    spread100.SetLineColor(40)
    spread100.Draw("SAME HIST")

    leg = ROOT.TLegend(0.55, 0.55, 0.8, 0.9)
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
    mnv.AddHistoTitle("Fake Data == Orig. Data: Stat Universe 0", 0.04, 1)
    can.Print(fileName + "check.png")