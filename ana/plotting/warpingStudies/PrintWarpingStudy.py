import ROOT
from ROOT import PlotUtils
import sys

chi2SummaryDir = "Chi2_Iteration_Dists"
chi2SummaryName = "h_chi2_modelData_trueData_iter_chi2"
medianHistName = "h_median_chi2_modelData_trueData_iter_chi2"
meanChi2ProfileName = "m_avg_chi2_modelData_trueData_iter_chi2_truncated"

lineWidth = 3
iterChosen = 3

can = ROOT.TCanvas("chi2")
mnv = PlotUtils.MnvPlotter()
#mnv.SetRedHeatPalette()
mnv.SetROOT6Palette(57) #kBird

fileName = sys.argv[1]
warp = sys.argv[2] 
var = sys.argv[3]

myFile = ROOT.TFile.Open(fileName)

#Try to infer a useful universe name from the file name
univName = "TODO"
#univName = fileName[fileName.find("merged") + len("merged") + 1:fileName.find(".root")]
#if fileName.find("SuSA") != -1: #The SuSA warp is a stand-alone CV, so it needs special treatment
#  univName = "SuSA"

spread = myFile.Get(chi2SummaryDir).Get(chi2SummaryName)

#Infer number of degrees of freedom from y axis title
axisTitle = spread.GetYaxis().GetTitle()
yNDF = int(axisTitle[axisTitle.find("ndf=") + 4:axisTitle.find(")")])

spread.SetTitle("Universe: " + univName)
spread.GetXaxis().CenterTitle()
spread.SetTitleOffset(1.1, "X")
spread.GetYaxis().CenterTitle()
spread.SetTitleOffset(1.1, "Y")
spread.GetYaxis().SetRangeUser(1, 60)
#spread.GetXaxis().SetRangeUser(1, 30)
spread.Draw("colz")

profile = myFile.Get(chi2SummaryDir).Get(meanChi2ProfileName)
profile.SetTitle("Mean Chi2")
profile.SetLineWidth(lineWidth)
profile.SetLineColor(ROOT.kRed)
profile.SetMarkerStyle(0)
profile.Draw("SAME")

median = myFile.Get(chi2SummaryDir).Get(medianHistName)
median.SetTitle("Median Chi2")
median.SetLineWidth(lineWidth)
median.SetLineColor(ROOT.kBlack)
median.Draw("HIST SAME")

#Draw lines at number of degrees of freedom and 2x NDF
ndfLine = ROOT.TLine(1, yNDF, spread.GetXaxis().GetXmax(), yNDF)
ndfLine.SetLineWidth(lineWidth)
ndfLine.SetLineStyle(ROOT.kDashed)
ndfLine.Draw()

#doubleNDFLine = ROOT.TLine(1, 2*yNDF, spread.GetXaxis().GetXmax(), 2*yNDF)
#doubleNDFLine.SetLineColor(ROOT.kRed)
#doubleNDFLine.SetLineWidth(lineWidth)
#doubleNDFLine.SetLineStyle(ROOT.kDashed)
#doubleNDFLine.Draw()

#Draw a line at the chosen number of iterations.
#iterLine = ROOT.TLine(iterChosen + 0.5, 0, iterChosen + 0.5, spread.GetYaxis().GetXmax())
#iterLine.SetLineWidth(lineWidth)
#iterLine.SetLineStyle(ROOT.kDotted)
#iterLine.Draw()

#Make a custom legend because I don't want to include the 2D histogram
#while I must Draw() it first to set the right axis limits.
leg = ROOT.TLegend(0.55, 0.7, 0.8, 0.9)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.035)
leg.SetTextFont(42)
leg.AddEntry(profile)
leg.AddEntry(median)
leg.AddEntry(ndfLine, "Number of bins", "l")
#leg.AddEntry(doubleNDFLine, "2x Number of Bins", "l")
#leg.AddEntry(iterLine, str(iterChosen) + " iterations", "l")
leg.Draw()
can.SetLogx()


if var == "Enu":
    mnv.AddHistoTitle("Fake Data == %s Warp Antineutrino Energy"%(warp), 0.032, 1)
if var == "x":
    mnv.AddHistoTitle("Fake Data == %s Warp Bjorken x" %(warp), 0.032, 1)
if var == "pTmu1D":
    mnv.AddHistoTitle("Fake Data == %s Warp Muon Transverse Momentum" %(warp), 0.032, 1)
if var == "pZmu1D":
    mnv.AddHistoTitle("Fake Data == %s Warp Muon Longitudinal Momentum" %(warp), 0.032, 1)
if var == "ThetamuDeg":
    mnv.AddHistoTitle("Fake Data == %s Warp Muon Angle" %(warp), 0.032, 1)

can.Print(fileName + ".png")