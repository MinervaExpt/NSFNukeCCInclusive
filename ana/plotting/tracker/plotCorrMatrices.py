import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle, gPad


ROOT.TH1.AddDirectory(False)

infile= ROOT.TFile("CrossSection_t99_z99_minervame6A_NEW.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()
mnv.SetCorrelationPalette() #kLightTenperature
#mnv.SetCorrelationPalette()

steps = ["total_corr_matrix", "stat_corr_matrix"]
vars = ["Enu", "x"]
prefixes = ['data']#, 'mc']

for step in steps:
        for prefix in prefixes:
                for var in vars:
                        print("---------------------------------------")
                        print(var)

                        mc_hist = infile.Get("%s_%s_%s"%(step,prefix, var))

                        if var == "Enu":
                                mc_hist.GetXaxis().SetTitle("E_{#bar{#nu}} bin number")
                                mc_hist.GetYaxis().SetTitle("E_{#bar{#nu}} bin number")
                                mc_hist.SetMarkerSize(1.2)

                        if var == "x":
                                mc_hist.GetXaxis().SetTitle("Bjorken x bin number")
                                mc_hist.GetYaxis().SetTitle("Bjorken x bin number")
                                mc_hist.SetMarkerSize(1.5)


                        mc_hist.GetXaxis().CenterTitle()
                        mc_hist.GetYaxis().CenterTitle()
                        mc_hist.GetZaxis().CenterTitle()


                        mc_hist.SetMaximum((1*mc_hist.GetMaximum()))
                        mc_hist.SetMinimum(-1.*(mc_hist.GetMaximum()))

                        gStyle.SetNdivisions(505,"xy") # less ticks and numbers on axis (this is VERY handy!)
                        gStyle.SetNdivisions(506,"z") 

                        mc_hist.GetZaxis().SetTitle("Bin-to-bin correlation")
                        mc_hist.Draw("COLZ")
                        print (mc_hist.GetListOfFunctions())

                        mc_hist.GetXaxis().SetTitleSize(0.05)
                        mc_hist.GetYaxis().SetTitleSize(0.05)
                        mc_hist.GetYaxis().SetTitleOffset(1.)
                        mc_hist.GetZaxis().SetTitleSize(0.045)
                        mc_hist.GetZaxis().SetTitleOffset(1.2)
                        mc_hist.GetZaxis().SetLabelSize(0.04)
                        #hist.SetMaximum(10**(-6))
                        #gStyle.SetLabelFont(42)

                        gPad.SetRightMargin(0.2)
                        '''
                        palette = mc_hist.GetListOfFunctions().FindObject("palette")
                        palette.SetX1NDC(0.81)
                        palette.SetX2NDC(0.85)
                        '''
                        #palette.SetY1NDC(0.2)
                        #palette.SetY2NDC(0.8)
                        gPad.Modified()
                        mnv.AddHistoTitle("Cross-section %s: %s (Tracker)"%(step, prefix), 0.038, 1)
                        canvas1.Print("ME6A_Tracker_%s_%s_%s.png"%(step, prefix,var))


raw_input("Done")
