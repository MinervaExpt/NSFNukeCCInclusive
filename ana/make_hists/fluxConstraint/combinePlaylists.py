import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import TParameter

ROOT.TH1.AddDirectory(False)
mcPOTAll = 0

pwd = "/pnfs/minerva/persistent/users/anezkak/default_analysis_loc/AntiNu_flux/"

infile = ROOT.TFile(str(pwd)+"Hists_EventSelection_minervame6A_FluxConstraint_optimPetal_sys_t99_z99_AntiNu.root ")
out = ROOT.TFile("CombinedPlaylists_AntiNu_sys_Enu_vs_Petal.root","RECREATE")

# declare histograms
histogram = infile.Get("selected_mc_truth_trackerC_Enu")

histogram1 = infile.Get("selected_mc_truth_trackerC_Enu")
histogram1.ClearAllErrorBands()
histogram1.AddMissingErrorBandsAndFillWithCV(histogram)
histogram1.Reset()

histogram2 = infile.Get("selected_mc_truth_waterO_Enu")
histogram2.ClearAllErrorBands()
histogram2.AddMissingErrorBandsAndFillWithCV(histogram)
histogram2.Reset()

histogram3 = infile.Get("selected_mc_truth_t1fe_Enu")
histogram3.ClearAllErrorBands()
histogram3.AddMissingErrorBandsAndFillWithCV(histogram)
histogram3.Reset()

histogram4 = infile.Get("selected_mc_truth_t1pb_Enu")
histogram4.ClearAllErrorBands()
histogram4.AddMissingErrorBandsAndFillWithCV(histogram)
histogram4.Reset()

histogram5 = infile.Get("selected_mc_truth_t2fe_Enu")
histogram5.ClearAllErrorBands()
histogram5.AddMissingErrorBandsAndFillWithCV(histogram)
histogram5.Reset()

histogram6 = infile.Get("selected_mc_truth_t2pb_Enu")
histogram6.ClearAllErrorBands()
histogram6.AddMissingErrorBandsAndFillWithCV(histogram)
histogram6.Reset()

histogram7 = infile.Get("selected_mc_truth_t3fe_Enu")
histogram7.ClearAllErrorBands()
histogram7.AddMissingErrorBandsAndFillWithCV(histogram)
histogram7.Reset()

histogram8 = infile.Get("selected_mc_truth_t3pb_Enu")
histogram8.ClearAllErrorBands()
histogram8.AddMissingErrorBandsAndFillWithCV(histogram)
histogram8.Reset()

histogram9 = infile.Get("selected_mc_truth_t3c_Enu")
histogram9.ClearAllErrorBands()
histogram9.AddMissingErrorBandsAndFillWithCV(histogram)
histogram9.Reset()

histogram10 = infile.Get("selected_mc_truth_t4pb_Enu")
histogram10.ClearAllErrorBands()
histogram10.AddMissingErrorBandsAndFillWithCV(histogram)
histogram10.Reset()

histogram11 = infile.Get("selected_mc_truth_t5fe_Enu")
histogram11.ClearAllErrorBands()
histogram11.AddMissingErrorBandsAndFillWithCV(histogram)
histogram11.Reset()

histogram12 = infile.Get("selected_mc_truth_t5pb_Enu")
histogram12.ClearAllErrorBands()
histogram12.AddMissingErrorBandsAndFillWithCV(histogram)
histogram12.Reset()

histogram13 = infile.Get("selected_mc_truth_t15fe_Enu")
histogram13.ClearAllErrorBands()
histogram13.AddMissingErrorBandsAndFillWithCV(histogram)
histogram13.Reset()

histogram14 = infile.Get("selected_mc_truth_t25fe_Enu")
histogram14.ClearAllErrorBands()
histogram14.AddMissingErrorBandsAndFillWithCV(histogram)
histogram14.Reset()

histogram15 = infile.Get("selected_mc_truth_t15pb_Enu")
histogram15.ClearAllErrorBands()
histogram15.AddMissingErrorBandsAndFillWithCV(histogram)
histogram15.Reset()

histogram16 = infile.Get("selected_mc_truth_t25pb_Enu")
histogram16.ClearAllErrorBands()
histogram16.AddMissingErrorBandsAndFillWithCV(histogram)
histogram16.Reset()

histogram17 = infile.Get("selected_mc_truth_trackerC_Enu_Petal")
histogram17.ClearAllErrorBands()
histogram17.AddMissingErrorBandsAndFillWithCV(histogram)
histogram17.Reset()

histogram18 = infile.Get("selected_mc_truth_trackerC_Petal_Enu")
histogram18.ClearAllErrorBands()
histogram18.AddMissingErrorBandsAndFillWithCV(histogram)
histogram18.Reset()



for filename in os.listdir(str(pwd)):
    if filename.endswith(".root"): 
            #print (filename)
            print(str(filename))
            infile = ROOT.TFile(str(pwd)+str(filename))
            histogram1.Add(infile.Get("selected_mc_truth_trackerC_Enu"))
            histogram2.Add(infile.Get("selected_mc_truth_waterO_Enu"))
            histogram3.Add(infile.Get("selected_mc_truth_t1fe_Enu"))
            histogram4.Add(infile.Get("selected_mc_truth_t1pb_Enu"))
            histogram5.Add(infile.Get("selected_mc_truth_t2fe_Enu"))
            histogram6.Add(infile.Get("selected_mc_truth_t2pb_Enu"))
            histogram7.Add(infile.Get("selected_mc_truth_t3fe_Enu"))
            histogram8.Add(infile.Get("selected_mc_truth_t3pb_Enu"))
            histogram9.Add(infile.Get("selected_mc_truth_t3c_Enu"))
            histogram10.Add(infile.Get("selected_mc_truth_t4pb_Enu"))
            histogram11.Add(infile.Get("selected_mc_truth_t5fe_Enu"))
            histogram12.Add(infile.Get("selected_mc_truth_t5pb_Enu"))
            histogram13.Add(infile.Get("selected_mc_truth_t15fe_Enu"))
            histogram14.Add(infile.Get("selected_mc_truth_t25fe_Enu"))
            histogram15.Add(infile.Get("selected_mc_truth_t15pb_Enu"))
            histogram16.Add(infile.Get("selected_mc_truth_t25pb_Enu"))
            histogram17.Add(infile.Get("selected_mc_truth_trackerC_Enu_Petal"))
            histogram18.Add(infile.Get("selected_mc_truth_trackerC_Petal_Enu"))
            mcPOTAll += infile.Get("MCPOT").GetVal()
            print(infile.Get("MCPOT").GetVal())
            print(mcPOTAll)
            print("Added Histograms")
        #continue

out.cd()
histogram1.Write()
histogram2.Write()
histogram3.Write()
histogram4.Write()
histogram5.Write()
histogram6.Write()
histogram7.Write()
histogram8.Write()
histogram9.Write()
histogram10.Write()
histogram11.Write()
histogram12.Write()
histogram13.Write()
histogram14.Write()
histogram15.Write()
histogram16.Write()
histogram17.Write()
histogram18.Write()
mcPOTout = TParameter(float)("MCPOT", mcPOTAll)
mcPOTout.Write()

raw_input("Done")

