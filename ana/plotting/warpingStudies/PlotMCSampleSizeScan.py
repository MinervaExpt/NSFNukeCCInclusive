import os,sys
import ROOT
import PlotUtils
from ROOT import TLine

ROOT.TH1.AddDirectory(False)
mnv = PlotUtils.MnvPlotter()
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))


def sortInputFiles(filelist):
    newlist = []
    factors = []
    for f in filelist:
        fn = f.rstrip("\n")
        mcfactor = float(fn.split(".txt")[-2].split("_")[-1])
        factors.append(mcfactor)

    factors = sorted(factors)
    for f in factors:
        for fi in filelist:
            mcfactor = float(fi.split(".txt")[-2].split("_")[-1])
            if(f==mcfactor and fi not in newlist): newlist.append(fi)

    return newlist

colors = [1,20,30,40,50,60,70,80,90,100,51,1,20,30,40,50,60,70,80,90,100,51]#22 allowed inputs
styles = [1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2]#22 allowed inputs

inputfiles = sys.argv[1:]
print("Input files")
print(inputfiles)
newlist = sortInputFiles(inputfiles)
print("Sorted input files")
print(newlist)
for el in newlist:
    print el.rstrip("\n")
legend = ROOT.TLegend(0.2,0.7,0.82,0.9)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetNColumns(3)
legend.SetTextFont(42)

yNFD = 0
name = ""


chi2_by_iter_by_file = []
for f in newlist:
    fn = f.rstrip("\n")
    temp_file = open(fn,"r").readlines()
    n_uni = 0
    current_iter = 0
    chi2_by_iter = []
    for l in temp_file:
        temp_line = l.split()
        chi2 = float(temp_line[0])
        iteration = int(temp_line[1])
        universe = int(temp_line[2])
        if(universe>n_uni):
            n_uni=universe+1
        if(iteration!=current_iter):
            if(current_iter!=0):
                chi2_by_iter[-1][1]/=n_uni
            chi2_by_iter.append([iteration,chi2])
            current_iter=iteration
        else:
            chi2_by_iter[-1][1]+=chi2
        if l == temp_file[-1]:
            chi2_by_iter[-1][1]/=n_uni
#            chi2_by_iter.append([iteration,chi2])

    print chi2_by_iter
    chi2_by_iter_by_file.append(chi2_by_iter)

mygraphs = []
for i in range(0,len(chi2_by_iter_by_file)):
    tmpgraph = ROOT.TGraph()
    for j,el in enumerate(chi2_by_iter_by_file[i]):
        tmpgraph.SetPoint(j,el[0],el[1])
        
    tmpgraph.SetLineColor(colors[i])
    tmpgraph.SetLineStyle(styles[i])
    tmpgraph.SetMarkerColor(colors[i])
    mygraphs.append(tmpgraph)
        


for i,g in enumerate(mygraphs):
    legend.AddEntry(g,"%s"%(float(newlist[i].split(".txt")[-2].split("_")[-2])),"l")
    legend.Draw("SAME")

    if(i==0):
        g.Draw("ALP")
        g.SetMinimum(0)
    else:
        g.Draw("SAMELP")

    # variable
    print(newlist[i].split(".txt")[-2].split("_")[-3].split(".root")[-2])
    var = newlist[i].split(".txt")[-2].split("_")[-3].split(".root")[-2]
    # f factor
    print(float(newlist[i].split(".txt")[-2].split("_")[-1]))
    f = float(newlist[i].split(".txt")[-2].split("_")[-1])
    if var == "Enu":
        mnv.AddHistoTitle("Antineutrino Energy (GeV), f = %f"%(f),  0.04, 1)
        yNDF = 10
        g.SetMaximum(13)
        g.SetMinimum(7)
        #g.SetMaximum(20)
        #g.SetMinimum(0)
        g.GetXaxis().SetRangeUser(0, 10)
        g.GetYaxis().SetNdivisions(505)
        ndfLine = ROOT.TLine(0, yNDF, 11, yNDF)
        ndfLine.SetLineWidth(3)
        ndfLine.SetLineStyle(ROOT.kDashed)
        ndfLine.Draw("SAME") 
    if var == "x":
        mnv.AddHistoTitle("Bjorken x, f = %f"%(f),  0.04, 1)
        yNDF = 6
        g.SetMaximum(9)
        g.SetMinimum(3)
        g.GetXaxis().SetRangeUser(0, 10)
        ndfLine = ROOT.TLine(0, yNDF, 11, yNDF)
        ndfLine.SetLineWidth(3)
        ndfLine.SetLineStyle(ROOT.kDashed)
        ndfLine.Draw("SAME") 
    if var == "pZmu1D":
        mnv.AddHistoTitle("Muon p_{Z} (GeV/c), f = %f"%(f),  0.04, 1)
        yNDF = 13
        g.SetMaximum(16)
        g.SetMinimum(10)
        g.GetXaxis().SetRangeUser(0, 10)
        ndfLine = ROOT.TLine(0, yNDF, 11, yNDF)
        ndfLine.SetLineWidth(3)
        ndfLine.SetLineStyle(ROOT.kDashed)
        ndfLine.Draw("SAME") 
    if var == "pTmu1D":
        mnv.AddHistoTitle("Muon p_{T} (GeV/c), f = %f"%(f),  0.04, 1)
        yNDF = 13
        g.SetMaximum(16)
        g.SetMinimum(10)
        g.GetXaxis().SetRangeUser(0, 10)
        ndfLine = ROOT.TLine(0, yNDF, 11, yNDF)
        ndfLine.SetLineWidth(3)
        ndfLine.SetLineStyle(ROOT.kDashed)
        ndfLine.Draw("SAME") 
    if var == "ThetamuDeg":
        mnv.AddHistoTitle("Muon Angle (Deg), f = %f"%(f),  0.04, 1)
        yNDF = 9
        g.SetMaximum(12)
        g.SetMinimum(6)
        g.GetXaxis().SetRangeUser(0, 10)
        ndfLine = ROOT.TLine(0, yNDF, 11, yNDF)
        ndfLine.SetLineWidth(3)
        ndfLine.SetLineStyle(ROOT.kDashed)
        ndfLine.Draw("SAME") 
 
    g.GetYaxis().SetTitle("Average #chi^{2} (NDF =%d)"%yNDF) # number of bins
    g.GetYaxis().CenterTitle()

    g.GetXaxis().SetTitle("(Unfolded Data:True Data) # of Iterations") # number of bins
    g.GetXaxis().CenterTitle()

    name = newlist[i].split(".root")[-2]

canvas1.Modified()
canvas1.Print("%s_fStudy.png"%name)
