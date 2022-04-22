import os,sys
import ROOT

infile = ROOT.TFile(sys.argv[1])
inhist = infile.Get(sys.argv[2])
threshold = float(sys.argv[3])


#n_binsx = inhist.GetNbinsX()
#n_binsy = n_binsx
n_binsx = 8
n_binsy = n_binsx

bin_boundary = []
bin_boundary.append(0)
bin_boundary.append(1)
diag = []
counter = 1
while bin_boundary[-1]<n_binsy and counter<n_binsy:
    print counter
    integral_row = inhist.Integral(0,-1,bin_boundary[-1],counter)
    diagonal_val = inhist.Integral(bin_boundary[-1],counter,bin_boundary[-1],counter)
    print counter,integral_row,diagonal_val,inhist.GetYaxis().GetBinCenter(counter)
    if(integral_row<1e-5): 
        counter+=1
        continue
    if(diagonal_val/integral_row>threshold):
        bin_boundary.append(counter)
        counter+=1
        diag.append(diagonal_val/integral_row)
        continue
    else:
        counter+=1
        continue
    
print "Optimized"
print bin_boundary
print diag
bin_edges = [0]
for el in bin_boundary[2:]:
    bin_edges.append(inhist.GetXaxis().GetBinLowEdge(el))

print bin_edges
