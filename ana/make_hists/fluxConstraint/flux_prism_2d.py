import os
import sys
import math
import ROOT
ROOT.gROOT.SetBatch(True)
#import makexsec_new
import PlotUtils
import numpy as np
from matplotlib import pyplot as plt

# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)


c = ROOT.TCanvas("c", "c", 800, 800)

plotter = PlotUtils.MnvPlotter()
plotter.SetRedHeatPalette()


LIMIT_ENU = False

def copy(hist, isROOT = False):
    #print "copy: hist n vert errors bars", hist.GetNVertErrorBands()
    #print "copy hist n lat errors bars", hist.GetNLatErrorBands()
    #new = type(hist)(hist)
    #new.TransferErrorBands(hist, False)
    new = make_clean_copy(hist, not isROOT)
    new.Add(hist)
    #print "copy: new n vert errors bars", new.GetNVertErrorBands()
    #print "copy: new n lat errors bars", new.GetNLatErrorBands()
    return new
    
def make_clean_copy(hist, add_error_bands = True):
    # Try to copy the histogram and make a blank one
    new = type(hist)(hist)
    # Removes seg fault by taking direct control of this histogram.
    # See https://root.cern.ch/root/htmldoc/guides/users-guide/ObjectOwnership.html
    new.SetDirectory(0)
    #new = copy(hist)
    new.Reset()
    if add_error_bands:
        new.AddMissingErrorBandsAndFillWithCV(hist)
    return new
    

def main(hist_2d, target_hist, reg_lambda = 0, material = 'lead_t15'):
    number_of_parameters = hist_2d.GetNbinsY()
    
    nbins = target_hist.GetNbinsX()
    assert target_hist.GetNbinsX() == hist_2d.GetNbinsX(), "The number of bins is different"
    
    stats_hist = hist_2d.ProjectionY()
    total_stats = float(stats_hist.Integral())
    stats_vector = []
    for b in range(1, number_of_parameters + 1):
        relative_stats = stats_hist.GetBinContent(b) / total_stats
        stats_vector.append(relative_stats)
    assert sum(stats_vector) > 0 and total_stats > 0, "Something's not right with the regularization stats"
        
    x_bins = range(1, hist_2d.GetNbinsX() + 1)
    y_bins = range(1, hist_2d.GetNbinsY() + 1)
    if LIMIT_ENU:
        bin_limit_start = hist_2d.GetXaxis().FindBin(4)
        bin_limit_end = hist_2d.GetXaxis().FindBin(20)
        x_bins = x_bins[bin_limit_start:bin_limit_end + 1]
    
    def utility_function(npar, gin, f, par, iflag, verbose = False):
        chi2 = 0
        for x, bx in enumerate(x_bins):
            target = target_hist.GetBinContent(bx)
            total_content = 0
            total_error = target_hist.GetBinError(bx)**2
            for y, by in enumerate(y_bins):
                global_bin = hist_2d.GetBin(bx, by)
                current_content = hist_2d.GetBinContent(global_bin)
                current_error = hist_2d.GetBinError(global_bin)
                p = par[y]
                total_content += current_content * p
                total_error += current_error**2 * p
            chi2_bin = 0
            if total_error > 0:
                chi2_bin = (target - total_content)**2 / total_error
            else:
                assert abs(total_content) + abs(target) == 0, "Found zero error bin with non-zero content. total_content: {} target: {} error {}".format(total_content, target, total_error)
            chi2 += chi2_bin
        
        
        regularization = 0
        for b, relative_stats in enumerate(stats_vector):
            # Stats error scales as dist from 1 squared
            p = par[b]
            foo = relative_stats * (p - 1)**2
            regularization += foo * 1000
        
        
        output = chi2 + reg_lambda * regularization
        f[0] = output
        return chi2, regularization


    m = ROOT.TMinuit(number_of_parameters)
    m.SetPrintLevel()
    m.SetFCN(utility_function)
    import array
    ierflg = array.array("i", [0])
    CALL_LIMIT = 500000 #original 500000
    arglist = array.array("d", [CALL_LIMIT, 1])
    
    for parameter in range(number_of_parameters):
        parameter_name = "p_%s" % parameter
        # Start with an even weighting
        default = 0 # originally 1
        default = 2 if parameter < 5 or parameter > 11 else 0 # original
        default_min = 0
        if material in {'carbon', 'lead', 'iron_t15'}: # this is empirical, inherent to the data set
            default_min = -1 # originally zero
        default_max = 5.0
        default_step = 0.001
        m.mnparm(parameter, parameter_name, default, default_step, default_min, default_max, ierflg)
        
    
    fit_result = m.mnexcm("MIGRAD", arglist, number_of_parameters, ierflg)
    print "Error flags", ierflg
    print "Fit result:", fit_result
    if ierflg[0] != 0:
        print "Fit failed!"
        return None
    parameters = list()
    parameter_errors = list()
    for i in range(number_of_parameters):
        value = array.array("d", [0.0])
        error = array.array("d", [0.0])
        m.GetParameter(i, value, error)
        parameters.append(value[0])
        parameter_errors.append(error[0])
    
    chi2, reg = utility_function(len(parameters), None, [0], parameters, None)
    print("From the fitter:")
    print reg_lambda, chi2, reg, parameters, parameter_errors
    
    # --------------------------------------------------------------------------------
    # NEW: make param zero if negative Anezka February 2022
    if material in {'carbon', 'lead', 'iron_t15'}:
        new_parameters = []
        new_parameter_errors = [ ]
        for par,err in zip(parameters,parameter_errors):
            if par < 0.0:
                err = math.sqrt(par**2 + err**2)
                par = 0.0
                new_parameter_errors.append(err)
                new_parameters.append(par)
            else:
                new_parameter_errors.append(err)
                new_parameters.append(par)
            
        print reg_lambda, chi2, reg, new_parameters, new_parameter_errors
        return new_parameters, new_parameter_errors, chi2, reg

        #new_parameters = map(lambda x: max(x,0),parameters)
        #print reg_lambda, chi2, reg, new_parameters, new_parameter_errors
    else:
        print("Final:")
        print reg_lambda, chi2, reg, parameters, parameter_errors
        return parameters, parameter_errors, chi2, reg
    
def GetResultingHist(param_hist, hist_2d):
    temp = copy(hist_2d)
    x_bins = range(1, temp.GetNbinsX() + 1)
    y_bins = range(1, temp.GetNbinsY() + 1)
    for bx in x_bins:
        for by in y_bins:
            target_value = param_hist.GetBinContent(by)
            target_error = param_hist.GetBinError(by)
            b = temp.GetBin(bx, by)
            content = temp.GetBinContent(b)
            error = temp.GetBinError(b)
            new_content = target_value * content
            new_error = target_value * error
            temp.SetBinContent(b, new_content)
            temp.SetBinError(b, new_error)
    target_hist = temp.ProjectionX()
    return target_hist
    
def GetResultingHistFull(param_hist, hist_2d): 
    #hist_2d_old = hist_2d
    #hist_2d = makexsec_new.make_clean_copy(hist_2d)
    #hist_2d.Add(hist_2d_old)
    #hist_2d.AddMissingErrorBandsAndFillWithCV(param_hist)
    out = GetResultingHist(param_hist, hist_2d)
    out.AddMissingErrorBandsAndFillWithCV(param_hist)
    for name in param_hist.GetVertErrorBandNames():
        vert = param_hist.GetVertErrorBand(name)
        vert_out = out.GetVertErrorBand(name)
        n = vert.GetNHists()
        for i in range(n):
            vert_param_hist = vert.GetHist(i)
            assert vert_param_hist != None
            hist = GetResultingHist(vert_param_hist, hist_2d)
            temp = vert_out.GetHist(i)
            nbins = temp.GetNbinsX()
            for j in range(0, nbins + 1):
                temp.SetBinContent(j, hist.GetBinContent(j))
                temp.SetBinError(j, hist.GetBinError(j))
            #vert_out.SetHist(i, hist)
    return out
    
    
def make_param_hist(parameters, parameter_errors, name="param_hist", param_hist = None):
    if param_hist == None:
        param_hist = PlotUtils.MnvH1D("param_hist", "param_hist;Daisy Bin;Weight", len(parameters), 0, len(parameters))
    for i in range(len(parameters)):
        b = i + 1
        p = parameters[i]
        perr = parameter_errors[i]
        param_hist.SetBinContent(b, p)
        param_hist.SetBinError(b, perr)
    return param_hist

'''    
def unit_test():
    tf = ROOT.TFile("../plots/2021-05-17/daisy_vs_flux/daisy_vs_flux/separated/tracker.root")
    assert not tf.IsZombie(), "Can't read file"
    tf.ls()
    hist_2d = tf.Get("truecarbon")
    assert hist_2d != None, "Can't get hist_2d"
    
    target_param_hist = ROOT.TH1D("target_param_hist", "target_param_hist", hist_2d.GetNbinsY(), 0, hist_2d.GetNbinsY())
    target_param_hist.SetBinContent(1, 1)
    target_param_hist.SetBinContent(2, 1)
    target_param_hist.SetBinContent(3, 1)
    target_param_hist.SetBinContent(4, 1)
    target_param_hist.SetBinContent(5, 1)
    target_param_hist.SetBinContent(6, 1.5)
    target_param_hist.SetBinContent(7, 1)
    
    target_hist = GetResultingHist(target_param_hist, hist_2d)
    
    parameters, parameter_errors, chi2, reg = main(hist_2d, target_hist, 10000)
    param_hist = make_param_hist(parameters, parameter_errors)
        
        
    hist = GetResultingHist(param_hist, hist_2d)
    target_hist.SetLineColor(ROOT.kRed)
    target_hist.Draw("hist")
    foo = max(target_hist.GetMaximum(), hist.GetMaximum()) * 1.2
    target_hist.GetYaxis().SetRangeUser(0, foo)
    target_hist.GetXaxis().SetRangeUser(0, 20)
    hist.Draw("same e")
    c.Print("flux_prism_2d_output/unit_test_fluxes.png")
    
        
    
    target_param_hist.SetLineColor(ROOT.kRed)
    target_param_hist.Draw("hist")
    target_param_hist.GetYaxis().SetRangeUser(0, 2)
    param_hist.Draw("same e")
    c.Print("flux_prism_2d_output/unit_test_parameters.png")
'''
'''    
def unit_test_with_target(reg_lambda = 0):
    
    name = "{:06d}".format(reg_lambda)
    
    input_directory = "../plots/2021-05-18/daisy_vs_flux/daisy_vs_flux_hd_no_circle"
    
    tf = ROOT.TFile("%s/separated/tracker.root" % input_directory)
    assert not tf.IsZombie(), "Can't read file"
    tf.ls()
    hist_2d = tf.Get("truecarbon")
    assert hist_2d != None, "Can't get hist_2d"
    target_hist = hist_2d.ProjectionX()
    
    tf2 = ROOT.TFile("%s/flux_carbon.root" % input_directory)
    assert not tf2.IsZombie(), "Can't read second file"
    tf2.ls()
    flux1 = tf2.Get("flux")
    assert flux1 != None, "Can't read flux1"
    
    tf3 = ROOT.TFile("%s/flux_tracker.root" % input_directory)
    assert not tf3.IsZombie(), "Can't read third file"
    tf3.ls()
    flux2 = tf3.Get("flux")
    assert flux2 != None, "Can't read flux2"
    
    ratio = make_clean_copy(flux1)
    ratio.Divide(flux1, flux2)
    
    target_hist.Multiply(target_hist, ratio) # Target is in numerator, tracker in denom, so multiply
    
    foo = main(hist_2d, target_hist, reg_lambda / 1000.0)
    if foo == None:
        print "Couldn't fit lambda =", reg_lambda
        return
    parameters, parameter_errors, chi2, reg = foo
    param_hist = make_param_hist(parameters, parameter_errors)
    
    hist = GetResultingHist(param_hist, hist_2d)
    target_hist.SetLineColor(ROOT.kRed)
    target_hist.Draw("hist")
    foo = max(target_hist.GetMaximum(), hist.GetMaximum()) * 1.2
    target_hist.GetYaxis().SetRangeUser(0, foo)
    target_hist.GetXaxis().SetRangeUser(0, 20)
    hist.Draw("same e")
    c.Print("flux_prism_2d_output/unit_test_with_target_fluxes_{}.png".format(name))
    
    flux_ratio = make_clean_copy(flux1)
    flux_ratio.Divide(hist, target_hist)
    flux_ratio.GetYaxis().SetRangeUser(0.5, 1.5)
    flux_ratio.Draw("e")
    c.Print("flux_prism_2d_output/unit_test_with_target_flux_ratio_{}.png".format(name))
    flux_ratio.GetXaxis().SetRangeUser(0, 100)
    flux_ratio.Draw("e")
    c.Print("flux_prism_2d_output/unit_test_with_target_flux_ratio_zoom_out_{}.png".format(name))
    
    flux_ratio_no_tune = make_clean_copy(flux1)
    no_tune = hist_2d.ProjectionX()
    flux_ratio_no_tune.Divide(no_tune, target_hist)
    flux_ratio_no_tune.Draw("e")
    c.Print("flux_prism_2d_output/unit_test_with_target_flux_ratio_no_tune_{}.png".format(name))
    
        
    
    param_hist.Draw("e")
    param_hist.GetYaxis().SetRangeUser(0, 2)
    c.Print("flux_prism_2d_output/unit_test_with_target_parameters_{}.png".format(name))
'''    
'''
def unit_test_l_curve():
    
    
    input_directory = "../plots/2021-05-17/daisy_vs_flux/daisy_vs_flux"
    input_directory += "_hd"
    
    tf = ROOT.TFile("%s/separated/tracker.root" % input_directory)
    assert not tf.IsZombie(), "Can't read file"
    tf.ls()
    hist_2d = tf.Get("truecarbon")
    assert hist_2d != None, "Can't get hist_2d"
    target_hist = hist_2d.ProjectionX()
    
    tf2 = ROOT.TFile("%s/flux_carbon.root" % input_directory)
    assert not tf2.IsZombie(), "Can't read second file"
    tf2.ls()
    flux1 = tf2.Get("flux")
    assert flux1 != None, "Can't read flux1"
    
    tf3 = ROOT.TFile("%s/flux_tracker.root" % input_directory)
    assert not tf3.IsZombie(), "Can't read third file"
    tf3.ls()
    flux2 = tf3.Get("flux")
    assert flux2 != None, "Can't read flux2"
    
    ratio = make_clean_copy(flux1)
    ratio.Divide(flux1, flux2)
    
    target_hist.Multiply(target_hist, ratio) # Target is in numerator, tracker in denom, so multiply
    
    output = dict()
    for reg_lambda in [0, 1, 10, 50, 75, 100, 150, 200, 250, 300, 400, 500, 750, 1000, 10000]:
        foo = main(hist_2d, target_hist, reg_lambda / 1000.0)
        if foo == None: 
            print "Failed on lambda =", reg_lambda
            continue
        parameters, parameter_errors, chi2, reg = foo
        output[reg_lambda] = chi2, reg
    for r in sorted(output):
        print r, output[r][0], output[r][1]
'''
'''        
def unit_test_poisson():
    
    
    input_directory = "../plots/2021-05-18/daisy_vs_flux/daisy_vs_flux_hd_no_circle"
    
    tf = ROOT.TFile("%s/separated/tracker.root" % input_directory)
    assert not tf.IsZombie(), "Can't read file"
    tf.ls()
    hist_2d = tf.Get("truecarbon")
    assert hist_2d != None, "Can't get hist_2d"
    target_hist = hist_2d.ProjectionX()
    
    tf2 = ROOT.TFile("%s/flux_carbon.root" % input_directory)
    assert not tf2.IsZombie(), "Can't read second file"
    tf2.ls()
    flux1 = tf2.Get("flux")
    assert flux1 != None, "Can't read flux1"
    
    tf3 = ROOT.TFile("%s/flux_tracker.root" % input_directory)
    assert not tf3.IsZombie(), "Can't read third file"
    tf3.ls()
    flux2 = tf3.Get("flux")
    assert flux2 != None, "Can't read flux2"
    
    
    rebin_hists = True
    if rebin_hists:
        hist_2d = rebin(hist_2d, True)
        target_hist = rebin(target_hist)
        flux1 = rebin(flux1)
        flux2 = rebin(flux2)
    
    ratio = make_clean_copy(flux1)
    ratio.Divide(flux1, flux2)
    
    target_hist.Multiply(target_hist, ratio) # Target is in numerator, tracker in denom, so multiply
    
    myrand = ROOT.TRandom3(1234)
    poisson = lambda x: myrand.PoissonD(x)
    N_THROWS = 10000
    THROW_TARGET = True
    THROW_HIST_2D = True
    output = dict()
    n_parameters = hist_2d.GetNbinsY()
    param_hist_2d = ROOT.TH2D("param_hist", "param_hist;Flux Bin;Weight", n_parameters, 0, n_parameters, 500, 0, 2.5)
    i = 0
    done = False
    import time
    start_time = time.clock()
    max_time = 0.5 # Minutes
    end_time = start_time + 60 * max_time
    fails = 0
    tries = 0
    successes = 0
    while not done:
        if i > N_THROWS:
            done = True
        current_time = time.clock()
        if current_time > end_time:
            done = True
        if not done:
            print "\n\n\nOn throw number %s" % i
            new_target_hist = copy(target_hist)
            new_hist_2d = copy(hist_2d)
            if THROW_TARGET:
                for b in range(1, target_hist.GetNbinsX() + 1):
                    old_content = target_hist.GetBinContent(b)
                    if i == 0:
                        new_content = old_content
                    else:
                        new_content = poisson(old_content)
                    # This would be right if I had not scaled the target hist by the ratio
                    #new_error = math.sqrt(new_content)
                    old_error = target_hist.GetBinError(b)
                    scale_factor = 1.0
                    if old_content > 0:
                        scale_factor = new_content / float(old_content)
                    new_error = scale_factor * old_error
                    new_target_hist.SetBinContent(b, new_content)
                    new_target_hist.SetBinError(b, new_error)
            if THROW_HIST_2D:
                nbins = (hist_2d.GetNbinsX() + 2) * (hist_2d.GetNbinsY() + 2)
                for b in range(1, nbins + 1):
                    old_content = hist_2d.GetBinContent(b)
                    if old_content == 0: continue # A small speedup by skipping over/underflow bins
                    if i == 0:
                        new_content = old_content
                    else:
                        new_content = poisson(old_content)
                    new_error = math.sqrt(new_content)
                    new_hist_2d.SetBinContent(b, new_content)
                    new_hist_2d.SetBinError(b, new_error)
            tries += 1
            foobar = main(new_hist_2d, new_target_hist, 0 / 1000.0)
            if foobar == None:
                fails += 1
            if foobar != None:
                successes += 1
                parameters, parameter_errors, chi2, reg = foobar
                output[i] = parameters, parameter_errors, chi2, reg
                for p, value in enumerate(parameters):
                    param_hist_2d.Fill(p, value)
                i += 1
            
            
    for r in sorted(output):
        print r, 
        for foo in output[r]:
            print foo,
        print
    
    parameters, parameter_errors, chi2, reg = output[0]
    param_hist = make_param_hist(parameters, parameter_errors)
    param_hist_2d.Draw("colz")
    param_hist.Draw("e same")
    c.Print("flux_prism_2d_output/unit_test_poisson_params.png")
    
    print "Tries", tries
    print "Successes", successes
    print "Fails", fails
    if tries > 0:
        print "Fail rate:", (fails / float(tries))
'''

def l_curve(input_directory, material):
    
    tf = ROOT.TFile("%s/separated/flux_tracker.root" % input_directory) # file with 2D histogram
    assert not tf.IsZombie(), "Can't read file"
    tf.ls()
    hist_2d = tf.Get("flux_2d")  # name of the 2D histogram!
    assert hist_2d != None, "Can't get hist_2d"
    target_hist = hist_2d.ProjectionX()
    
    tf2 = ROOT.TFile("%s/flux/flux_%s.root" % (input_directory, material))
    assert not tf2.IsZombie(), "Can't read second file"
    tf2.ls()
    flux1 = tf2.Get("flux")
    assert flux1 != None, "Can't read flux1"
    
    tf3 = ROOT.TFile("%s/flux/flux_tracker.root" % input_directory)
    assert not tf3.IsZombie(), "Can't read third file"
    tf3.ls()
    flux2 = tf3.Get("flux")
    assert flux2 != None, "Can't read flux2"
    
    '''playlist = "minervame1d1m1nweightedave"
    frw = PlotUtils.FluxReweighter(14, True, playlist, 2, 1, 500) 
    target_name = None
    if material == "lead_t25": target_name = "lead_t25"
    if material == "iron_t25": target_name = "iron_t25"
    if material == "lead_t15": target_name = "lead_t15"
    if material == "iron_t15": target_name = "iron_t15"
    if material == "water_radius": target_name = "water_radius"
    if material == "water_apothem": target_name = "water_apothem"
    if material == "tracker": target_name = "tracker"
    if material == "carbon": target_name = "carbon"
    assert target_name != None, "Wasn't able to get a target name from material name: %s" % material
    flux1 = frw.GetNuclearTargetFluxReweighted(14, target_name)
    flux2 = frw.GetNuclearTargetFluxReweighted(14, "tracker")'''
    
    rebin_hists = True
    if rebin_hists:
        hist_2d = rebin(hist_2d, True)
        target_hist = rebin(target_hist)
        flux1 = rebin(flux1)
        flux2 = rebin(flux2)
    
    ndf = target_hist.GetNbinsX() - 1
    ratio = make_clean_copy(flux1)
    ratio.Divide(flux1, flux2)
    
    target_hist.Multiply(target_hist, ratio) # Target is in numerator, tracker in denom, so multiply

    regs = []
    chis = []
    
    output = dict()
    #for reg_lambda in [0, 0.001, 0.01, 0.1, 1, 10, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000, 10000]:
    for reg_lambda in [80, 90, 100, 110, 120, 130, 140, 150]:
    #for reg_lambda in [0, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 500, 750, 1000]:

        foo = main(hist_2d, target_hist, reg_lambda / 1000.0, material)
        if foo == None: 
            print "Failed on lambda =", reg_lambda
            continue
        parameters, parameter_errors, chi2, reg = foo
        output[reg_lambda] = chi2, reg
    for r in sorted(output):
        print r, output[r][0], output[r][1]
        regs.append(r)
        chis.append(output[r][0])
    print(regs)
    print (chis)
    print(ndf) # nbins-1
    chi2perNdf = [ ]
    chi2perNdf = [chi2 / float(ndf) for chi2 in chis]
    print(chi2perNdf) 

        
def calculate_chi2_from_one(hist):
    chi2 = 0
    nbins = hist.GetNbinsX()
    for b in range(1, nbins + 1):
        content = hist.GetBinContent(b)
        error = hist.GetBinError(b)
        if error == 0: continue
        chi2_bin = (content - 1)**2 / error**2
        chi2 += chi2_bin
    return chi2, nbins, chi2/float(nbins)
    
def vector_as_array(v):
    import array
    return array.array("d", v)
    
def rebin(hist, is2d = False):
    #bins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 30, 40, 50, 75, 100]
    if False:
        dn = 4
        bins = [i / float(dn) for i in range(20 * dn + 1)]
        start = bins[-1]
        assert start == 20
        bins.extend([i * 10 + start for i in range(9)])
        assert bins[-1] == 100, "%s" % bins[-1]
    else:
        bins = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 100]
        bins = [0, 4, 6, 8, 10, 15, 20]
    binarray = vector_as_array(bins)
    nbins = len(binarray) - 1
    if is2d:
        binsy = [i for i in range(hist.GetNbinsY() + 1)]
        binarrayy = vector_as_array(binsy)
        nbinsy = len(binarrayy) - 1
        out = PlotUtils.MnvH2D(hist.GetName() + "_rebinned", hist.GetTitle() + "_rebinned", nbins, binarray, nbinsy, binarrayy)
        
        assert hist.GetNbinsY() == out.GetNbinsY()
    else:
        out = PlotUtils.MnvH1D(hist.GetName() + "_rebinned", hist.GetTitle() + "_rebinned", nbins, binarray)
        
    for b in range(1, hist.GetNbinsX() + 1):
        if is2d:
            bx = b
            for by in range(1, hist.GetNbinsY() + 1):
                b = hist.GetBin(bx, by)
                content = hist.GetBinContent(b)
                error = hist.GetBinError(b)
                
                
                target_value_x = hist.GetXaxis().GetBinCenter(bx)
                target_bin_x = out.GetXaxis().FindBin(target_value_x)
                target_value_y = hist.GetYaxis().GetBinCenter(by)
                target_bin_y = out.GetYaxis().FindBin(target_value_y)
                target_bin = out.GetBin(target_bin_x, target_bin_y)
                
                old_content = out.GetBinContent(target_bin)
                old_error = out.GetBinError(target_bin)
                new_content = old_content + content
                new_error = math.sqrt(old_error**2 + error**2)
                out.SetBinContent(target_bin, new_content)
                out.SetBinError(target_bin, new_error)
        else:
            content = hist.GetBinContent(b)
            error = hist.GetBinError(b)
            target_value = hist.GetBinCenter(b)
            target_bin = out.FindBin(target_value)
            old_content = out.GetBinContent(target_bin)
            old_error = out.GetBinError(target_bin)
            new_content = old_content + content
            new_error = math.sqrt(old_error**2 + error**2)
            out.SetBinContent(target_bin, new_content)
            out.SetBinError(target_bin, new_error)
    print "Old Integral:", hist.Integral(), hist, is2d
    print "New Integral:", out.Integral()
    #assert abs(hist.Integral() - out.Integral()) < 0.01
    return out
    
def final(input_directory, material, reg_lambda = 0):

    # Anezka 2022 comments
    # some Assumptions about the input directory
    # two separate directories: separated (for the root file with 2D histogram), flux (with calculated fluxes)
    
    name = "{}_{:06d}".format(material, int(reg_lambda))
    print name
    
    tf = ROOT.TFile("%s/separated/flux_tracker.root" % input_directory) # file with 2D histogram
    assert not tf.IsZombie(), "Can't read file"
    tf.ls()
    hist_2d = tf.Get("flux_2d")  # name of the 2D histogram!
    assert hist_2d != None, "Can't get hist_2d"
    target_hist = hist_2d.ProjectionX()
    
    tf2 = ROOT.TFile("%s/flux/flux_%s.root" % (input_directory, material))
    assert not tf2.IsZombie(), "Can't read second file"
    tf2.ls()
    flux1 = tf2.Get("flux")
    assert flux1 != None, "Can't read flux1"
    
    tf3 = ROOT.TFile("%s/flux/flux_tracker.root" % input_directory)
    assert not tf3.IsZombie(), "Can't read third file"
    tf3.ls()
    flux2 = tf3.Get("flux")
    assert flux2 != None, "Can't read flux2"
    
    '''playlist = "minervame1d1m1nweightedave"
    frw = PlotUtils.FluxReweighter(14, True, playlist, 2, 1, 500) 
    target_name = None
    if material == "lead_t25": target_name = "lead_t25"
    if material == "iron_t25": target_name = "iron_t25"
    if material == "lead_t15": target_name = "lead_t15"
    if material == "iron_t15": target_name = "iron_t15"
    if material == "water_radius": target_name = "water_radius"
    if material == "water_apothem": target_name = "water_apothem"
    if material == "tracker": target_name = "tracker"
    if material == "carbon": target_name = "carbon"
    assert target_name != None, "Wasn't able to get a target name from material name: %s" % material
    flux1 = frw.GetNuclearTargetFluxReweighted(14, target_name)
    flux2 = frw.GetNuclearTargetFluxReweighted(14, "tracker")'''
    
    rebin_hists = True
    if rebin_hists:
        hist_2d = rebin(hist_2d, True)
        target_hist = rebin(target_hist)
        flux1 = rebin(flux1)
        flux2 = rebin(flux2)
    

    ratio = make_clean_copy(flux1)
    ratio.Divide(flux1, flux2)
    
    target_hist.Multiply(target_hist, ratio) # Target is in numerator, tracker in denom, so multiply
    '''for b in range(1, target_hist.GetNbinsX() + 1):
        center = target_hist.GetBinCenter(b)
        center_bin = ratio.FindBin(center)
        r = ratio.GetBinContent(center_bin)
        v = target_hist.GetBinContent(b) * r
        target_hist.SetBinContent(b, v)
        print b, center, center_bin, r, v'''
    
    foo = main(hist_2d, target_hist, reg_lambda / 1000.0, material)
    
    if foo == None:
        print "Couldn't fit lambda =", reg_lambda
        return False
    parameters, parameter_errors, chi2, reg = foo
    param_hist = make_param_hist(parameters, parameter_errors)
        
    
    myrand = ROOT.TRandom3(1234)
    poisson = lambda x: myrand.PoissonD(x)
    N_THROWS_MAX = 0
    THROW_TARGET = True
    THROW_HIST_2D = True
    output = dict()
    n_parameters = hist_2d.GetNbinsY()
    max_y = 5
    param_hist_2d = ROOT.TH2D("param_hist", "param_hist;Flux Bin;Weight", n_parameters, 0, n_parameters, 100 * max_y, 0, max_y)
    i = 0
    done = False
    import time
    start_time = time.clock()
    max_time = 1 # Minutes
    end_time = start_time + 60 * max_time
    fails = 0
    tries = 0
    successes = 0
    while not done and False:
        if i > N_THROWS_MAX:
            done = True
        current_time = time.clock()
        if current_time > end_time:
            done = True
        if not done:
            print "\n\n\nOn throw number %s" % i
            new_target_hist = copy(target_hist)
            new_hist_2d = copy(hist_2d)
            if THROW_TARGET:
                for b in range(1, target_hist.GetNbinsX() + 1):
                    old_content = target_hist.GetBinContent(b)
                    if i == 0:
                        new_content = old_content
                    else:
                        new_content = poisson(old_content)
                    # This would be right if I had not scaled the target hist by the ratio
                    #new_error = math.sqrt(new_content)
                    old_error = target_hist.GetBinError(b)
                    scale_factor = 1.0
                    if old_content > 0:
                        scale_factor = new_content / float(old_content)
                    new_error = scale_factor * old_error
                    new_target_hist.SetBinContent(b, new_content)
                    new_target_hist.SetBinError(b, new_error)
            if THROW_HIST_2D:
                nbins = (hist_2d.GetNbinsX() + 2) * (hist_2d.GetNbinsY() + 2)
                for b in range(1, nbins + 1):
                    old_content = hist_2d.GetBinContent(b)
                    if old_content == 0: continue # A small speedup by skipping over/underflow bins
                    if i == 0:
                        new_content = old_content
                    else:
                        new_content = poisson(old_content)
                    new_error = math.sqrt(new_content)
                    new_hist_2d.SetBinContent(b, new_content)
                    new_hist_2d.SetBinError(b, new_error)
            tries += 1
            foobar = main(new_hist_2d, new_target_hist, reg_lambda / 1000., material)
            if foobar == None:
                fails += 1
            if foobar != None:
                successes += 1
                parameters, parameter_errors, chi2, reg = foobar
                output[i] = parameters, parameter_errors, chi2, reg
                for p, value in enumerate(parameters):
                    param_hist_2d.Fill(p, value)
                i += 1
    
    print(output)
    hist = GetResultingHist(param_hist, hist_2d)
    target_hist.SetLineColor(ROOT.kRed)
    target_hist.Draw("hist")
    foo = max(target_hist.GetMaximum(), hist.GetMaximum()) * 1.2
    target_hist.GetYaxis().SetRangeUser(0, foo)
    target_hist.GetXaxis().SetRangeUser(0, 20)
    target_hist.GetXaxis().SetTitle("Neutrino Energy (GeV)")
    target_hist.GetXaxis().CenterTitle()
    target_hist.GetYaxis().SetTitle("#nu / m^{2} / P.O.T/ GeV")
    target_hist.GetYaxis().CenterTitle()
    target_hist.GetYaxis().SetTitleOffset(1.3)

    hist.Draw("same e")
    c.Print("%s/flux_prism_2d_output/final_fluxes_{}.png".format(name) %input_directory)
    
    flux_ratio = make_clean_copy(flux1)
    flux_ratio.Divide(hist, target_hist)
    flux_ratio.GetXaxis().SetRangeUser(0, 20)
    flux_ratio.GetYaxis().SetRangeUser(0.8, 1.2)
    flux_ratio.GetXaxis().SetTitle("Neutrino Energy (GeV)")
    flux_ratio.GetXaxis().CenterTitle()
    flux_ratio.GetYaxis().SetTitle("Tracker/Target")
    flux_ratio.GetYaxis().CenterTitle()
    flux_ratio.GetYaxis().SetTitleOffset(1.25)
    flux_ratio.Draw("e")
    c.Print("%s/flux_prism_2d_output/final_ratio_{}.png".format(name) %input_directory)
    #flux_ratio.GetXaxis().SetRangeUser(0, 100)
    #flux_ratio.Draw("e")
    #c.Print("flux_prism_2d_output/final_ratio_zoom_out_{}.png".format(name))
    
    flux_ratio_no_tune = make_clean_copy(flux1) # material flux
    no_tune = hist_2d.ProjectionX()  # tracker Enu distributiom
    flux_ratio_no_tune.Divide(no_tune, target_hist) 
    flux_ratio_no_tune.GetYaxis().SetRangeUser(0.8, 1.2)
    flux_ratio_no_tune.GetXaxis().SetTitle("Neutrino Energy (GeV)")
    flux_ratio_no_tune.GetXaxis().CenterTitle()
    flux_ratio_no_tune.GetYaxis().SetTitle("Tracker/Target")
    flux_ratio_no_tune.GetYaxis().CenterTitle()
    flux_ratio_no_tune.GetYaxis().SetTitleOffset(1.25)
    flux_ratio_no_tune.Draw("e")
    c.Print("%s/flux_prism_2d_output/final_no_tune_{}.png".format(name) %input_directory)
    
    
    print "Chi2 with tune:", calculate_chi2_from_one(flux_ratio)
    print "Chi2 with no tune:", calculate_chi2_from_one(flux_ratio_no_tune)
        
    
    param_hist_2d.Draw("COLZ")
    param_hist_2d.GetYaxis().SetRangeUser(0, max_y)
    param_hist_2d.GetXaxis().CenterTitle()
    param_hist_2d.GetYaxis().CenterTitle()
    param_hist.Draw("e same")
    param_hist.GetXaxis().CenterTitle()
    param_hist.GetYaxis().CenterTitle()
    c.Print("%s/flux_prism_2d_output/final_parameters_{}.png".format(name) %input_directory)
    
    tf = ROOT.TFile("%s/flux_prism_2d_output/out_{}.root".format(name) %input_directory, "recreate")
    param_hist.Write()
    flux1.Write("flux")
    flux2.Write("unweighted_flux")
    ratio.Write("flux_ratio")
    tf.Close()
    
    return True
    
def measure_max_difference_from_target(target, hist):
    out = 0
    nbins = hist.GetNbinsX()
    for b in range(1, nbins + 1):
        content = hist.GetBinContent(b)
        diff = abs(target.GetBinContent(b) - hist.GetBinContent(b)) / content
        sigma = math.sqrt(target.GetBinError(b)**2 + hist.GetBinError(b)) / content
        diff_above_one_sigma = diff - sigma
        out = max(out, diff_above_one_sigma)
        print b, diff, sigma, diff_above_one_sigma
    return out
    
def measure_chi2(target, hist):
    out = target.Chi2Test(hist, "WW")
    return out, float(target.GetNbinsX())
    
def measure_chi2_alt(target, hist, scale = 1):
    import array
    ndf = array.array("i", [0])
    chi2 = plotter.Chi2DataMC(target, hist, ndf, scale, True, False, True)
    return chi2, float(ndf[0])
    
def find_plus_one_chi2(target, hist, direction = 1):
    target_chi2, ndf = measure_chi2_alt(target, hist)
    target_chi2 += 1
    upper_scale = 100
    lower_scale = 0
    done = False
    while not done:
        scale = (upper_scale + lower_scale) * 0.5
        #new_hist = copy(hist)
        #new_hist.Scale(direction * scale + 1)
        chi2, ndf = measure_chi2_alt(target, hist, direction * scale + 1)
        diff = chi2 - target_chi2
        if abs(diff) < 0.01:
            done = True
        if abs(diff) < 0.1:
            print upper_scale, lower_scale, scale, diff
        if diff > 0:
            upper_scale = scale
        else:
            lower_scale = scale
            
    return scale
    
def make_chi2_plot(target, hist):
    n = 101
    diff = 0.005
    out = ROOT.TH1D("", "", n, -diff, diff)
    for i in range(1, n + 1):
        scale = out.GetBinCenter(i)
        chi2, ndf = measure_chi2_alt(target, hist, scale + 1)
        out.SetBinContent(i, chi2)
    return out

    
def measure_total_difference_from_target(target, hist):
    out = abs(target.Integral() - hist.Integral()) / hist.Integral()
    
    return out
    
'''    
def unit_test_simple_carbon():
    input_directory = "../plots/2021-05-18/daisy_vs_flux/daisy_vs_flux_hd_no_circle/"

    tf = ROOT.TFile("%s/separated/tracker.root" % input_directory)
    assert not tf.IsZombie(), "Can't read file"
    tf.ls()
    hist_2d = tf.Get("truecarbon")
    assert hist_2d != None, "Can't get hist_2d"
    target_hist = hist_2d.ProjectionX()
    
    tf2 = ROOT.TFile("%s/flux_carbon.root" % input_directory)
    assert not tf2.IsZombie(), "Can't read second file"
    tf2.ls()
    flux1 = tf2.Get("flux")
    assert flux1 != None, "Can't read flux1"
    
    tf3 = ROOT.TFile("%s/flux_tracker.root" % input_directory)
    assert not tf3.IsZombie(), "Can't read third file"
    tf3.ls()
    flux2 = tf3.Get("flux")
    assert flux2 != None, "Can't read flux2"
    
    rebin_hists = False
    if rebin_hists:
        hist_2d = rebin(hist_2d, True)
        target_hist = rebin(target_hist)
        flux1 = rebin(flux1)
        flux2 = rebin(flux2)
    
    ratio = make_clean_copy(flux1)
    ratio.Divide(flux1, flux2)
    
    target_hist.Multiply(target_hist, ratio) # Target is in numerator, tracker in denom, so multiply
    
    n_poisson_throws = 0
    
    mode = "fit"
    if mode == "carbon":
        parameters = [2] * hist_2d.GetNbinsY()
        assert len(parameters) == 12
        parameters[5] = 0
        parameters[6] = 0
        parameters[7] = 0
        parameters[8] = 0
        parameters[9] = 0
        parameters[10] = 0
        parameter_errors = [0] * len(parameters)
        param_hist = make_param_hist(parameters, parameter_errors)
    if mode == None:
        parameters = [1] * hist_2d.GetNbinsY()
        assert len(parameters) == 12
        parameter_errors = [0] * len(parameters)
        param_hist = make_param_hist(parameters, parameter_errors)
    if mode == "fit":
        reg_lambda = 100
        foo = main(hist_2d, target_hist, reg_lambda / 1000.0)
        assert foo != None
        parameters, parameter_errors, chi2, reg = foo
        parameter_errors = [0] * len(parameter_errors)
        param_hist = make_param_hist(parameters, parameter_errors)
        
        n_poisson_throws = 2
        param_hist.AddVertErrorBandAndFillWithCV("Daisy Poisson Throws", n_poisson_throws)
        n_throws = 0
        n_fails = 0
        myrand = ROOT.TRandom3(1234)
        poisson = lambda x: myrand.PoissonD(x)
        while n_throws < n_poisson_throws:
            new_target_hist = copy(target_hist)
            new_hist_2d = copy(hist_2d)
            THROW_TARGET = True
            THROW_HIST_2D = True
            if THROW_TARGET:
                for b in range(1, target_hist.GetNbinsX() + 1):
                    old_content = target_hist.GetBinContent(b)
                    new_content = poisson(old_content)
                    # This would be right if I had not scaled the target hist by the ratio
                    #new_error = math.sqrt(new_content)
                    old_error = target_hist.GetBinError(b)
                    scale_factor = 1.0
                    if old_content > 0:
                        scale_factor = new_content / float(old_content)
                    new_error = scale_factor * old_error
                    new_target_hist.SetBinContent(b, new_content)
                    new_target_hist.SetBinError(b, new_error)
            if THROW_HIST_2D:
                nbins = (hist_2d.GetNbinsX() + 2) * (hist_2d.GetNbinsY() + 2)
                for b in range(1, nbins + 1):
                    old_content = hist_2d.GetBinContent(b)
                    if old_content == 0: continue # A small speedup by skipping over/underflow bins
                    new_content = poisson(old_content)
                    new_error = math.sqrt(new_content)
                    new_hist_2d.SetBinContent(b, new_content)
                    new_hist_2d.SetBinError(b, new_error)
        
        
            foo = main(new_hist_2d, new_target_hist, reg_lambda / 1000.0)
            if foo == None:
                n_fails += 1
            else:
                parameters, parameter_errors, chi2, reg = foo
                modify_this_hist = param_hist.GetVertErrorBand("Daisy Poisson Throws").GetHist(n_throws)
                make_param_hist(parameters, parameter_errors, "param_hist", modify_this_hist)
                n_throws += 1
        print "N fails:", n_fails
            
        
    
    
    hist = GetResultingHist(param_hist, hist_2d)
    flux_ratio = make_clean_copy(flux1)
    flux_ratio.Divide(hist, target_hist)
    flux_ratio.GetXaxis().SetRangeUser(0, 20)
    flux_ratio.GetYaxis().SetRangeUser(0.8, 1.2)
    flux_ratio.Draw("e")
    c.Print("flux_prism_2d_output/unit_test_simple_carbon_ratio.png")
    
    param_hist.AddVertErrorBandAndFillWithCV("Daisy Total Diff", 2)
    #param_hist.AddVertErrorBandAndFillWithCV("Daisy Max Diff", 2)
    param_hist.AddVertErrorBandAndFillWithCV("Daisy Scale +/- 1 Chi2", 2)
    param_hist.AddVertErrorBandAndFillWithCV("Daisy Cut Eff", 2)
    
    print "Max diff"
    max_diff = measure_max_difference_from_target(target_hist, hist)
    print max_diff
    print "Total diff"
    total_diff = measure_total_difference_from_target(target_hist, hist)
    print total_diff
    print "simple chi2"
    print measure_chi2(target_hist, hist)
    print "alt chi2"
    print measure_chi2_alt(target_hist, hist)
    print "find_plus_one_chi2"
    chi2_plus_one = find_plus_one_chi2(target_hist, hist)
    print chi2_plus_one
    print "find_plus_one_chi2 alternate"
    chi2_minus_one = find_plus_one_chi2(target_hist, hist, -1)
    print chi2_minus_one
    chi2_plot = make_chi2_plot(target_hist, hist)
    chi2_plot.Draw("hist")
    c.Print("flux_prism_2d_output/unit_test_simple_carbon_chi2_plot.png")
    
    vert = param_hist.GetVertErrorBand("Daisy Scale +/- 1 Chi2")
    #vert.GetHist(0).Scale(1 - chi2_minus_one)
    #vert.GetHist(1).Scale(1 + chi2_plus_one)
    chi2_diff = max(abs(chi2_plus_one), abs(chi2_minus_one))
    vert.GetHist(0).Scale(1 - chi2_diff)
    vert.GetHist(1).Scale(1 + chi2_diff)
    
    
    vert = param_hist.GetVertErrorBand("Daisy Total Diff")
    vert.GetHist(0).Scale(1 - total_diff)
    vert.GetHist(1).Scale(1 + total_diff)
    
    #vert = param_hist.GetVertErrorBand("Daisy Max Diff")
    #vert.GetHist(0).Scale(1 - max_diff)
    #vert.GetHist(1).Scale(1 + max_diff)
    
    # Data / MC cut effs as a percent
    # Where cut eff = N events in daisy bin i / N events total
    cut_eff = [0.9, -3.3, -1.8, -2.4, -4.3, -0.01, 0.7, 1.2, 1.4, 4.1, 1.3, 2.5]
    for i in range(len(cut_eff)):
        b = i + 1
        cut_eff_scaled = cut_eff[i] / 100.0
        vert = param_hist.GetVertErrorBand("Daisy Cut Eff")
        content = param_hist.GetBinContent(b)
        vert.GetHist(0).SetBinContent(b, content * (1 - cut_eff_scaled))
        vert.GetHist(1).SetBinContent(b, content * (1 + cut_eff_scaled))
        
    
    import error_band_by_name
    
    hist = GetResultingHistFull(param_hist, hist_2d)
    flux_ratio = make_clean_copy(flux1)
    flux_ratio.AddMissingErrorBandsAndFillWithCV(hist)
    target_hist.AddMissingErrorBandsAndFillWithCV(hist)
    print "\n\n\nFixing"
    #flux_ratio.PopVertErrorBand("NuMI Flux")
    flux_ratio.PopVertErrorBand("Flux")
    print "Done fixing\n\n\n"
    flux_ratio.Divide(hist, target_hist)
    flux_ratio.GetXaxis().SetRangeUser(0, 20)
    flux_ratio.GetYaxis().SetRangeUser(0.8, 1.2)
    flux_ratio.Draw("e1")
    flux_ratio_errors = flux_ratio.GetCVHistoWithError(True, False)
    flux_ratio_errors.Draw("e1 same")
    c.Print("flux_prism_2d_output/unit_test_simple_carbon_ratio_full.png")
    error_band_by_name.draw_error_band(flux_ratio)
    c.Print("flux_prism_2d_output/unit_test_simple_carbon_ratio_full_errors.png")
    
    error_band_by_name.draw_error_band(param_hist)
    c.Print("flux_prism_2d_output/unit_test_simple_carbon_errors.png")
    
    flux_ratio_no_tune = make_clean_copy(flux1)
    no_tune = hist_2d.ProjectionX()
    flux_ratio_no_tune.Divide(no_tune, target_hist)
    flux_ratio_no_tune.GetXaxis().SetRangeUser(0, 20)
    flux_ratio_no_tune.GetYaxis().SetRangeUser(0.8, 1.2)
    flux_ratio_no_tune.Draw("e")
    c.Print("flux_prism_2d_output/unit_test_simple_carbon_no_tune.png")
    
    
    print "Chi2 with tune:", calculate_chi2_from_one(flux_ratio)
    print "Chi2 with no tune:", calculate_chi2_from_one(flux_ratio_no_tune)
    
    param_hist.GetYaxis().SetRangeUser(0, 2.5)
    param_hist.Draw("e")
    c.Print("flux_prism_2d_output/unit_test_simple_carbon_parameters.png")
    
    tfout = ROOT.TFile("flux_prism_2d_output/unit_test_simple_carbon_parameters.root", "recreate")
    tfout.cd()
    flux_ratio_no_tune.Write("flux_ratio_no_tune")
    flux_ratio.Write("flux_ratio")
    hist.Write("flux_hist")
    target_hist.Write("target_hist")
    flux1.Write("flux_initial")
    flux2.Write("flux_target")
    hist_2d.Write("hist_2d")
    param_hist.Write("param_hist")
    tfout.Close()
'''  
    
#unit_test()
# 0, 500
#for reg_lambda in [0, 10, 100, 750]:
#    unit_test_with_target(reg_lambda)
    
#unit_test_l_curve()
#unit_test_poisson()
#unit_test_simple_carbon()
#exit()

if __name__ == "__main__":
    '''
    material = sys.argv[1]
    l_curve("/minerva/data2/users/anezkak/fluxFinal/Antineutrinos", material)
    exit()
    '''
    filename = sys.argv[1]
    material = sys.argv[2]
    reg_lambda = None
    if len(sys.argv) > 3:
        reg_lambda = int(sys.argv[3])
    if reg_lambda == None:
      couldnt_fit = []
      for reg_lambda in [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 300, 500]:
          print "\n\n\n\n\nOn", reg_lambda
          foo = final(filename, material, reg_lambda)
          if not foo:
              print "\n\n\nCouldn't fit"
              couldnt_fit.append(reg_lambda)
      print "The following didn't fit for some reason", couldnt_fit
      exit()
    final(filename, material, reg_lambda)
    
    
    
    
    
    
    
    
    
    
    
    
    
    


