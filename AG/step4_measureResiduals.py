import ROOT
import sys
import os
import argparse
import json
from math import sqrt
import math
import array
import numpy as np
parser=argparse.ArgumentParser()
parser.add_argument("-d",  "--inputdir",  default="./",  help="Directory where derived corrections with type 2 only are located.")
parser.add_argument("-o",  "--outputdir", default="brilhist/", help="The output directory")
parser.add_argument("-f",  "--file", default="",  help="File of derived corrections is located.")
parser.add_argument("-l",  "--label",default="",  help="Append the names of output files with this label.")
parser.add_argument("-j",  "--json", default="",  help="Certification JSON file for selecting runs.")
args=parser.parse_args()
verbose=1

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)

BXLength=3564

fill2lumi={}
with open('supertable_22.txt') as f:
    lines = f.readlines()
    for iline in lines:
        #print(iline.split("\t"))
        fill2lumi[iline.split("\t")[0]]=(iline.split("\t")[1])

if args.outputdir is not None :
    if not os.path.isdir( args.outputdir ) :
        os.makedirs( args.outputdir )

filenames=[]
if args.inputdir!="":
    shortfilenames=os.listdir(args.inputdir)
    for shortfilename in shortfilenames:
        if shortfilename.find("Overall_")!=-1:
            filenames.append(args.inputdir+"/"+shortfilename)
if args.file!="":
    filenames.append(args.file)
print('filenames:',filenames)

hists={}
if args.label!="":
    outRootFileName="systematicHistograms_"+args.label+".root"
else:
    outRootFileName="systematicHistograms_"+".root"


oFilehistograms=ROOT.TFile(args.outputdir+"/"+outRootFileName,"RECREATE")
newfile=ROOT.TFile(args.outputdir+"/type2minitree_"+args.label+".root", "RECREATE")
newtree=ROOT.TTree("minitree", "minitree")

FILL = array.array('I', [0])
RUN = array.array( 'I', [0])
LS  = array.array( 'I', [0])
NActive  = array.array( 'I', [0])

sum_lead     = array.array( 'd', [0])
sum_train     = array.array( 'd', [0])
AveSBIL_leadBX     = array.array( 'd', [0])

# raw
AveSBIL_raw     = array.array( 'd', [0])
overallLumi_raw = array.array( 'd', [0])
type1frac_raw   = array.array( 'd', [0])
type1abs_raw    = array.array( 'd', [0])
type1frac_raw_e   = array.array( 'd', [0])
type1abs_raw_e    = array.array( 'd', [0])
type2frac_raw   = array.array( 'd', [0])
type2abs_raw    = array.array( 'd', [0])
type2fracRMS_raw   = array.array( 'd', [0])
# raw - noise
AveSBIL_rawSubnoise     = array.array( 'd', [0])
overallLumi_rawSubnoise = array.array( 'd', [0])
type1frac_rawSubnoise   = array.array( 'd', [0])
type1abs_rawSubnoise    = array.array( 'd', [0])
type1frac_rawSubnoise_e   = array.array( 'd', [0])
type1abs_rawSubnoise_e    = array.array( 'd', [0])
type2frac_rawSubnoise   = array.array( 'd', [0])
type2abs_rawSubnoise    = array.array( 'd', [0])
type2fracRMS_rawSubnoise   = array.array( 'd', [0])
# raw - localtype1
AveSBIL_rawSublocalT1     = array.array( 'd', [0])
overallLumi_rawSublocalT1 = array.array( 'd', [0])
type1frac_rawSublocalT1   = array.array( 'd', [0])
type1abs_rawSublocalT1    = array.array( 'd', [0])
type1frac_rawSublocalT1_e   = array.array( 'd', [0])
type1abs_rawSublocalT1_e    = array.array( 'd', [0])
type2frac_rawSublocalT1   = array.array( 'd', [0])
type2abs_rawSublocalT1    = array.array( 'd', [0])
type2fracRMS_rawSublocalT1   = array.array( 'd', [0])


newtree.Branch("FILL", FILL, "FILL/I")
newtree.Branch("RUN", RUN, "RUN/I")
newtree.Branch("LS", LS, "LS/I")
newtree.Branch("NActive", NActive, "NActive/I")
newtree.Branch("sum_lead", sum_lead, "sum_lead/D")
newtree.Branch("sum_train",sum_train , "sum_train/D")
newtree.Branch("AveSBIL_leadBX", AveSBIL_leadBX, "AveSBIL_leadBX/D")
# raw - localtype1
newtree.Branch("AveSBIL_rawSublocalT1", AveSBIL_rawSublocalT1, "AveSBIL_rawSublocalT1/D")
newtree.Branch("overallLumi_rawSublocalT1",overallLumi_rawSublocalT1 , "overallLumi_rawSublocalT1/D")
newtree.Branch("type1frac_rawSublocalT1", type1frac_rawSublocalT1, "type1frac_rawSublocalT1/D")
newtree.Branch("type1abs_rawSublocalT1", type1abs_rawSublocalT1, "type1abs_rawSublocalT1/D")
newtree.Branch("type1frac_rawSublocalT1_e", type1frac_rawSublocalT1_e, "type1frac_rawSublocalT1_e/D")
newtree.Branch("type1abs_rawSublocalT1_e", type1abs_rawSublocalT1_e, "type1abs_rawSublocalT1_e/D")
newtree.Branch("type2frac_rawSublocalT1", type2frac_rawSublocalT1, "type2frac_rawSublocalT1/D")
newtree.Branch("type2abs_rawSublocalT1", type2abs_rawSublocalT1, "type2abs_rawSublocalT1/D")
newtree.Branch("type2fracRMS_rawSublocalT1", type2fracRMS_rawSublocalT1, "type2fracRMS_rawSublocalT1/D")
# raw - noise
newtree.Branch("AveSBIL_rawSubnoise", AveSBIL_rawSubnoise, "AveSBIL_rawSubnoise/D")
newtree.Branch("overallLumi_rawSubnoise",overallLumi_rawSubnoise , "overallLumi_rawSubnoise/D")
newtree.Branch("type1frac_rawSubnoise", type1frac_rawSubnoise, "type1frac_rawSubnoise/D")
newtree.Branch("type1abs_rawSubnoise", type1abs_rawSubnoise, "type1abs_rawSubnoise/D")
newtree.Branch("type1frac_rawSubnoise_e", type1frac_rawSubnoise_e, "type1frac_rawSubnoise_e/D")
newtree.Branch("type1abs_rawSubnoise_e", type1abs_rawSubnoise_e, "type1abs_rawSubnoise_e/D")
newtree.Branch("type2frac_rawSubnoise", type2frac_rawSubnoise, "type2frac_rawSubnoise/D")
newtree.Branch("type2abs_rawSubnoise", type2abs_rawSubnoise, "type2abs_rawSubnoise/D")
newtree.Branch("type2fracRMS_rawSubnoise", type2fracRMS_rawSubnoise, "type2fracRMS_rawSubnoise/D")
# raw
newtree.Branch("AveSBIL_raw", AveSBIL_raw, "AveSBIL_raw/D")
newtree.Branch("overallLumi_raw",overallLumi_raw , "overallLumi_raw/D")
newtree.Branch("type1frac_raw", type1frac_raw, "type1frac_raw/D")
newtree.Branch("type1abs_raw", type1abs_raw, "type1abs_raw/D")
newtree.Branch("type1frac_raw_e", type1frac_raw_e, "type1frac_raw_e/D")
newtree.Branch("type1abs_raw_e", type1abs_raw_e, "type1abs_raw_e/D")
newtree.Branch("type2frac_raw", type2frac_raw, "type2frac_raw/D")
newtree.Branch("type2abs_raw", type2abs_raw, "type2abs_raw/D")
newtree.Branch("type2fracRMS_raw", type2fracRMS_raw, "type2fracRMS_raw/D")

LB2Fill={}

# Make some plots
# Raw t1 residual (frac) vs SBIL   -> Profile/Per-fill
tgraph_t1_frac_sbil=ROOT.TGraphErrors()
tgraph_t1_frac_sbil.SetTitle("type1FracVsSBIL;Average SBIL (Hz/#muB);Type 1 residual (Fraction)")
thist_t1_frac_sbil={}
thist_lumi_correction={}
thist_lumi_correction['total']=ROOT.TH2F("thist_lumi_correction","Overall Correction;Average SBIL (Hz/#muB);Correction Factor",100,0,10,100,0.9,1.1)
# Raw t1 residual (frac) vs Lumi
tgraph_t1_frac_lumi=ROOT.TGraphErrors()
tgraph_t1_frac_lumi.SetTitle("type1FracVsLumi;Lumi(fb^{-1});Type 1 residual (Fraction)")
# Raw t1 residual (abs) vs SBIL    -> Profile/Per-fill
tgraph_t1_abs_sbil=ROOT.TGraphErrors()
tgraph_t1_abs_sbil.SetTitle("type1FracVsSBIL;Average SBIL (Hz/#muB);Type 1 residual (SBIL,Hz/#mub)")
thist_t1_abs_sbil={}
# Raw t1 residual (abs) vs Lumi
tgraph_t1_abs_lumi=ROOT.TGraphErrors()
tgraph_t1_abs_lumi.SetTitle("type1FracVsLumi;Lumi(fb^{-1});Type 1 residual (SBIL,Hz/#mub)")

tgraph_t2_frac_sbil=ROOT.TGraphErrors()
tgraph_t2_frac_sbil.SetTitle("type2FracVsSBIL;Average SBIL (Hz/#muB);Type 2 residual (Fraction)")
tgraph_t2_frac_lumi=ROOT.TGraphErrors()
tgraph_t2_frac_lumi.SetTitle("type2FracVsLumi;Lumi(fb^{-1});Type 2 residual (Fraction)")
tgraph_t2_abs_sbil=ROOT.TGraphErrors()
tgraph_t2_abs_sbil.SetTitle("type2FracVsSBIL;Average SBIL (Hz/#muB);Type 2 residual (SBIL,Hz/#mub)")
tgraph_t2_abs_lumi=ROOT.TGraphErrors()
tgraph_t2_abs_lumi.SetTitle("type2FracVsLumi;Lumi(fb^{-1});Type 2 residual (SBIL,Hz/#mub)")
tgraph_t2_rms_sbil=ROOT.TGraphErrors()
tgraph_t2_rms_sbil.SetTitle("type2FracVsSBIL;Average SBIL (Hz/#muB);Type 2 residual RMS (Fraction)")
tgraph_t2_rms_lumi=ROOT.TGraphErrors()
tgraph_t2_rms_lumi.SetTitle("type2FracVsLumi;Lumi(fb^{-1});Type 2 residual RMS (Fraction)")
thist_t2_frac_sbil={}
thist_t2_abs_sbil={}
thist_t2_rms_sbil={}

# Overall correction vs SBIL   -> Profile/Per-fill   -> Fit with pol1
overallcor_sbil={}
# po and p1 vs Fill, Lumi
nevent=0
lumitotal=0
for filename in filenames:
    try:
        tfile=ROOT.TFile.Open(filename)
        if verbose>10: print("Open ", filename)
        f1keys=tfile.GetListOfKeys()
    except:
        continue
    fHistNames=[]
    hists={}
    for f1key in f1keys:
        if f1key.GetName().find("hist_info_")!=-1:
            fHistNames.append(f1key.GetName())
    if verbose>10: print(fHistNames)
    fHistNames.sort()
    for fHistName in fHistNames:
        #if int(fHistName.split("_")[5].lstrip("Fill"))!=8313: continue
        #if int(fHistName.split("_")[5].lstrip("Fill")) >8313: exit()
        tfile.cd()
        thisLB=fHistName.split("_")[2]+"_"+fHistName.split("_")[3]+"_"+fHistName.split("_")[4]
        allLumiPerBX_oldSBR = "allLumiPerBX_oldSBR_"+fHistName.split("_")[2]+"_"+fHistName.split("_")[3]+"_"+fHistName.split("_")[4]+"_"+fHistName.split("_")[5]
        allLumiPerBX_oldSBR_cor = "allLumiPerBX_oldSBR_cor_"+fHistName.split("_")[2]+"_"+fHistName.split("_")[3]+"_"+fHistName.split("_")[4]+"_"+fHistName.split("_")[5]
        hist_info = "hist_info_"+fHistName.split("_")[2]+"_"+fHistName.split("_")[3]+"_"+fHistName.split("_")[4]+"_"+fHistName.split("_")[5]

        LB2Fill[thisLB] = int(fHistName.split("_")[5].lstrip("Fill"))
        if verbose>10: print("thisLB, LB2Fill[thisLB] ", thisLB, LB2Fill[thisLB])
        hists[allLumiPerBX_oldSBR]=tfile.Get(allLumiPerBX_oldSBR)
        hists[allLumiPerBX_oldSBR_cor]=tfile.Get(allLumiPerBX_oldSBR_cor)
        hists[hist_info]=tfile.Get(hist_info)
        # Raw-noise, raw, raw-noise-t1
        # Average SBIL, Total, type1 fraction, type1 Hz/ub
        FILL[0] = LB2Fill[thisLB]
        RUN[0] = int(thisLB.split("_")[0])
        LS[0] = int(thisLB.split("_")[1].lstrip("LS"))
        NActive[0] = int(hists[hist_info][22])
        #print(FILL[0], RUN[0],  LS[0])
        AveSBIL_rawSubnoise[0]        = hists[hist_info][1]
        overallLumi_rawSubnoise[0]    = hists[hist_info][2] /25.0
        type1frac_rawSubnoise[0]      = hists[hist_info][3]
        type1abs_rawSubnoise[0]       = hists[hist_info][4]
        type1frac_rawSubnoise_e[0]      = hists[hist_info].GetBinError(3)
        type1abs_rawSubnoise_e[0]       = hists[hist_info].GetBinError(4)
        type2frac_rawSubnoise[0]      = hists[hist_info][5]
        type2abs_rawSubnoise[0]       = hists[hist_info][6]
        type2fracRMS_rawSubnoise[0]   = hists[hist_info][7]

        AveSBIL_raw[0]        = hists[hist_info][8]
        overallLumi_raw[0]    = hists[hist_info][9]/25.0
        type1frac_raw[0]      = hists[hist_info][10]
        type1abs_raw[0]       = hists[hist_info][11]
        type1frac_raw_e[0]      = hists[hist_info].GetBinError(10)
        type1abs_raw_e[0]       = hists[hist_info].GetBinError(11)
        type2frac_raw[0]      = hists[hist_info][12]
        type2abs_raw[0]       = hists[hist_info][13]
        type2fracRMS_raw[0]   = hists[hist_info][14]

        AveSBIL_rawSublocalT1[0]        = hists[hist_info][15]
        overallLumi_rawSublocalT1[0]    = hists[hist_info][16]/25.0
        type1frac_rawSublocalT1[0]      = hists[hist_info][17]
        type1abs_rawSublocalT1[0]       = hists[hist_info][18]
        type1frac_rawSublocalT1_e[0]      = hists[hist_info].GetBinError(17)
        type1abs_rawSublocalT1_e[0]       = hists[hist_info].GetBinError(18)
        type2frac_rawSublocalT1[0]      = hists[hist_info][19]
        type2abs_rawSublocalT1[0]       = hists[hist_info][20]
        type2fracRMS_rawSublocalT1[0]   = hists[hist_info][21]

        sum_lead[0]   = hists[hist_info][23]
        sum_train[0]   = hists[hist_info][24]
        AveSBIL_leadBX[0]   = hists[hist_info][25]

        newtree.Fill()
        if math.isnan(overallLumi_rawSublocalT1[0]): 
            print('nan ',FILL[0],RUN[0],LS[0],NActive[0])
        if not math.isnan(overallLumi_rawSublocalT1[0]): lumitotal+=overallLumi_rawSublocalT1[0]*1e-9
        # Use integrated Lumi
        if str(LB2Fill[thisLB]) in fill2lumi.keys(): lumitotal=1e-6*float(fill2lumi[str(LB2Fill[thisLB])])
        else:
            print(thisLB," not in fill2lumi" )
            exit()
        tgraph_t1_frac_sbil.SetPoint(nevent,AveSBIL_raw[0],type1frac_raw[0])
        tgraph_t1_frac_sbil.SetPointError(nevent,0,type1frac_raw_e[0])
        tgraph_t1_frac_lumi.SetPoint(nevent,lumitotal,type1frac_raw[0])
        tgraph_t1_frac_lumi.SetPointError(nevent,0,type1frac_raw_e[0])
        tgraph_t1_abs_sbil.SetPoint(nevent,AveSBIL_raw[0],type1abs_raw[0])
        tgraph_t1_abs_sbil.SetPointError(nevent,0,type1abs_raw_e[0])
        tgraph_t1_abs_lumi.SetPoint(nevent,lumitotal,type1abs_raw[0])
        tgraph_t1_abs_lumi.SetPointError(nevent,0,type1abs_raw_e[0])
        if FILL[0] not in thist_t1_frac_sbil.keys(): thist_t1_frac_sbil[FILL[0]]=[ROOT.TGraphErrors(),0]
        thist_t1_frac_sbil[FILL[0]][0].SetPoint(thist_t1_frac_sbil[FILL[0]][1],type1abs_rawSubnoise[0],type1frac_rawSubnoise[0])
        thist_t1_frac_sbil[FILL[0]][0].SetPointError(thist_t1_frac_sbil[FILL[0]][1],type1abs_rawSubnoise_e[0],type1frac_rawSubnoise_e[0])
        thist_t1_frac_sbil[FILL[0]][1]+=1
        if FILL[0] not in thist_lumi_correction.keys(): thist_lumi_correction[FILL[0]]=[ROOT.TGraphErrors(),0,0,0,0]
        thist_lumi_correction[FILL[0]][0].SetPoint(thist_lumi_correction[FILL[0]][1],AveSBIL_raw[0], AveSBIL_rawSublocalT1[0]/AveSBIL_raw[0])
        thist_lumi_correction[FILL[0]][2] += overallLumi_rawSublocalT1[0]
        thist_lumi_correction[FILL[0]][3] += overallLumi_raw[0]
        thist_lumi_correction[FILL[0]][4] = NActive[0]
        thist_lumi_correction[FILL[0]][1]+=1
        thist_lumi_correction['total'].Fill(AveSBIL_raw[0], AveSBIL_rawSublocalT1[0]/AveSBIL_raw[0])
        tgraph_t2_frac_sbil.SetPoint(nevent,AveSBIL_raw[0],type2frac_raw[0])
        tgraph_t2_frac_lumi.SetPoint(nevent,lumitotal,type2frac_raw[0])
        tgraph_t2_abs_sbil.SetPoint(nevent,AveSBIL_raw[0],type2abs_raw[0])
        tgraph_t2_abs_lumi.SetPoint(nevent,lumitotal,type2abs_raw[0])
        tgraph_t2_rms_sbil.SetPoint(nevent,AveSBIL_raw[0],type2fracRMS_raw[0])
        tgraph_t2_rms_lumi.SetPoint(nevent,lumitotal,type2fracRMS_raw[0])
        nevent+=1

newfile.WriteTObject(newtree, "newtree")

text=ROOT.TLatex(0.72,0.92,"2022  (13.6 TeV)")
text.SetNDC()
text.SetTextFont(62)
text.SetTextSize(0.05)
text2=ROOT.TLatex(0.15,0.92,"CMS #bf{#scale[0.75]{#it{Preliminary}}}")
text2.SetNDC()
text2.SetTextSize(0.05)
text2.SetTextFont(62)

skip_someplots=True
if not skip_someplots:
    tgraph_t1_frac_sbil.GetYaxis().SetTitleOffset(1.5)
    tgraph_t1_frac_sbil.SetMarkerStyle(23)
    tgraph_t1_frac_sbil.SetMarkerSize(0.5)
    tgraph_t1_frac_sbil.SetMarkerColor(ROOT.kBlue)
    tgraph_t1_frac_sbil.SetLineColor(ROOT.kBlue)
    tgraph_t1_frac_lumi.GetYaxis().SetTitleOffset(1.5)
    tgraph_t1_frac_lumi.SetMarkerStyle(23)
    tgraph_t1_frac_lumi.SetMarkerSize(0.5)
    tgraph_t1_frac_lumi.SetMarkerColor(ROOT.kBlue)
    tgraph_t1_frac_lumi.SetLineColor(ROOT.kBlue)
    tgraph_t1_abs_sbil.GetYaxis().SetTitleOffset(1.5)
    tgraph_t1_abs_sbil.SetMarkerStyle(23)
    tgraph_t1_abs_sbil.SetMarkerSize(0.5)
    tgraph_t1_abs_sbil.SetMarkerColor(ROOT.kBlue)
    tgraph_t1_abs_sbil.SetLineColor(ROOT.kBlue)
    tgraph_t1_abs_lumi.GetYaxis().SetTitleOffset(1.5)
    tgraph_t1_abs_lumi.SetMarkerStyle(23)
    tgraph_t1_abs_lumi.SetMarkerSize(0.5)
    tgraph_t1_abs_lumi.SetMarkerColor(ROOT.kBlue)
    tgraph_t1_abs_lumi.SetLineColor(ROOT.kBlue)
    
    can=ROOT.TCanvas("can","can",800,800)
    can.Divide(2,2)
    can.cd(1)
    tgraph_t1_frac_sbil.Draw("AP")
    tgraph_t1_frac_sbil.GetYaxis().SetRangeUser(-0.1, 0.1)
    can.cd(2)
    tgraph_t1_frac_lumi.Draw("AP")
    tgraph_t1_frac_lumi.GetYaxis().SetRangeUser(-0.1, 0.1)
    can.cd(3)
    tgraph_t1_abs_sbil.Draw("AP")
    tgraph_t1_abs_sbil.GetYaxis().SetRangeUser(-0.1, 0.1)
    can.cd(4)
    tgraph_t1_abs_lumi.Draw("AP")
    tgraph_t1_abs_lumi.GetYaxis().SetRangeUser(-0.1, 0.1)
    can.Update()
    can.SaveAs(args.outputdir+"/type1_raw"+args.label+".png")
    can.SaveAs(args.outputdir+"/type1_raw"+args.label+".C")
    
    
    tgraph_t2_frac_sbil.GetYaxis().SetTitleOffset(1.5)
    tgraph_t2_frac_sbil.SetMarkerStyle(23)
    tgraph_t2_frac_sbil.SetMarkerSize(0.5)
    tgraph_t2_frac_sbil.SetMarkerColor(ROOT.kBlue)
    tgraph_t2_frac_sbil.SetLineColor(ROOT.kBlue)
    tgraph_t2_frac_lumi.GetYaxis().SetTitleOffset(1.5)
    tgraph_t2_frac_lumi.SetMarkerStyle(23)
    tgraph_t2_frac_lumi.SetMarkerSize(0.5)
    tgraph_t2_frac_lumi.SetMarkerColor(ROOT.kBlue)
    tgraph_t2_frac_lumi.SetLineColor(ROOT.kBlue)
    tgraph_t2_abs_sbil.GetYaxis().SetTitleOffset(1.5)
    tgraph_t2_abs_sbil.SetMarkerStyle(23)
    tgraph_t2_abs_sbil.SetMarkerSize(0.5)
    tgraph_t2_abs_sbil.SetMarkerColor(ROOT.kBlue)
    tgraph_t2_abs_sbil.SetLineColor(ROOT.kBlue)
    tgraph_t2_abs_lumi.GetYaxis().SetTitleOffset(1.5)
    tgraph_t2_abs_lumi.SetMarkerStyle(23)
    tgraph_t2_abs_lumi.SetMarkerSize(0.5)
    tgraph_t2_abs_lumi.SetMarkerColor(ROOT.kBlue)
    tgraph_t2_abs_lumi.SetLineColor(ROOT.kBlue)
    tgraph_t2_rms_sbil.GetYaxis().SetTitleOffset(1.5)
    tgraph_t2_rms_sbil.SetMarkerStyle(23)
    tgraph_t2_rms_sbil.SetMarkerSize(0.5)
    tgraph_t2_rms_sbil.SetMarkerColor(ROOT.kBlue)
    tgraph_t2_rms_sbil.SetLineColor(ROOT.kBlue)
    tgraph_t2_rms_lumi.GetYaxis().SetTitleOffset(1.5)
    tgraph_t2_rms_lumi.SetMarkerStyle(23)
    tgraph_t2_rms_lumi.SetMarkerSize(0.5)
    tgraph_t2_rms_lumi.SetMarkerColor(ROOT.kBlue)
    tgraph_t2_rms_lumi.SetLineColor(ROOT.kBlue)
    
    can2=ROOT.TCanvas("can2","can2",1200,800)
    can2.Divide(3,2)
    can2.cd(1)
    tgraph_t2_frac_sbil.Draw("AP")
    tgraph_t2_frac_sbil.GetYaxis().SetRangeUser(-0.1, 0.1)
    can2.cd(4)
    tgraph_t2_frac_lumi.Draw("AP")
    tgraph_t2_frac_lumi.GetYaxis().SetRangeUser(-0.1, 0.1)
    can2.cd(2)
    tgraph_t2_abs_sbil.Draw("AP")
    tgraph_t2_abs_sbil.GetYaxis().SetRangeUser(-0.1, 0.1)
    can2.cd(5)
    tgraph_t2_abs_lumi.Draw("AP")
    tgraph_t2_abs_lumi.GetYaxis().SetRangeUser(-0.1, 0.1)
    can2.cd(3)
    tgraph_t2_rms_sbil.Draw("AP")
    tgraph_t2_rms_sbil.GetYaxis().SetRangeUser(-0.1, 0.1)
    can2.cd(6)
    tgraph_t2_rms_lumi.Draw("AP")
    tgraph_t2_rms_lumi.GetYaxis().SetRangeUser(-0.1, 0.1)
    can2.Update()
    can2.SaveAs(args.outputdir+"/type2_raw"+args.label+".png")
    can2.SaveAs(args.outputdir+"/type2_raw"+args.label+".C")


# Overall Lumi change vs Average SBIL  -> Per fill / profiled
# 
def find_fit_range(graph):
    x_list=[]
    y_list=[]
    nps = graph.GetN()
    for ip in range(nps):
        if graph.GetX()[ip]>0.1:
            x_list.append(graph.GetX()[ip])
            y_list.append(graph.GetY()[ip])

    x_list.sort()
    y_list.sort()
    counts={}
    for ix in x_list:
        if math.isnan(ix): continue
        if not int(ix/1.0) in counts.keys():
            counts[int(ix/1.0)]=1
        else:
            counts[int(ix/1.0)]+=1

    min_key = 40000
    max_key = -1
    miny_key = -0.1
    maxy_key = 0.1
    if len(y_list)>1:
        miny_key = np.min(y_list)
        maxy_key = np.max(y_list)
    #print(counts)
    for key in counts.keys():
        if counts[key]>1 and key>0.1:
            if key<min_key:
                min_key=key

            if key>max_key:
                max_key=key
    #print(min_key, max_key, miny_key, maxy_key)
    return 1.0*min_key, 1.0*(max_key+1), miny_key, maxy_key

def remove_outlier(graph):
    mean = graph.GetMean(2)
    np = graph.GetN()
    ip=0
    while ip<np:
        if abs(graph.GetY()[ip]-mean)>mean*0.20 or graph.GetX()[ip]<=0.001:
            print('remove ', ip, graph.GetX()[ip], graph.GetY()[ip])
            graph.RemovePoint(ip)
            ip-=1
            np-=1
        ip+=1

fill2p0p1={}
makeperfill=True
if makeperfill:
    P0_vs_Fill = ROOT.TGraphErrors()
    P1_vs_Fill = ROOT.TGraphErrors()
    lumicorrect_fill = ROOT.TGraphErrors()
    P0_vs_nactive = ROOT.TGraphErrors()
    P1_vs_nactive = ROOT.TGraphErrors()
    lumicorrect_nactive = ROOT.TGraphErrors()
    P0_vs_lumi = ROOT.TGraphErrors()
    P1_vs_lumi = ROOT.TGraphErrors()
    lumicorrect_lumi = ROOT.TGraphErrors()

    ifill=0
    cantemp=ROOT.TCanvas("cantemp","cantemp",1000,700)
    cantemp.SetTickx()
    cantemp.SetTicky()
    for ikey in thist_t1_frac_sbil.keys():
        print(ikey)
        if thist_t1_frac_sbil[ikey][1]<3 or ikey=='total': 
            print('not enough points in ', ikey)
            continue
        cantemp.Update()
        grtmp=thist_t1_frac_sbil[ikey][0]
        low_x,high_x,lowy,highy = find_fit_range(grtmp)
        if grtmp.GetN()<2:
            print('after remove_outlier not enough points in ', ikey)
            continue
        grtmp.SetTitle("Fill %s;IL (/#muB);Type1 fraction"%str(ikey))
        grtmp.SetMarkerColor(ROOT.kBlue)
        grtmp.GetYaxis().SetTitleOffset(1.5)
        grtmp.SetMarkerStyle(23)
        grtmp.SetMarkerSize(1.2)
        grtmp.SetMarkerColor(ROOT.kBlue)
        grtmp.GetYaxis().SetRangeUser(lowy,highy)
        grtmp.Draw("APE0Z")
        FitResult = grtmp.Fit("pol1", "Q","",max(0.2,low_x), high_x)
        FitResult = grtmp.GetFunction("pol1")
        try:
            p0 = FitResult.GetParameter(0)
            p0err = FitResult.GetParError(0)
            p1 = FitResult.GetParameter(1)
            p1err = FitResult.GetParError(1)
            if p0err < 1 and p1err < 1:
                P0_vs_Fill.SetPoint(ifill, float(ikey), p0)
                P0_vs_Fill.SetPointError(ifill, 0, p0err)
                P1_vs_Fill.SetPoint(ifill, float(ikey), p1)
                P1_vs_Fill.SetPointError(ifill, 0, p1err)
                lumicorrect_fill.SetPoint(ifill, float(ikey), thist_lumi_correction[ikey][2]/thist_lumi_correction[ikey][3])
                P0_vs_nactive.SetPoint(ifill, thist_lumi_correction[ikey][4], p0)
                P0_vs_nactive.SetPointError(ifill, 0, p0err)
                P1_vs_nactive.SetPoint(ifill, thist_lumi_correction[ikey][4], p1)
                P1_vs_nactive.SetPointError(ifill, 0, p1err)
                lumicorrect_nactive.SetPoint(ifill, thist_lumi_correction[ikey][4], thist_lumi_correction[ikey][2]/thist_lumi_correction[ikey][3])
                lumitotal=0
                if str(int(ikey)) in fill2lumi.keys(): lumitotal=1e-6*float(fill2lumi[str(int(ikey))])
                else:
                    print(ikey," not in fill2lumi" )
                    exit()
                P0_vs_lumi.SetPoint(ifill, lumitotal, p0)
                P0_vs_lumi.SetPointError(ifill, 0, p0err)
                P1_vs_lumi.SetPoint(ifill, lumitotal, p1)
                P1_vs_lumi.SetPointError(ifill, 0, p1err)
                fill2p0p1[str(int(ikey))]=[p0,p0err,p1,p1err]
                lumicorrect_lumi.SetPoint(ifill, lumitotal, thist_lumi_correction[ikey][2]/thist_lumi_correction[ikey][3])
                ifill+=1
            else:
                print(ikey, " has large uncertianry")
            text.Draw("same")
            text2.Draw("same")
            f1 = ROOT.TF1("f1","[0]+[1]*x",max(0.2,low_x), high_x)
            f1.SetParameters(p0,p1)
            f1.Draw("same")
            fitfuc="Fit: p0+p1*x"
            fitfuc1="P0=%0.9f #pm%0.5f"%(p0,p0err)
            fitfuc2="P1=%0.9f #pm%0.5f"%(p1,p1err)
            #if p1<0: fitfuc="Fit: %0.5f%0.5f*x"%(p0,p1)
            text3=ROOT.TLatex(0.4,0.8,fitfuc)
            text3.SetNDC()
            text3.SetTextFont(62)
            text3.SetTextSize(0.05)
            text3.Draw("same")
            text4=ROOT.TLatex(0.4,0.7,fitfuc1)
            text4.SetNDC()
            text4.SetTextFont(62)
            text4.SetTextSize(0.05)
            text4.Draw("same")
            text5=ROOT.TLatex(0.4,0.6,fitfuc2)
            text5.SetNDC()
            text5.SetTextFont(62)
            text5.SetTextSize(0.05)
            text5.Draw("same")
            cantemp.SaveAs(args.outputdir+"/OverallLumiCor_Fill_"+str(ikey)+args.label+".png")
            cantemp.SaveAs(args.outputdir+"/OverallLumiCor_Fill_"+str(ikey)+args.label+".C")
        except:
            print("Problem with the Fill:", ikey)
        #input()
    fill2p0p1keys=fill2p0p1.keys()
    fill2p0p1keys.sort()
    for ikey in fill2p0p1keys:
        print( "'%s':[%0.8f,%0.8f,%0.8f,%0.8f],"%(ikey,fill2p0p1[ikey][0],fill2p0p1[ikey][1],fill2p0p1[ikey][2],fill2p0p1[ikey][3],) )
    cantemp.Update()
    P0_vs_Fill.GetXaxis().SetTitle("Fill")
    P0_vs_Fill.GetYaxis().SetTitle("p0")
    P0_vs_Fill.GetYaxis().SetTitleOffset(1.2)
    P0_vs_Fill.GetYaxis().SetRangeUser(0.98, 1.02)
    P0_vs_Fill.SetMarkerStyle(23)
    P0_vs_Fill.SetMarkerSize(0.7)
    P0_vs_Fill.SetMarkerColor(ROOT.kBlue)
    P0_vs_Fill.Draw("APE0Z")
    cantemp.SaveAs(args.outputdir+"/P0_vs_Fill_"+args.label+".png")
    cantemp.SaveAs(args.outputdir+"/P0_vs_Fill_"+args.label+".C")
    cantemp.Update()
    P1_vs_Fill.GetXaxis().SetTitle("Fill")
    P1_vs_Fill.GetYaxis().SetTitle("p1")
    P1_vs_Fill.GetYaxis().SetTitleOffset(1.2)
    P1_vs_Fill.GetYaxis().SetRangeUser(-0.002, 0.002)
    P1_vs_Fill.SetMarkerStyle(23)
    P1_vs_Fill.SetMarkerSize(0.7)
    P1_vs_Fill.SetMarkerColor(ROOT.kBlue)
    P1_vs_Fill.Draw("APE0Z")
    cantemp.SaveAs(args.outputdir+"/P1_vs_Fill_"+args.label+".png")
    cantemp.SaveAs(args.outputdir+"/P1_vs_Fill_"+args.label+".C")
    cantemp.Update()
    thist_lumi_correction_prof=thist_lumi_correction['total'].ProfileX("_pfx",1,-1,"S")
    thist_lumi_correction_prof.GetYaxis().SetRangeUser(0.95,1.05)
    thist_lumi_correction_prof.Draw()
    cantemp.SaveAs(args.outputdir+"/OverallLumiCor_profile"+args.label+".png")
    cantemp.SaveAs(args.outputdir+"/OverallLumiCor_profile"+args.label+".C")
    cantemp.Update()
    lumicorrect_fill.GetXaxis().SetTitle("Fill")
    lumicorrect_fill.GetYaxis().SetTitle("Correction factor")
    lumicorrect_fill.GetYaxis().SetTitleOffset(1.2)
    lumicorrect_fill.GetYaxis().SetRangeUser(0.98, 1.02)
    lumicorrect_fill.SetMarkerStyle(23)
    lumicorrect_fill.SetMarkerSize(0.7)
    lumicorrect_fill.SetMarkerColor(ROOT.kBlue)
    lumicorrect_fill.Draw("APE0Z")
    cantemp.SaveAs(args.outputdir+"/OverallLumiCor_fill_"+args.label+".png")
    cantemp.SaveAs(args.outputdir+"/OverallLumiCor_fill_"+args.label+".C")
    make_lumi=True
    if make_lumi:
        cantemp.Update()
        P0_vs_lumi.GetXaxis().SetTitle("Lumi (fb^{-1})")
        P0_vs_lumi.GetYaxis().SetTitle("p0")
        P0_vs_lumi.GetYaxis().SetTitleOffset(1.2)
        P0_vs_lumi.GetYaxis().SetRangeUser(0.98, 1.02)
        P0_vs_lumi.SetMarkerStyle(23)
        P0_vs_lumi.SetMarkerSize(0.7)
        P0_vs_lumi.SetMarkerColor(ROOT.kBlue)
        P0_vs_lumi.Draw("APE0Z")
        cantemp.SaveAs(args.outputdir+"/P0_vs_lumi_"+args.label+".png")
        cantemp.SaveAs(args.outputdir+"/P0_vs_lumi_"+args.label+".C")
        cantemp.Update()
        P1_vs_lumi.GetXaxis().SetTitle("Lumi (fb^{-1})")
        P1_vs_lumi.GetYaxis().SetTitle("p1")
        P1_vs_lumi.GetYaxis().SetTitleOffset(1.2)
        P1_vs_lumi.GetYaxis().SetRangeUser(-0.002, 0.002)
        P1_vs_lumi.SetMarkerStyle(23)
        P1_vs_lumi.SetMarkerSize(0.7)
        P1_vs_lumi.SetMarkerColor(ROOT.kBlue)
        P1_vs_lumi.Draw("APE0Z")
        cantemp.SaveAs(args.outputdir+"/P1_vs_lumi_"+args.label+".png")
        cantemp.SaveAs(args.outputdir+"/P1_vs_lumi_"+args.label+".C")
        cantemp.Update()
        lumicorrect_lumi.GetXaxis().SetTitle("Lumi (fb^{-1})")
        lumicorrect_lumi.GetYaxis().SetTitle("Correction factor")
        lumicorrect_lumi.GetYaxis().SetTitleOffset(1.2)
        lumicorrect_lumi.GetYaxis().SetRangeUser(0.98, 1.02)
        lumicorrect_lumi.SetMarkerStyle(23)
        lumicorrect_lumi.SetMarkerSize(0.7)
        lumicorrect_lumi.SetMarkerColor(ROOT.kBlue)
        lumicorrect_lumi.Draw("APE0Z")
        cantemp.SaveAs(args.outputdir+"/OverallLumiCor_lumi_"+args.label+".png")
        cantemp.SaveAs(args.outputdir+"/OverallLumiCor_lumi_"+args.label+".C")
    make_nactive=False
    if make_nactive:
        cantemp.Update()
        P0_vs_nactive.GetXaxis().SetTitle("NBX active")
        P0_vs_nactive.GetYaxis().SetTitle("p0")
        P0_vs_nactive.GetYaxis().SetTitleOffset(1.2)
        P0_vs_nactive.GetYaxis().SetRangeUser(0.9, 1.1)
        P0_vs_nactive.SetMarkerStyle(23)
        P0_vs_nactive.SetMarkerSize(0.7)
        P0_vs_nactive.SetMarkerColor(ROOT.kBlue)
        P0_vs_nactive.Draw("APE0Z")
        cantemp.SaveAs(args.outputdir+"/P0_vs_nactive_"+args.label+".png")
        cantemp.SaveAs(args.outputdir+"/P0_vs_nactive_"+args.label+".C")
        cantemp.Update()
        P1_vs_nactive.GetXaxis().SetTitle("NBX active")
        P1_vs_nactive.GetYaxis().SetTitle("p1")
        P1_vs_nactive.GetYaxis().SetTitleOffset(1.2)
        P1_vs_nactive.GetYaxis().SetRangeUser(-0.02, 0.02)
        P1_vs_nactive.SetMarkerStyle(23)
        P1_vs_nactive.SetMarkerSize(0.7)
        P1_vs_nactive.SetMarkerColor(ROOT.kBlue)
        P1_vs_nactive.Draw("APE0Z")
        cantemp.SaveAs(args.outputdir+"/P1_vs_nactive_"+args.label+".png")
        cantemp.SaveAs(args.outputdir+"/P1_vs_nactive_"+args.label+".C")
        cantemp.Update()
        lumicorrect_nactive.GetXaxis().SetTitle("NBX active")
        lumicorrect_nactive.GetYaxis().SetTitle("Correction factor")
        lumicorrect_nactive.GetYaxis().SetTitleOffset(1.2)
        lumicorrect_nactive.GetYaxis().SetRangeUser(0.9, 1.1)
        lumicorrect_nactive.SetMarkerStyle(23)
        lumicorrect_nactive.SetMarkerSize(0.7)
        lumicorrect_nactive.SetMarkerColor(ROOT.kBlue)
        lumicorrect_nactive.Draw("APE0Z")
        cantemp.SaveAs(args.outputdir+"/OverallLumiCor_nactive_"+args.label+".png")
        cantemp.SaveAs(args.outputdir+"/OverallLumiCor_nactive_"+args.label+".C")


