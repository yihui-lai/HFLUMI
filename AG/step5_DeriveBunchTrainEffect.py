import ROOT
import sys
import os
import argparse
import json
from math import sqrt
import array
import math
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument("--hfoc", default="", help="The file for the hfoc luminosity")
parser.add_argument("--hfet", default="", help="The file for the hfet luminosity")
parser.add_argument("--outputdir", default="./", help="outputdir ")

args = parser.parse_args()

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)

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
    #print('x_list:',x_list)
    #print('y_list:',y_list)
    counts={}
    for ix in x_list:
        if math.isnan(ix): continue
        if not int(ix/1.0) in counts.keys():
            counts[int(ix/1.0)]=1
        else:
            counts[int(ix/1.0)]+=1
    min_key = 40000
    max_key = -1
    if len(y_list)>0:
        miny_key = np.min(y_list)
        maxy_key = np.max(y_list)
    #print(counts)
    for key in counts.keys():
        if counts[key]>1:
            if key<min_key:
                min_key=key

            if key>max_key:
                max_key=key
    return 1.0*min_key, 1.0*(max_key+1), miny_key, maxy_key



ocfile = ROOT.TChain("newtree")
ocfile.Add(args.hfoc)
etfile = ROOT.TChain("newtree")
etfile.Add(args.hfet)

doubleratio_list={}
for iev in ocfile:
   LB = '%0.0f_%0.0f_%0.0f'%(iev.FILL, iev.RUN, iev.LS) 
   if LB not in doubleratio_list.keys(): doubleratio_list[LB]={}
   if iev.sum_train!=0: 
       doubleratio_list[LB]['oc'] = iev.sum_lead/iev.sum_train
   else: doubleratio_list[LB]['oc'] = 0 
   doubleratio_list[LB]['sbil'] = iev.AveSBIL_leadBX
   doubleratio_list[LB]['oc_raw'] = iev.overallLumi_raw*1e-9
   doubleratio_list[LB]['oc_rawnoise'] = iev.overallLumi_rawSubnoise*1e-9
   doubleratio_list[LB]['oc_rawnoiselocalt1'] = iev.overallLumi_rawSublocalT1*1e-9

for iev in etfile:
   LB = '%0.0f_%0.0f_%0.0f'%(iev.FILL, iev.RUN, iev.LS)
   if LB not in doubleratio_list.keys(): doubleratio_list[LB]={}
   if iev.sum_train!=0: doubleratio_list[LB]['et'] = iev.sum_lead/iev.sum_train
   else: doubleratio_list[LB]['et'] = 0
   doubleratio_list[LB]['et_raw'] = iev.overallLumi_raw*1e-9
   doubleratio_list[LB]['et_rawnoise'] = iev.overallLumi_rawSubnoise*1e-9
   doubleratio_list[LB]['et_rawnoiselocalt1'] = iev.overallLumi_rawSublocalT1*1e-9

del ocfile
del etfile

fill2lumi={}
with open('supertable_22.txt') as f:
    lines = f.readlines()
    for iline in lines:
        #print(iline.split("\t"))
        fill2lumi[iline.split("\t")[0]]=(iline.split("\t")[1])

newfile=ROOT.TFile(args.outputdir+"/oc_bunchtrain.root", "RECREATE")
newtree=ROOT.TTree("newtree", "newtree")

FILL = array.array('I', [0])
RUN = array.array('I', [0])
LS = array.array('I', [0])
AveSBIL = array.array('d', [0])
doubleratio = array.array('d', [0])
newtree.Branch("FILL", FILL, "FILL/I")
newtree.Branch("RUN", RUN, "RUN/I")
newtree.Branch("LS", LS, "LS/I")
newtree.Branch("AveSBIL", AveSBIL, "AveSBIL/D")
newtree.Branch("doubleratio",doubleratio , "doubleratio/D")
    
doubleratio_vs_Fill = ROOT.TGraphErrors()
doubleratio_vs_sbil = ROOT.TGraphErrors()
hist_db_sbil = ROOT.TH2F("hist_db_sbil",";Average leading SBIL (Hz/#muB);double ratio",50,0.5,8,200,0.9,1.1)
hist_db_sbil_1 = ROOT.TH2F("hist_db_sbil_1",";Average leading SBIL (Hz/#muB);double ratio",50,0.5,8,200,0.9,1.1)
hist_db_sbil_2 = ROOT.TH2F("hist_db_sbil_2",";Average leading SBIL (Hz/#muB);double ratio",50,0.5,8,200,0.9,1.1)
hfochfet_vs_lumi_raw = ROOT.TGraphErrors()
hfochfet_vs_lumi_rawsubnoise = ROOT.TGraphErrors()
hfochfet_vs_lumi_rawsubnoisesublocalt1 = ROOT.TGraphErrors()

histfochfet_vs_lumi_raw = ROOT.TH2F("histfochfet_vs_lumi_raw",";Lumi (fb^{-1});hfoc/hfet",40,0,40,200,0.96,1.04)
histfochfet_vs_lumi_rawsubnoise = ROOT.TH2F("histfochfet_vs_lumi_rawsubnoise",";Lumi (fb^{-1});hfoc/hfet",40,0,40,200,0.96,1.04)
histfochfet_vs_lumi_rawsubnoisesublocalt1 = ROOT.TH2F("histfochfet_vs_lumi_rawsubnoisesublocalt1",";Lumi (fb^{-1});hfoc/hfet",40,0,40,200,0.96,1.04)

ifill=0
totallumi=0
doublerratio_perfill={}
LBs=doubleratio_list.keys()
LBs.sort()
print(LBs)
for LB in LBs:
        #print("LumiBlock: ", LB)
        if 'et' in doubleratio_list[LB].keys() and 'oc' in doubleratio_list[LB].keys():
            FILL[0] = int(LB.split("_")[0])
            RUN[0] = int(LB.split("_")[1])
            LS[0] = int(LB.split("_")[2])
            AveSBIL[0] =  doubleratio_list[LB]['sbil']
            if doubleratio_list[LB]['et']!=0: 
                #print(LB,  doubleratio_list[LB]['oc'],doubleratio_list[LB]['et'],doubleratio_list[LB]['oc']/doubleratio_list[LB]['et'])
                doubleratio[0] = doubleratio_list[LB]['oc']/doubleratio_list[LB]['et']
            else: doubleratio[0] = -1
            doubleratio_vs_Fill.SetPoint(ifill, FILL[0], doubleratio[0])
            doubleratio_vs_sbil.SetPoint(ifill, AveSBIL[0], doubleratio[0])
            hist_db_sbil.Fill(AveSBIL[0], doubleratio[0])
            if FILL[0]<8115: hist_db_sbil_1.Fill(AveSBIL[0], doubleratio[0])
            else: hist_db_sbil_2.Fill(AveSBIL[0], doubleratio[0])
            #if AveSBIL[0]>3 and doubleratio[0]<0.985:
                #print(LB,  doubleratio_list[LB]['oc'],doubleratio_list[LB]['et'],doubleratio_list[LB]['oc']/doubleratio_list[LB]['et'])
                #print(AveSBIL[0], doubleratio[0])
                #input()
            totallumi+=doubleratio_list[LB]['et_raw']
            #print(totallumi)
            if FILL[0] not in doublerratio_perfill.keys():
                doublerratio_perfill[FILL[0]] = [ROOT.TGraphErrors(),0]
                if doubleratio[0]>=0:
                    doublerratio_perfill[FILL[0]][0].SetPoint(doublerratio_perfill[FILL[0]][1],AveSBIL[0],doubleratio[0])
                    doublerratio_perfill[FILL[0]][1]+=1
            else:
                if doubleratio[0]>=0:
                    doublerratio_perfill[FILL[0]][0].SetPoint(doublerratio_perfill[FILL[0]][1],AveSBIL[0],doubleratio[0])
                    doublerratio_perfill[FILL[0]][1]+=1
            #print(LB.split("_")[0],fill2lumi[LB.split("_")[0]])
            if LB.split("_")[0] not in fill2lumi.keys(): 
                print(LB.split("_")[0],' not in fill2lumi')
                exit()
            else:
                hfochfet_vs_lumi_raw.SetPoint(ifill, 1e-6*float(fill2lumi[LB.split("_")[0]]), doubleratio_list[LB]['oc_raw']/doubleratio_list[LB]['et_raw'] )
                hfochfet_vs_lumi_rawsubnoise.SetPoint(ifill, 1e-6*float(fill2lumi[LB.split("_")[0]]), doubleratio_list[LB]['oc_rawnoise']/doubleratio_list[LB]['et_rawnoise'] )
                hfochfet_vs_lumi_rawsubnoisesublocalt1.SetPoint(ifill, 1e-6*float(fill2lumi[LB.split("_")[0]]), doubleratio_list[LB]['oc_rawnoiselocalt1']/doubleratio_list[LB]['et_rawnoiselocalt1'] )
            histfochfet_vs_lumi_raw.Fill(totallumi,doubleratio_list[LB]['oc_raw']/doubleratio_list[LB]['et_raw'])
            histfochfet_vs_lumi_rawsubnoise.Fill(totallumi,doubleratio_list[LB]['oc_rawnoise']/doubleratio_list[LB]['et_rawnoise'])
            histfochfet_vs_lumi_rawsubnoisesublocalt1.Fill(totallumi,doubleratio_list[LB]['oc_rawnoiselocalt1']/doubleratio_list[LB]['et_rawnoiselocalt1'])
            ifill+=1
            newtree.Fill()
newfile.WriteTObject(newtree, "newtree")

cantemp=ROOT.TCanvas("cantemp","cantemp",1000,700)
cantemp.SetTickx()
cantemp.SetTicky()

P0_vs_Fill = ROOT.TGraphErrors()
P1_vs_Fill = ROOT.TGraphErrors()
lumicorrect_lumi = ROOT.TGraphErrors()
lumicorrect_fill = ROOT.TGraphErrors()

text=ROOT.TLatex(0.72,0.92,"2022  (13.6 TeV)")
text.SetNDC()
text.SetTextFont(62)
text.SetTextSize(0.05)
text2=ROOT.TLatex(0.15,0.92,"CMS #bf{#scale[0.75]{#it{Preliminary}}}")
text2.SetNDC()
text2.SetTextSize(0.05)
text2.SetTextFont(62)
ifill=0
Fills=doublerratio_perfill.keys()
Fills.sort()
print(Fills)
for ikey in Fills:
        cantemp.Update()
        if doublerratio_perfill[ikey][0].GetN()<=0: continue
        #print('Fill ', ikey, doublerratio_perfill[ikey][0].GetN())
        low_x,high_x,lowy,highy=find_fit_range(doublerratio_perfill[ikey][0])
        doublerratio_perfill[ikey][0].SetTitle("Fill %s;Average SBIL (Hz/#muB);Correction Factor"%str(ikey))
        doublerratio_perfill[ikey][0].SetMarkerColor(ROOT.kBlue)
        doublerratio_perfill[ikey][0].GetYaxis().SetTitleOffset(1.5)
        doublerratio_perfill[ikey][0].SetMarkerStyle(23)
        doublerratio_perfill[ikey][0].SetMarkerSize(1.2)
        doublerratio_perfill[ikey][0].SetMarkerColor(ROOT.kBlue)
        doublerratio_perfill[ikey][0].GetYaxis().SetRangeUser(lowy,highy)
        doublerratio_perfill[ikey][0].GetXaxis().SetRangeUser(low_x,high_x)
        doublerratio_perfill[ikey][0].Draw("APE0Z")
        #FitResult = doublerratio_perfill[ikey][0].Fit("pol1", "Q","",max(0.2,low_x), high_x)
        #FitResult = doublerratio_perfill[ikey][0].GetFunction("pol1")
        #p0 = FitResult.GetParameter(0)
        #p0err = FitResult.GetParError(0)
        #p1 = FitResult.GetParameter(1)
        #p1err = FitResult.GetParError(1)
        ##print(p0, p0err, p1, p1err)
        #if p0err < 0.01 and p1err < 0.001:
        #    P0_vs_Fill.SetPoint(ifill, float(ikey), p0)
        #    P0_vs_Fill.SetPointError(ifill, 0, p0err)
        #    P1_vs_Fill.SetPoint(ifill, float(ikey), p1)
        #    P1_vs_Fill.SetPointError(ifill, 0, p1err)
        #    #lumicorrect_fill.SetPoint(ifill, float(ikey), doublerratio_perfill[ikey][2]/doublerratio_perfill[ikey][3])
        #    #lumitotal=0
        #    #if str(int(ikey)) in fill2lumi.keys(): lumitotal=1e-6*float(fill2lumi[str(int(ikey))])
        #    #else:
        #    #    print(ikey," not in fill2lumi" )
        #    #    exit()
        #    #lumicorrect_lumi.SetPoint(ifill, lumitotal, doublerratio_perfill[ikey][2]/doublerratio_perfill[ikey][3])
        #    #print(ifill, lumitotal,)
        #    ifill+=1
        #else:
        #    print(ikey, " has large uncertianry")
        text.Draw("same")
        text2.Draw("same")
        #fitfuc="Fit: p0+p1*x"
        #fitfuc1="P0=%0.5f#pm%0.5f"%(p0,p0err)
        #fitfuc2="P1=%0.5f#pm%0.5f"%(p1,p1err)
        #text3=ROOT.TLatex(0.4,0.8,fitfuc)
        #text3.SetNDC()
        #text3.SetTextFont(62)
        #text3.SetTextSize(0.05)
        #text3.Draw("same")
        #text4=ROOT.TLatex(0.4,0.7,fitfuc1)
        #text4.SetNDC()
        #text4.SetTextFont(62)
        #text4.SetTextSize(0.05)
        #text4.Draw("same")
        #text5=ROOT.TLatex(0.4,0.6,fitfuc2)
        #text5.SetNDC()
        #text5.SetTextFont(62)
        #text5.SetTextSize(0.05)
        #text5.Draw("same")
        cantemp.SaveAs(args.outputdir+"/doubleratio_Fill_"+str(ikey)+".png")
        cantemp.SaveAs(args.outputdir+"/doubleratio_Fill_"+str(ikey)+".C")


cantemp.Update()
doubleratio_vs_Fill.GetXaxis().SetTitle("Fill")
doubleratio_vs_Fill.GetYaxis().SetTitle("double ratio")
doubleratio_vs_Fill.GetYaxis().SetTitleOffset(1.2)
doubleratio_vs_Fill.GetYaxis().SetRangeUser(0.9, 1.1)
doubleratio_vs_Fill.SetMarkerStyle(23)
doubleratio_vs_Fill.SetMarkerSize(0.7)
doubleratio_vs_Fill.SetMarkerColor(ROOT.kBlue)
doubleratio_vs_Fill.Draw("APE0Z")
cantemp.SaveAs(args.outputdir+"/doubleratio_vs_Fill.png")
cantemp.SaveAs(args.outputdir+"/doubleratio_vs_Fill.C")
cantemp.Update()
doubleratio_vs_sbil.GetXaxis().SetTitle("Average leading SBIL (Hz/#muB)")
doubleratio_vs_sbil.GetYaxis().SetTitle("double ratio")
doubleratio_vs_sbil.GetYaxis().SetTitleOffset(1.2)
doubleratio_vs_sbil.GetYaxis().SetRangeUser(0.9, 1.1)
doubleratio_vs_sbil.SetMarkerStyle(23)
doubleratio_vs_sbil.SetMarkerSize(0.7)
doubleratio_vs_sbil.SetMarkerColor(ROOT.kBlue)
doubleratio_vs_sbil.Draw("APE0Z")
cantemp.SaveAs(args.outputdir+"/doubleratio_vs_sbil.png")
cantemp.SaveAs(args.outputdir+"/doubleratio_vs_sbil.C")
cantemp.Update()
thist_prof=hist_db_sbil.ProfileX()
thist_prof.GetXaxis().SetTitle("Average leading SBIL (Hz/#muB)")
thist_prof.GetYaxis().SetTitle("double ratio")
thist_prof.GetYaxis().SetTitleOffset(1.2)
thist_prof.GetYaxis().SetRangeUser(0.986, 1.002)
thist_prof.SetMarkerStyle(23)
thist_prof.SetMarkerSize(0.7)
thist_prof.SetLineColor(ROOT.kBlue)
thist_prof.SetMarkerColor(ROOT.kBlue)
thist_prof.Draw("E")
dofit=False
if dofit:
        FitResult = thist_prof.Fit("pol1", "MF", "")
        FitResult = thist_prof.GetFunction("pol1")
        try:
            p0 = FitResult.GetParameter(0)
            p0err = FitResult.GetParError(0)
            p1 = FitResult.GetParameter(1)
            p1err = FitResult.GetParError(1)
        except:
            print("Problem with the Fill:", ikey)

cantemp.SaveAs(args.outputdir+"/doubleratio_vs_sbil_profile.png")
cantemp.SaveAs(args.outputdir+"/doubleratio_vs_sbil_profile.C")

cantemp.Update()
thist_prof=hist_db_sbil_1.ProfileX()
thist_prof.GetXaxis().SetTitle("Average leading SBIL (Hz/#muB)")
thist_prof.GetYaxis().SetTitle("double ratio")
thist_prof.GetYaxis().SetTitleOffset(1.2)
thist_prof.GetYaxis().SetRangeUser(0.986, 1.002)
thist_prof.SetMarkerStyle(23)
thist_prof.SetMarkerSize(0.7)
thist_prof.SetLineColor(ROOT.kBlue)
thist_prof.SetMarkerColor(ROOT.kBlue)
thist_prof.Draw("E")
if dofit:
        FitResult = thist_prof.Fit("pol1", "MF", "")
        FitResult = thist_prof.GetFunction("pol1")
        try:
            p0 = FitResult.GetParameter(0)
            p0err = FitResult.GetParError(0)
            p1 = FitResult.GetParameter(1)
            p1err = FitResult.GetParError(1)
        except:
            print("Problem with the Fill:", ikey)
cantemp.SaveAs(args.outputdir+"/doubleratio_vs_sbil_profile_before8115.png")
cantemp.SaveAs(args.outputdir+"/doubleratio_vs_sbil_profile_before8115.C")


cantemp.Update()
thist_prof2=hist_db_sbil_2.ProfileX()
thist_prof2.GetXaxis().SetTitle("Average leading SBIL (Hz/#muB)")
thist_prof2.GetYaxis().SetTitle("double ratio")
thist_prof2.GetYaxis().SetTitleOffset(1.2)
thist_prof2.GetYaxis().SetRangeUser(0.986, 1.002)
thist_prof2.SetMarkerStyle(23)
thist_prof2.SetMarkerSize(0.7)
thist_prof2.SetLineColor(ROOT.kRed)
thist_prof2.SetMarkerColor(ROOT.kRed)
thist_prof2.Draw("E")
if dofit:
        FitResult = thist_prof2.Fit("pol1", "MF", "")
        FitResult = thist_prof2.GetFunction("pol1")
        try:
            p0 = FitResult.GetParameter(0)
            p0err = FitResult.GetParError(0)
            p1 = FitResult.GetParameter(1)
            p1err = FitResult.GetParError(1)
        except:
            print("Problem with the Fill:", ikey)
cantemp.SaveAs(args.outputdir+"/doubleratio_vs_sbil_profile_after8115.png")
cantemp.SaveAs(args.outputdir+"/doubleratio_vs_sbil_profile_after8115.C")

cantemp.Update()
thist_prof.Draw("E")
thist_prof2.Draw("E same")
cantemp.SaveAs(args.outputdir+"/doubleratio_vs_sbil_profile_overlay.png")
cantemp.SaveAs(args.outputdir+"/doubleratio_vs_sbil_profile_overlay.C")

cantemp.Update()
thist_prof=histfochfet_vs_lumi_raw.ProfileX()
thist_prof2=histfochfet_vs_lumi_rawsubnoise.ProfileX()
thist_prof3=histfochfet_vs_lumi_rawsubnoisesublocalt1.ProfileX()
thist_prof.GetXaxis().SetTitle("Lumi (fb^{-1})")
thist_prof.GetYaxis().SetTitle("hfoc/hfet")
thist_prof.GetYaxis().SetTitleOffset(1.2)
thist_prof.GetYaxis().SetRangeUser(0.96, 1.04)
thist_prof.SetMarkerStyle(23)
thist_prof.SetMarkerSize(0.7)
thist_prof.SetLineColor(ROOT.kRed)
thist_prof.SetMarkerColor(ROOT.kRed)
thist_prof.Draw("E")

thist_prof2.SetMarkerStyle(23)
thist_prof2.SetMarkerSize(0.7)
thist_prof2.SetLineColor(ROOT.kBlue)
thist_prof2.SetMarkerColor(ROOT.kBlue)
thist_prof2.Draw("E same")

thist_prof3.SetMarkerStyle(23)
thist_prof3.SetMarkerSize(0.7)
thist_prof3.SetLineColor(ROOT.kGreen)
thist_prof3.SetMarkerColor(ROOT.kGreen)
thist_prof3.Draw("E same")

cantemp.SaveAs(args.outputdir+"/doubleratio_hfoc_hfet_ratio.png")
cantemp.SaveAs(args.outputdir+"/doubleratio_hfoc_hfet_ratio.C")

#cantemp.Update()
#P0_vs_Fill.GetXaxis().SetTitle("Fill")
#P0_vs_Fill.GetYaxis().SetTitle("p0")
#P0_vs_Fill.GetYaxis().SetTitleOffset(1.2)
#P0_vs_Fill.GetYaxis().SetRangeUser(0.9, 1.1)
#P0_vs_Fill.SetMarkerStyle(23)
#P0_vs_Fill.SetMarkerSize(0.7)
#P0_vs_Fill.SetMarkerColor(ROOT.kBlue)
#P0_vs_Fill.Draw("APE0Z")
#cantemp.SaveAs(args.outputdir+"/doubleratio_P0_vs_Fill_.png")
#cantemp.SaveAs(args.outputdir+"/doubleratio_P0_vs_Fill_.C")
#cantemp.Update()
#P1_vs_Fill.GetXaxis().SetTitle("Fill")
#P1_vs_Fill.GetYaxis().SetTitle("p1")
#P1_vs_Fill.GetYaxis().SetTitleOffset(1.2)
#P1_vs_Fill.GetYaxis().SetRangeUser(-0.02, 0.02)
#P1_vs_Fill.SetMarkerStyle(23)
#P1_vs_Fill.SetMarkerSize(0.7)
#P1_vs_Fill.SetMarkerColor(ROOT.kBlue)
#P1_vs_Fill.Draw("APE0Z")
#cantemp.SaveAs(args.outputdir+"/doubleratio_P1_vs_Fill_.png")
#cantemp.SaveAs(args.outputdir+"/doubleratio_P1_vs_Fill_.C")
#cantemp.Update()
#lumicorrect_fill.GetXaxis().SetTitle("Fill")
#lumicorrect_fill.GetYaxis().SetTitle("Correction factor")
#lumicorrect_fill.GetYaxis().SetTitleOffset(1.2)
#lumicorrect_fill.GetYaxis().SetRangeUser(0.9, 1.1)
#lumicorrect_fill.SetMarkerStyle(23)
#lumicorrect_fill.SetMarkerSize(0.7)
#lumicorrect_fill.SetMarkerColor(ROOT.kBlue)
#lumicorrect_fill.Draw("APE0Z")
#cantemp.SaveAs(args.outputdir+"/doubleratio_OverallLumiCor_fill_.png")
#cantemp.SaveAs(args.outputdir+"/doubleratio_OverallLumiCor_fill_.C")
#cantemp.Update()
#lumicorrect_lumi.GetXaxis().SetTitle("Lumi (fb^{-1})")
#lumicorrect_lumi.GetYaxis().SetTitle("Correction factor")
#lumicorrect_lumi.GetYaxis().SetTitleOffset(1.2)
#lumicorrect_lumi.GetYaxis().SetRangeUser(0.9, 1.1)
#lumicorrect_lumi.SetMarkerStyle(23)
#lumicorrect_lumi.SetMarkerSize(0.7)
#lumicorrect_lumi.SetMarkerColor(ROOT.kBlue)
#lumicorrect_lumi.Draw("APE0Z")
#cantemp.SaveAs(args.outputdir+"/doubleratio_OverallLumiCor_lumi_.png")
#cantemp.SaveAs(args.outputdir+"/doubleratio_OverallLumiCor_lumi_.C")


cantemp2=ROOT.TCanvas("cantemp","cantemp",1000,700)
cantemp2.SetTickx()
cantemp2.SetTicky()
tmultig=ROOT.TMultiGraph()
hfochfet_vs_lumi_raw.SetLineColor(ROOT.kRed)
hfochfet_vs_lumi_raw.SetMarkerColor(ROOT.kRed)
hfochfet_vs_lumi_rawsubnoise.SetLineColor(ROOT.kBlue)
hfochfet_vs_lumi_rawsubnoise.SetMarkerColor(ROOT.kBlue)
hfochfet_vs_lumi_rawsubnoisesublocalt1.SetLineColor(ROOT.kGreen)
hfochfet_vs_lumi_rawsubnoisesublocalt1.SetMarkerColor(ROOT.kGreen)
hfochfet_vs_lumi_raw.SetMarkerStyle(23)
hfochfet_vs_lumi_raw.SetMarkerSize(0.7)
hfochfet_vs_lumi_rawsubnoise.SetMarkerStyle(23)
hfochfet_vs_lumi_rawsubnoise.SetMarkerSize(0.7)
hfochfet_vs_lumi_rawsubnoisesublocalt1.SetMarkerStyle(23)
hfochfet_vs_lumi_rawsubnoisesublocalt1.SetMarkerSize(0.7)
tmultig.Add(hfochfet_vs_lumi_raw)
#tmultig.Add(hfochfet_vs_lumi_rawsubnoise)
tmultig.Add(hfochfet_vs_lumi_rawsubnoisesublocalt1)
tmultig.Draw("ap")
tmultig.GetXaxis().SetTitle("Lumi (fb^{-1})")
tmultig.GetYaxis().SetTitle("hfoc/hfet")
tmultig.GetYaxis().SetTitleOffset(1.2)
tmultig.GetYaxis().SetRangeUser(0.5, 1.5)
legend=ROOT.TLegend(0.5,0.7,0.9,0.9)
legend.AddEntry(hfochfet_vs_lumi_raw,"raw","l")
legend.AddEntry(hfochfet_vs_lumi_rawsubnoisesublocalt1,"raw + local T1 correction","l")
legend.Draw()
cantemp2.SaveAs(args.outputdir+"/doubleratio_hfoc_hfet_ratiov2.png")
cantemp2.SaveAs(args.outputdir+"/doubleratio_hfoc_hfet_ratiov2.C")



newfile.Close()


