# program to generate LUT for HF lumi that translates
# floating-point codes from the QIE8 to a fixed point
# linear scale

import numpy as np
import matplotlib.pyplot as plt

def canvas_margin(c1, c1_up, c1_down):
  c1_up.SetTopMargin( 0.07 )
  c1_up.SetBottomMargin( 0.01 )
  c1_up.SetLeftMargin( 0.15 )
  c1_up.SetRightMargin( 0.03 )

  c1_down.SetTopMargin( 0.01 )
  c1_down.SetBottomMargin( 0.4 )
  c1_down.SetLeftMargin( 0.15 )
  c1_down.SetRightMargin( 0.03 )

  c1.SetTopMargin( 0.05 )
  c1.SetBottomMargin( 0.13 )
  c1.SetRightMargin( 0.05 )
  c1.SetLeftMargin( 0.16 )

NewLUT = np.zeros(256)
nline=0
hnewLUT={}
for line in open('outLUT_2022.txt','r').readlines() : 
    nline+=1
    #if nline!=10: continue
    vals = line.split()
    if len(vals) !=263: 
        print(len(vals))
        exit()
    name=''
    for i in range(7) :
        name+='%s '%(vals[i])
    hnewLUT[name]=[]
    for i in range(256) :
        NewLUT[i] += float(vals[i+7])
        hnewLUT[name].append(float(vals[i+7]))

OldLUT = np.zeros(256)
nline=0
holdLUT={}
for line in open('luts_hf2017_v0.txt','r').readlines() :
    nline+=1
    #if nline!=10: continue
    vals = line.split()
    if len(vals) !=263:
        print(len(vals))
        exit()
    name=''
    for i in range(7) :
        name+='%s '%(vals[i])
    holdLUT[name]=[]
    for i in range(256) :
        OldLUT[i] += float(vals[i+7])
        holdLUT[name].append(float(vals[i+7]))

print(NewLUT[10:20])
print(OldLUT[10:20])

import ROOT as r
#r.gROOT.SetBatch(True)

channellist=hnewLUT.keys()

for ic in range(600):
    hnew = r.TH1F("",";ADC;Charge (fC)",256,0,256)
    hratio = r.TH1F("",";ADC;Charge (fC)",256,0,256)
    hold = r.TH1F("","old LUT;ADC;Charge (fC)",256,0,256)
    if ic==0:
        for i in range(256):
            hnew.SetBinContent(i+1,NewLUT[i])
            hratio.SetBinContent(i+1,NewLUT[i])
            hold.SetBinContent(i+1,OldLUT[i])
    else:
        for i in range(256):
            hnew.SetBinContent(i+1,hnewLUT[channellist[ic+1]][i])
            hratio.SetBinContent(i+1,hnewLUT[channellist[ic+1]][i])
            hold.SetBinContent(i+1,holdLUT[channellist[ic+1]][i])
    Drawratio=True
    if Drawratio:
        c = r.TCanvas("c","c",800,800) ;
        r.gStyle.SetOptStat(0)
        c.SetFillColor(0)
        c.SetBorderMode(0)
        c.SetBorderSize(2)
        c.SetFrameBorderMode(0)
        #------------>Primitives in pad: toppad
    
        toppad = r.TPad('toppad','toppad',0,0.32 ,1.0,1.0)
        bottompad = r.TPad('bottompad','bottompad',0,0.0,1.0,0.30)
        canvas_margin(c,toppad,bottompad)
        toppad.SetFillStyle(4000)
        toppad.SetFrameFillStyle(1000)
        toppad.SetFrameFillColor(0)
        toppad.SetFillColor(0)
        toppad.SetBorderMode(0)
        toppad.SetBorderSize(2)
        #toppad.SetLogy()
        toppad.SetFrameBorderMode(0)
        toppad.SetFrameBorderMode(0)
        toppad.SetLeftMargin(0.15)
        bottompad.SetFillStyle(4000)
        bottompad.SetFrameFillStyle(1000)
        bottompad.SetFrameFillColor(0)
        bottompad.SetFillColor(0)
        bottompad.SetBorderMode(0)
        bottompad.SetBorderSize(2)
        bottompad.SetFrameBorderMode(0)
        bottompad.SetFrameBorderMode(0)
        toppad.Draw()
        bottompad.Draw()
    
        c.cd()
        c.Update()
        c.RedrawAxis()
        cframe = c.GetFrame()
        cframe.Draw()
        toppad.cd()
        frame_4fa51a0__1 = r.TH1D("frame_4fa51a0__1","",256,0,256)
        frame_4fa51a0__1.GetXaxis().SetTitle("ADC")
        frame_4fa51a0__1.GetXaxis().SetLabelFont(42);
        frame_4fa51a0__1.GetXaxis().SetLabelSize(0.05);
        frame_4fa51a0__1.GetXaxis().SetTitleSize(0.05);
        frame_4fa51a0__1.GetXaxis().SetTitleOffset(1);
        frame_4fa51a0__1.GetXaxis().SetTitleFont(42);
        frame_4fa51a0__1.GetYaxis().SetTitle("Charge (fC)")
        frame_4fa51a0__1.GetYaxis().SetLabelFont(42);
        frame_4fa51a0__1.GetYaxis().SetLabelSize(0.05);
        frame_4fa51a0__1.GetYaxis().SetTitleSize(0.05);
        frame_4fa51a0__1.GetYaxis().SetTitleFont(42);
        frame_4fa51a0__1.GetYaxis().SetRangeUser(0.8e-2,1e4)
        frame_4fa51a0__1.GetXaxis().SetLabelOffset(999)
        frame_4fa51a0__1.GetXaxis().SetLabelSize(0)
        frame_4fa51a0__1.Draw("AXISSAME");
    
        hnew.SetLineColor(r.kRed)
        hnew.Draw("HIST")
        hold.Draw("same HIST")
        frame_4fa51a0__1.Draw("AXISSAME");
    
        leg = r.TLegend(0.25,0.73,0.4,0.85);
        leg.SetBorderSize(0);
        leg.SetLineStyle(1);
        leg.SetLineWidth(1);
        entry=leg.AddEntry(hnew,"New LUT","l");
        entry.SetFillStyle(1001);
        entry.SetMarkerStyle(8);
        entry.SetMarkerSize(1.5);
        entry.SetLineStyle(1);
        entry.SetLineWidth(2);
        entry.SetTextFont(42);
        entry.SetTextSize(0.04);
        entry=leg.AddEntry(hold,"Old LUT","l");
        entry.SetFillStyle(1001);
        entry.SetMarkerStyle(8);
        entry.SetMarkerSize(1.5);
        entry.SetLineStyle(1);
        entry.SetLineWidth(2);
        entry.SetTextFont(42);
        entry.SetTextSize(0.04);
        leg.Draw()
    
        if ic==0:
            tex = r.TLatex(0.45,0.54,'All channels');
        else:
            tex2 = r.TLatex(0.65,0.62,'crate, slot, fiber, fiberChan, iEta, iPhi, depth')
            tex2.SetNDC();
            tex2.SetTextAlign(31);
            tex2.SetTextFont(42);
            tex2.SetTextSize(0.04);
            tex2.SetLineWidth(2);
            tex2.Draw();
            tex = r.TLatex(0.5,0.54,channellist[ic+1].replace(' ',', '));
        tex.SetNDC();
        tex.SetTextAlign(31);
        tex.SetTextFont(42);
        tex.SetTextSize(0.05);
        tex.SetLineWidth(2);
        tex.Draw();
    
        bottompad.cd()
        frame_4fa51a0__2 = r.TH1D("frame_4fa51a0__2","",256,0,256)
        frame_4fa51a0__2.GetXaxis().SetTitle("ADC");
        frame_4fa51a0__2.GetXaxis().SetLabelFont(42);
        frame_4fa51a0__2.GetXaxis().SetLabelSize(0.1);
        frame_4fa51a0__2.GetXaxis().SetTitleSize(0.1);
        frame_4fa51a0__2.GetXaxis().SetTitleOffset(1);
        frame_4fa51a0__2.GetXaxis().SetTitleFont(42);
        frame_4fa51a0__2.GetYaxis().SetTitle("New/Old")
        frame_4fa51a0__2.GetYaxis().CenterTitle()
        frame_4fa51a0__2.GetYaxis().SetLabelFont(42);
        frame_4fa51a0__2.GetYaxis().SetLabelSize(0.1);
        frame_4fa51a0__2.GetYaxis().SetTitleSize(0.1);
        frame_4fa51a0__2.GetYaxis().SetTitleFont(42);
        frame_4fa51a0__2.GetYaxis().SetRangeUser(0,2)
        frame_4fa51a0__2.GetYaxis().SetNdivisions(5, 2, 0, r.kTRUE)
        frame_4fa51a0__2.GetYaxis().SetTitleOffset(0.4)
        frame_4fa51a0__2.Draw("AXISSAME")
        hratio.Divide(hold)
        hratio.Draw("same")
        line = r.TLine(0,1,256,1);
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
        frame_4fa51a0__2.Draw("AXISSAME")
    if ic==0:
        c.SaveAs('LUT_allchan.png')
    else:
        c.SaveAs('LUT_iEta_'+str(channellist[ic+1].split(' ')[4])+'_iPhi_'+str(channellist[ic+1].split(' ')[5])+'_depth_'+str(channellist[ic+1].split(' ')[6])+'.png')
    input()
    del hnew
    del hold
    del hratio
exit()
    
    

