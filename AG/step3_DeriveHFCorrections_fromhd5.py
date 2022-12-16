import sys, os
from math import exp
import argparse
import subprocess
import ROOT
import array
import numpy as np
import pathlib
import pandas as pd
import tables

parser=argparse.ArgumentParser()
parser.add_argument("-f", "--file", default="", help="The path to a hd5 file")
parser.add_argument("-m", "--method", default="hfet", help="hfet")
parser.add_argument("--fills", default="", help="Comma separated list of fills.")
parser.add_argument("--auto", default=True, action="store_true", help="Determine the runs from the certtree")
parser.add_argument("-l", "--label", default="", help="The label for outputs")
parser.add_argument('-b', '--batch',   action='store_true', default=False, help="Batch mode (doesn't make GUI TCanvases)")
parser.add_argument('--nLSInLumiBlock', default=50, type=int, help="Number of LSs to group")
parser.add_argument('-a', '--all', action='store_true', default=False, help="Store all the data into a single hist")
parser.add_argument('-p', '--plot', action='store_true', default=False, help="Make plots")
parser.add_argument('-u', '--undoAG', action='store_true', default=False, help="undo after glow")

verbose=1
BXLength=3564
zeroes = array.array('d', [0.]*BXLength)
args=parser.parse_args()
stable_fills = [7920, 7921, 7923, 7960, 7963, 7965, 7966, 7967, 7969, 7978, 8007, 8016, 8017, 8018, 8019, 8020, 8022, 8023, 8027, 8030, 8033, 8043, 8057, 8058, 8059, 8062, 8063, 8067, 8068, 8072, 8073, 8076, 8078, 8079, 8081, 8083, 8087, 8088, 8091, 8094, 8098, 8099, 8100, 8102, 8103, 8106, 8108, 8111, 8112, 8113, 8115, 8118, 8120, 8121, 8124, 8125, 8128, 8132, 8136, 8142, 8143, 8144, 8146, 8147, 8148, 8149, 8151, 8177, 8178, 8181, 8184, 8210, 8211, 8212, 8214, 8216, 8220, 8221, 8222, 8223, 8225, 8226, 8228, 8230, 8232, 8233, 8236, 8238, 8245, 8247, 8248, 8253, 8260, 8263, 8267, 8269, 8272, 8273, 8274, 8276, 8289, 8293, 8295, 8297, 8299, 8301, 8302, 8304, 8305, 8306, 8307, 8309, 8311, 8312, 8313, 8314, 8315, 8316, 8317, 8319, 8320, 8321, 8322, 8324, 8327, 8330, 8331, 8333, 8334, 8335, 8342, 8345, 8347, 8379, 8381, 8383, 8385, 8387, 8389, 8392, 8395, 8398, 8399, 8401, 8402, 8412, 8413, 8456, 8471, 8474, 8479, 8484, 8489, 8491, 8496]


if args.method=='hfet':
    sigvis = 3364.0
elif args.method=='hfoc':
    sigvis = 962.8
else:
    print('unknown method')
    exit()
scale= 11245.6/sigvis

# This script reads hd5 files and calculate residuals, derive local type1 corrections ... 

def hdf_to_pandas(file_path, node):
    # Read a single hd5 file
    try:
        return pd.read_hdf(file_path, key=node, mode='r') #[['fillnum','runnum','lsnum','nbnum', 'avg', 'bx']]
    except KeyError as k:
        return pd.DataFrame()
    except ValueError as v:
        with tables.open_file(file_path, 'r') as table:
            data = pd.DataFrame(data=[table.root[node][:][col] for col in table.root[node].colnames]).T
            data.columns = table.root[node].colnames
        return data

def read_fill(path,fill,node):
    # Read multiple files and concat them
    return pd.concat([hdf_to_pandas(f, node) for f in pathlib.Path(path+'/'+fill+'/').glob('*.hd5')])

# group Lumi Block, create name list
# calculate redisual, type1/2
# np.roll(x2, -1, axis=1)
def makeplot(allLumiPerBX,allCorrLumiPerBX, rawdata_mask):
    print("make plot")
    can=ROOT.TCanvas("can","",1200,1200)
    can.Divide(1,2)
    can.cd(1)
    allLumiPerBX.GetXaxis().SetName("BXID")
    allCorrLumiPerBX.GetXaxis().SetName("BXID")
    #allLumiPerBX.SetMaximum(10)
    #allLumiPerBX.SetMinimum(-1)
    allLumiPerBX.SetTitle("Raw Data")
    #allLumiPerBX.SetStats(00000)
    allLumiPerBX.Draw("hist")
    rawdata_mask.SetLineColor(ROOT.kGreen)
    rawdata_mask.Draw("hist same")
    can.cd(2)
    allCorrLumiPerBX.SetTitle("+Noise/AfterGlow corrections")
    #allCorrLumiPerBX.SetStats(0000000)
    min_=allCorrLumiPerBX.GetBinContent(1)
    max_=allCorrLumiPerBX.GetBinContent(1)
    for i in range(3500):
        if allCorrLumiPerBX.GetBinContent(i) < min_: min_=allCorrLumiPerBX.GetBinContent(i)
        if allCorrLumiPerBX.GetBinContent(i) > max_: max_=allCorrLumiPerBX.GetBinContent(i)
    allCorrLumiPerBX.SetMaximum(max_+ (max_-min_))
    allCorrLumiPerBX.SetMinimum(min_ - (max_-min_))
    allCorrLumiPerBX.Draw("hist")
    rawdata_mask.Draw("hist same")
    input()


def makeplots(allLumiPerBX,allCorrLumiPerBX_raw, allLumiPerBX_newSBR):
    print("make plot")
    can=ROOT.TCanvas("can","",1200,1200)
    can.Divide(1,3)
    can.cd(1)
    allLumiPerBX.SetMaximum(10)
    allLumiPerBX.SetMinimum(-1)
    allLumiPerBX.SetTitle("Raw Data")
    allLumiPerBX.Draw("hist")
    can.cd(2)
    allCorrLumiPerBX_raw.SetTitle("undo AfterGlow corrections")
    allCorrLumiPerBX_raw.SetMaximum(10)
    allCorrLumiPerBX_raw.SetMinimum(-1)
    allCorrLumiPerBX_raw.Draw("hist")
    can.cd(3)
    allLumiPerBX_newSBR.SetTitle("+New AfterGlow corrections")
    allLumiPerBX_newSBR.SetMaximum(10)
    allLumiPerBX_newSBR.SetMinimum(-1)
    allLumiPerBX_newSBR.Draw("hist")
    input()

ROOT.gStyle.SetOptStat(0)
if args.batch is True:
    ROOT.gROOT.SetBatch(ROOT.kTRUE)

# Demo
#f='/brildata/22/8057/8057_356308_2207271707_2207271732.hd5'
#lumi = hdf_to_pandas(f, 'hfetlumi')
#print(lumi[['fillnum','runnum','lsnum','nbnum','avg']])
#lumi['ls:nb'] = lumi['lsnum'].astype(str) + ':' + lumi['nbnum'].astype(str)
#beam = hdf_to_pandas(f, 'beam')
#print(beam[['fillnum','runnum','lsnum','nbnum','collidable']])
#beam['ls:nb'] = beam['lsnum'].astype(str) + ':' + beam['nbnum'].astype(str)
#df = lumi.set_index('ls:nb').join(beam.set_index('ls:nb').collidable).dropna()
#print(df[['fillnum','runnum','lsnum','nbnum','collidable','avg']])
#df = df.reset_index().drop(columns=['ls:nb'])[['fillnum','runnum','lsnum','nbnum','collidable','avg']]
#print(df)
#exit()

direc_additional='/brildata/hf_reprocess/hfet/'
if args.method=='hfoc':
    direc_additional='/brildata/hf_reprocess/hfoc/'
if args.fills!='':
    if int(args.fills)>=8115:
        direc_additional='/brildata/22/'
direc='/brildata/22/'
filename=args.file
data=''
if args.fills!="":
    fills=args.fills.split(",")
    data=[]
    for fill in fills:
        if os.path.isdir(direc_additional+fill):
            for f in pathlib.Path(direc_additional+'/'+fill+'/').glob('*.hd5'):
                lumitemp=hdf_to_pandas(f, args.method+'lumi')
                beamtemp=[]
                if direc_additional=='/brildata/hf_reprocess/hfoc/' or direc_additional=='/brildata/hf_reprocess/hfet/':
                    print(str(f).split('/')[-1].split('.hd5')[0])
                    for g in pathlib.Path(direc+'/'+fill+'/').glob('*.hd5'):
                        if str(f).split('/')[-1].split('.hd5')[0] in str(g):
                            print(g)
                            beamtemp=hdf_to_pandas(g, 'beam')
                else:
                    beamtemp=hdf_to_pandas(f, 'beam')
                #print('lumi ',lumitemp)
                #print('beam ',beamtemp)
                if len(lumitemp)<=0 or len(beamtemp)<=0: continue
                lumitemp['fill:ls:nb'] = lumitemp['fillnum'].astype(str) + ':' +lumitemp['lsnum'].astype(str) + ':' + lumitemp['nbnum'].astype(str)
                #print(lumitemp[['fillnum','runnum','lsnum','nbnum','avg']])
                beamtemp['fill:ls:nb'] = beamtemp['fillnum'].astype(str) + ':' +beamtemp['lsnum'].astype(str) + ':' + beamtemp['nbnum'].astype(str)
                #print(beamtemp[['fillnum','runnum','lsnum','nbnum','collidable']])
                datatemp = lumitemp.set_index('fill:ls:nb').join(beamtemp.set_index('fill:ls:nb').collidable).dropna()
                #print(datatemp[['fillnum','runnum','lsnum','nbnum','avg','collidable']])
                datatemp = datatemp.reset_index().drop(columns=['fill:ls:nb']) #[['fillnum','runnum','lsnum','nbnum','collidable','avg']]
                #print(datatemp[['fillnum','runnum','lsnum','nbnum','avg','collidable']])
                if len(data)==0: 
                    data=datatemp
                else: data =  pd.concat([data,datatemp])
    #print(data)#[['fillnum','runnum','lsnum','nbnum','avg','collidable']])
    #lumi = pd.concat([read_fill(direc_additional,f,args.method+'lumi') for f in fills if os.path.isdir(direc_additional+f)])
    #beam = pd.concat([read_fill(direc,f,'beam') for f in fills if os.path.isdir(direc+f)])
    #lumi['fill:ls:nb'] = lumi['fillnum'].astype(str) + ':' +lumi['lsnum'].astype(str) + ':' + lumi['nbnum'].astype(str)
    #print(lumi[['fillnum','runnum','lsnum','nbnum','avg']])
    #beam['fill:ls:nb'] = beam['fillnum'].astype(str) + ':' +beam['lsnum'].astype(str) + ':' + beam['nbnum'].astype(str)
    #print(beam[['fillnum','runnum','lsnum','nbnum','collidable']])
    #data = lumi.set_index('fill:ls:nb').join(beam.set_index('fill:ls:nb').collidable).dropna()
    #print(data[['fillnum','runnum','lsnum','nbnum','avg','collidable']])
    #data = data.reset_index().drop(columns=['fill:ls:nb']) #[['fillnum','runnum','lsnum','nbnum','collidable','avg']]
    #print(data[['fillnum','runnum','lsnum','nbnum','avg','collidable']])

else:
    fills=[]
    #print('only a single file: ', filename)
    lumi=hdf_to_pandas(filename,args.method+'lumi')
    beam=hdf_to_pandas(filename,'beam')
    lumi['fill:ls:nb'] = lumi['fillnum'].astype(str) + ':' +lumi['lsnum'].astype(str) + ':' + lumi['nbnum'].astype(str)
    #print(lumi[['fillnum','runnum','lsnum','nbnum','avg']])
    beam['fill:ls:nb'] = beam['fillnum'].astype(str) + ':' +beam['lsnum'].astype(str) + ':' + beam['nbnum'].astype(str)
    #print(beam[['fillnum','runnum','lsnum','nbnum','collidable']])
    data = lumi.set_index('fill:ls:nb').join(beam.set_index('fill:ls:nb').collidable).dropna()
    #print(data[['fillnum','runnum','lsnum','nbnum','avg','collidable']])
    data = data.reset_index().drop(columns=['fill:ls:nb']) #[['fillnum','runnum','lsnum','nbnum','collidable','avg']]

label=args.label
if len(data)<=0: 
    print('no events left from total')
    exit()
data_info=data[['fillnum','runnum','lsnum','nbnum']]
data_info=data_info.to_numpy()
data_bx=data[['bxraw']]
data_bx=data_bx.to_numpy()
data_bx = np.vstack(data_bx[:,0])
data_bx=data_bx*scale
data_mask=data[['collidable']]
data_mask=data_mask.to_numpy()
data_mask = np.vstack(data_mask[:,0])
nevents=len(data_info)
print('Number of events: ', len(data_info))
print(np.shape(data_bx))
print(np.shape(data_mask))

# -------------------- Prepare histograms 
maxLSInRun={}
runs=[]
runtoFills={}
if args.auto:
    for iev in range(len(data_info)):
        if str(data_info[iev][1]) not in runs:
            print("Adding run ",data_info[iev][1])
            runs.append(str(data_info[iev][1]))
            runtoFills[str(data_info[iev][1])] = str(data_info[iev][0])
        if not str(data_info[iev][1]) in  maxLSInRun:
            maxLSInRun[str(data_info[iev][1])] = data_info[iev][2]
        elif maxLSInRun[str(data_info[iev][1])]<data_info[iev][2]:
            maxLSInRun[str(data_info[iev][1])]=data_info[iev][2]
    #print("auto: (runs, runtoFills, maxLSInRun)  -> ", runs, runtoFills, maxLSInRun)
allLumiPerBX_raw={}
allLumiPerBX_oldSBR={}
allLumiPerBX_newSBR={}
allLumiPerBX_oldSBR_cor={}
allLumiPerBX_newSBR_cor={}
lumimask={}
LBKeys=[]
hist_info={}
for run in runs:
    for iLB in range(int(maxLSInRun[run]/args.nLSInLumiBlock+1)):
        if maxLSInRun[run]<10:
            continue
        LBKey=run+"_LS"+str(iLB*args.nLSInLumiBlock+1)+"_LS"+str((iLB+1)*args.nLSInLumiBlock)+"_Fill"+runtoFills[run]
        if iLB==maxLSInRun[run]/args.nLSInLumiBlock:
            LBKey=run+"_LS"+str(iLB*args.nLSInLumiBlock+1)+"_LS"+str(maxLSInRun[run])+"_Fill"+runtoFills[run]
        LBKeys.append(LBKey)
        #allLumiPerBX_raw[LBKey]=ROOT.TH1F("allLumiPerBX_raw"+LBKey, "", BXLength, 0, BXLength)
        allLumiPerBX_oldSBR[LBKey]=ROOT.TH1F("allLumiPerBX_oldSBR"+LBKey, "", BXLength, 0, BXLength)
        #allLumiPerBX_newSBR[LBKey]=ROOT.TH1F("allLumiPerBX_newSBR"+LBKey, "", BXLength, 0, BXLength)
        allLumiPerBX_oldSBR_cor[LBKey]=ROOT.TH1F("allLumiPerBX_oldSBR_cor"+LBKey, "", BXLength, 0, BXLength)
        #allLumiPerBX_newSBR_cor[LBKey]=ROOT.TH1F("allLumiPerBX_newSBR_cor"+LBKey, "", BXLength, 0, BXLength)
        hist_info[LBKey] = ROOT.TH1F("hist_info"+LBKey, "", 25, 0, 25)
        lumimask[LBKey]=ROOT.TH1F("lumimask"+LBKey, "", BXLength, 0, BXLength)
        # Average lumi, raw, old SBR, new SBR; old residual type1 cor, new residual type1 cor

LBKeys.sort()
#print(LBKeys)
LBKeys2npindex={}
# -------------------- first scan of events, know the start and end of Lumi blocks
for iev in range(nevents):
    if iev%5000==101:
        print("event,fill,run,ls,nb,avg ", iev, data_info[iev][:])
    if str(data_info[iev][1]) not in runs:
        continue
    iLB = int( (data_info[iev][2]-1)/args.nLSInLumiBlock)
    LBKey = str(data_info[iev][1])+"_LS"+str(int(iLB*args.nLSInLumiBlock)+1)+"_LS"+str(int((iLB+1)*args.nLSInLumiBlock))+"_Fill"+str(data_info[iev][0])
    if iLB==(maxLSInRun[str(data_info[iev][1])]/args.nLSInLumiBlock):
        LBKey=str(data_info[iev][1])+"_LS"+str(int(iLB*args.nLSInLumiBlock)+1)+"_LS"+str(int(maxLSInRun[str(data_info[iev][1])]))+"_Fill"+str(data_info[iev][0])
    if not LBKey in LBKeys:
        continue
    if LBKey not in LBKeys2npindex: LBKeys2npindex[LBKey]=[iev,iev]
    elif iev> LBKeys2npindex[LBKey][1]:LBKeys2npindex[LBKey][1]=iev
    elif iev< LBKeys2npindex[LBKey][0]:LBKeys2npindex[LBKey][0]=iev

# -------------------- Fill histograms according to the Lumi blocks
newfile=ROOT.TFile("Overall_Correction_"+label+".root", "recreate")
#for LBKey in LBKeys2npindex.keys():
#    rawdata = data_bx[LBKeys2npindex[LBKey][0]:LBKeys2npindex[LBKey][1]+1]
#    rawdata_mask = data_mask[LBKeys2npindex[LBKey][0]:LBKeys2npindex[LBKey][1]+1]
#    totallumi = np.average( np.multiply(rawdata,rawdata_mask) , axis=1)
#    print(rawdata.shape, rawdata_mask.shape)
#    print(totallumi, np.average(totallumi) )
#    stable = np.where(totallumi>np.average(totallumi)*0.8)
#    print(stable)
#    totallumi=totallumi[stable]
#    g=ROOT.TGraph()
#    for i in range(min(len(totallumi),5000)):
#        g.SetPoint(i, i, totallumi[i])
#    g.SetMarkerSize(2)
#    g.SetMarkerColor(ROOT.kBlue)
#    g.Draw("AP")
#    input()
#
#exit()
#for LBKey in LBKeys2npindex.keys():
#    rawdata = np.average( data_bx[LBKeys2npindex[LBKey][0]:LBKeys2npindex[LBKey][1]+1], axis=0)
#    rawdata_mask = data_mask[LBKeys2npindex[LBKey][0]]
#    rawdata_mask = rawdata_mask.astype('int32')
#    rawdata=np.nan_to_num(rawdata, nan=0, posinf=0)  # remove nan
#    nactive=np.sum(rawdata_mask)
#    print('average,total: ', np.sum(np.multiply(rawdata,rawdata_mask))/nactive, np.sum(np.multiply(rawdata,rawdata_mask))*(LBKeys2npindex[LBKey][1]+1-LBKeys2npindex[LBKey][0]))
#    for i in range(LBKeys2npindex[LBKey][1]+1-LBKeys2npindex[LBKey][0]):
#        rawdata = data_bx[LBKeys2npindex[LBKey][0]+i]
#        rawdata_mask = data_mask[LBKeys2npindex[LBKey][0]]
#        rawdata_mask = rawdata_mask.astype('int32')
#        rawdata=np.nan_to_num(rawdata, nan=0, posinf=0)  # remove nan
#        nactive=np.sum(rawdata_mask)
#        print(i, ' : ', np.sum(np.multiply(rawdata,rawdata_mask))/nactive, np.sum(np.multiply(rawdata,rawdata_mask))*(LBKeys2npindex[LBKey][1]+1-LBKeys2npindex[LBKey][0]))
#    input()
#
#exit()

for LBKey in LBKeys2npindex.keys():
    if verbose>10: print(LBKey)
    rawdata = data_bx[LBKeys2npindex[LBKey][0]:LBKeys2npindex[LBKey][1]+1]
    rawdata_mask = data_mask[LBKeys2npindex[LBKey][0]:LBKeys2npindex[LBKey][1]+1]
    rawdata=np.nan_to_num(rawdata, nan=0, posinf=0)  # remove nan
    averagelumi = np.average( np.multiply(rawdata,rawdata_mask) , axis=1)
    unstable_nbs = np.where(averagelumi<np.average(averagelumi)*0.7)
    if len(unstable_nbs[0])>0: 
        rawdata=np.delete(rawdata, unstable_nbs, 0)
    if len(rawdata)<=0:
        print('no events left after delete')
        continue
    rawdata=np.nan_to_num(rawdata, nan=0, posinf=0)  # remove nan
    rawdata=np.average(rawdata, axis=0)
    rawdata_mask = data_mask[LBKeys2npindex[LBKey][0]]
    rawdata_mask = rawdata_mask.astype('int32')
    nactive=np.sum(rawdata_mask)
    if nactive<=0:
        print("nactive is ", nactive)
        continue
    if verbose>10: print("LBKey ", LBKey)
    if args.undoAG:
        # undo afterglow
        if args.plot:
            htempraw=ROOT.TH1F("htempraw", "", BXLength, 0, BXLength)
            htempraw.SetContent(array.array('d',[0])+array.array('d',rawdata))
        currentSBR='HFSBR_ET_run2.txt'
        newSBR='HFSBR_ET_22v0.txt'
        currentSBR='HFSBR_ET_22v0.txt'
        newSBR='HFSBR_ET_22v11.txt'
        if args.method=='hfoc':
            currentSBR='HFSBR_OC_run2.txt'
            newSBR='HFSBR_OC_22v0.txt'
        text_file = open(currentSBR, "r")
        HFSBR = text_file.read().split(',')
        text_file.close()
        for i in range(len(HFSBR)):
            HFSBR[i]=float(HFSBR[i])
        HFSBR=np.array(HFSBR)
        HFSBR[0]=0
        text_file = open(newSBR, "r")
        newHFSBR = text_file.read().split(',')
        text_file.close()
        for i in range(len(newHFSBR)):
            newHFSBR[i]=float(newHFSBR[i])
        newHFSBR=np.array(newHFSBR)
        newHFSBR[0]=0
        simpletest=False
        activeBXMask_list=np.argwhere(rawdata_mask==1)[0]
        if simpletest:
            rawdata_test = np.copy(rawdata)
            for ibx in activeBXMask_list:
                rawdata -= rawdata[ibx]*np.roll(HFSBR, ibx)
            for ibx in range(len(activeBXMask_list)):
                rawdata += rawdata[activeBXMask_list[-1*ibx-1]]*np.roll(HFSBR, activeBXMask_list[-1*ibx-1])
            for ibx in range(len(rawdata_test)):
                if (rawdata_test[ibx]-rawdata[ibx])>1e-6:
                    print('oops, why? ')
                    exit()
        # Undo afterglow correction
        for ibx in range(len(activeBXMask_list)):
            rawdata += rawdata[activeBXMask_list[-1*ibx-1]]*np.roll(HFSBR, activeBXMask_list[-1*ibx-1])
        for ibx in activeBXMask_list:
            rawdata -= rawdata[ibx]*np.roll(newHFSBR, ibx)
    train_ends=np.roll(rawdata_mask,1)-rawdata_mask
    lumi_plus1=rawdata[np.where(train_ends==1)]
    lumi=np.roll(rawdata,1)[np.where(train_ends==1)]
    type1_res_frac=np.divide(lumi_plus1,lumi)
    meanType1_frac = np.mean(type1_res_frac)
    meanType1_abs = np.mean(lumi_plus1)
    meanType1_frac_e = np.std(type1_res_frac)
    meanType1_abs_e = np.std(lumi_plus1)
    meanType2_frac=[]
    meanType2_abs=[]
    for iend in np.where(train_ends==1)[0]:
        for t2bx in range(30):
            if rawdata_mask[iend+t2bx+1]==1:
                break
            meanType2_abs.append(rawdata[iend+t2bx+1])
            meanType2_frac.append(rawdata[iend+t2bx+1]/rawdata[iend-1])
    # ---------------------------  Local correction: Noise + type1
    # Use the first 50 BX avrgae as noise, A naive noise
    noise=np.mean(np.multiply(rawdata,-1*rawdata_mask+1)[3500:3550])
    if verbose>10: print("The average noise is ", noise)
    rawdata_cor=np.copy(rawdata)-noise
    # Do type1 correction
    lumi_plus1_cor=rawdata_cor[np.where(train_ends==1)]
    lumi_cor=np.roll(rawdata_cor,1)[np.where(train_ends==1)]
    type1_res_frac_cor =np.divide(lumi_plus1_cor,lumi_cor)
    meanType1_frac_cor = np.mean(type1_res_frac_cor)
    meanType1_abs_cor  = np.mean(lumi_plus1_cor)
    meanType1_frac_cor_e = np.std(type1_res_frac_cor)
    meanType1_abs_cor_e  = np.std(lumi_plus1_cor)
    meanType2_frac_cor=[]
    meanType2_abs_cor=[]
    for iend in np.where(train_ends==1)[0]:
        for t2bx in range(30):
            if rawdata_mask[iend+t2bx+1]==1:
                break
            meanType2_abs_cor.append(rawdata_cor[iend+t2bx+1])
            meanType2_frac_cor.append(rawdata_cor[iend+t2bx+1]/rawdata_cor[iend-1])
    # save raw - noise
    hist_info[LBKey].SetBinContent(1,np.sum(np.multiply(rawdata_cor,rawdata_mask))/nactive) #Average SBIL
    hist_info[LBKey].SetBinContent(2,np.sum(np.multiply(rawdata_cor,rawdata_mask))*(LBKeys2npindex[LBKey][1]+1-LBKeys2npindex[LBKey][0]-len(unstable_nbs[0])))  # Total L
    hist_info[LBKey].SetBinContent(3,meanType1_frac_cor) # type1 fraction
    hist_info[LBKey].SetBinContent(4,meanType1_abs_cor)  # type1 Hz/ub
    hist_info[LBKey].SetBinError(3,meanType1_frac_cor_e) # type1 fraction
    hist_info[LBKey].SetBinError(4,meanType1_abs_cor_e)  # type1 Hz/ub
    hist_info[LBKey].SetBinContent(5,np.mean(meanType2_frac_cor) ) # type2 fraction
    hist_info[LBKey].SetBinContent(6,np.mean(meanType2_abs_cor))  # type2 Hz/ub
    hist_info[LBKey].SetBinContent(7,np.std(meanType2_frac_cor)) # type2 fraction RMS
    if verbose>10: print('raw-noise Lumi/average, meanType1_frac, meanType1_abs,meanType2_frac,meanType2_abs ', np.sum(np.multiply(rawdata_cor,rawdata_mask))/nactive, meanType1_frac_cor, meanType1_abs_cor, np.mean(meanType2_frac_cor), np.mean(meanType2_abs_cor))

    rawdata_cor=rawdata_cor-np.roll(np.multiply(rawdata_cor,rawdata_mask)*meanType1_frac_cor,1)
    # repeat type1 residual calculation
    lumi_plus1_cor=rawdata_cor[np.where(train_ends==1)]
    lumi_cor=np.roll(rawdata_cor,1)[np.where(train_ends==1)]
    type1_res_frac_cor=np.divide(lumi_plus1_cor,lumi_cor)
    meanType1_frac_cor = np.mean(type1_res_frac_cor)
    meanType1_abs_cor = np.mean(lumi_plus1_cor)
    meanType1_frac_cor_e = np.std(type1_res_frac_cor)
    meanType1_abs_cor_e = np.std(lumi_plus1_cor)
    meanType2_frac_cor=[]
    meanType2_abs_cor=[]
    for iend in np.where(train_ends==1)[0]:
        for t2bx in range(30):
            if rawdata_mask[iend+t2bx+1]==1:
                break
            meanType2_abs_cor.append(rawdata_cor[iend+t2bx+1])
            meanType2_frac_cor.append(rawdata_cor[iend+t2bx+1]/rawdata_cor[iend-1])

    if verbose>10: print('raw Lumi/average, meanType1_frac, meanType1_abs,meanType2_frac,meanType2_abs ', np.sum(np.multiply(rawdata,rawdata_mask))/nactive, meanType1_frac, meanType1_abs,np.mean(meanType2_frac),np.std(meanType2_abs))
    if verbose>10: print('cor Lumi/average, meanType1_frac, meanType1_abs,meanType2_frac,meanType2_abs ', np.sum(np.multiply(rawdata_cor,rawdata_mask))/nactive, meanType1_frac_cor, meanType1_abs_cor,np.mean(meanType2_frac_cor),np.std(meanType2_abs_cor))
    # ---------------------------  Bunch train effect 
    train_starts=np.roll(rawdata_mask,1) - rawdata_mask
    n_train=0
    sum_lead=0
    sum_train=0
    AveSBIL_leadBX=0
    for istart in np.where(train_starts==-1)[0]:
        if (rawdata_mask[istart+1:istart+5]==1).all():
                sum_lead += rawdata_cor[istart]
                sum_train += np.mean(rawdata_cor[istart+1:istart+5])
                n_train+=1
    if n_train>0: AveSBIL_leadBX = float(sum_lead)/n_train
    #print(n_train, sum_lead, sum_train, AveSBIL_leadBX)
    # ---------------------------  Undo afterglow, redo another afterglow

    # ---------------------------  Save results
    lumimask[LBKey].SetContent(array.array('d',[0])+array.array('d',rawdata_mask))
    allLumiPerBX_oldSBR[LBKey].SetContent(array.array('d',[0])+array.array('d',rawdata)) # start from Bin 1
    allLumiPerBX_oldSBR_cor[LBKey].SetContent(array.array('d',[0])+array.array('d',rawdata_cor)) # start from Bin 1
    if args.plot:
        print(np.max(rawdata))
        htemp=ROOT.TH1F("htemp", "", BXLength, 0, BXLength)
        htemp.SetContent(array.array('d',[0])+array.array('d',rawdata_mask))
        #makeplot(allLumiPerBX_oldSBR[LBKey], allLumiPerBX_oldSBR_cor[LBKey], htemp)
        makeplots(htempraw, allLumiPerBX_oldSBR[LBKey], allLumiPerBX_oldSBR_cor[LBKey])
    hist_info[LBKey].SetBinContent(8,np.sum(np.multiply(rawdata,rawdata_mask))/nactive) #Average SBIL
    hist_info[LBKey].SetBinContent(9,np.sum(np.multiply(rawdata,rawdata_mask))*(LBKeys2npindex[LBKey][1]+1-LBKeys2npindex[LBKey][0]-len(unstable_nbs[0])))  # Total L
    hist_info[LBKey].SetBinContent(10,meanType1_frac) # type1 fraction
    hist_info[LBKey].SetBinContent(11,meanType1_abs)  # type1 Hz/ub
    hist_info[LBKey].SetBinError(10,meanType1_frac_e) # type1 fraction
    hist_info[LBKey].SetBinError(11,meanType1_abs_e)  # type1 Hz/ub
    hist_info[LBKey].SetBinContent(12,np.mean(meanType2_frac) ) # type2 fraction
    hist_info[LBKey].SetBinContent(13,np.mean(meanType2_abs))  # type2 Hz/ub
    hist_info[LBKey].SetBinContent(14,np.std(meanType2_frac)) # type2 fraction RMS

    hist_info[LBKey].SetBinContent(15,np.sum(np.multiply(rawdata_cor,rawdata_mask))/nactive) #Average SBIL
    hist_info[LBKey].SetBinContent(16,np.sum(np.multiply(rawdata_cor,rawdata_mask))*(LBKeys2npindex[LBKey][1]+1-LBKeys2npindex[LBKey][0]-len(unstable_nbs[0])))  # Total L
    hist_info[LBKey].SetBinContent(17,meanType1_frac_cor) # type1 fraction
    hist_info[LBKey].SetBinContent(18,meanType1_abs_cor)  # type1 Hz/ub
    hist_info[LBKey].SetBinError(17,meanType1_frac_cor_e) # type1 fraction
    hist_info[LBKey].SetBinError(18,meanType1_abs_cor_e)  # type1 Hz/ub
    hist_info[LBKey].SetBinContent(19,np.mean(meanType2_frac_cor) ) # type2 fraction
    hist_info[LBKey].SetBinContent(20,np.mean(meanType2_abs_cor))  # type2 Hz/ub
    hist_info[LBKey].SetBinContent(21,np.std(meanType2_frac_cor)) # type2 fraction RMS
    hist_info[LBKey].SetBinContent(22,np.sum(rawdata_mask)) # Nactive BX
    hist_info[LBKey].SetBinContent(23,sum_lead)
    hist_info[LBKey].SetBinContent(24,sum_train)
    hist_info[LBKey].SetBinContent(25,AveSBIL_leadBX) 
    newfile.WriteTObject(allLumiPerBX_oldSBR[LBKey], "allLumiPerBX_oldSBR_"+LBKey)
    newfile.WriteTObject(allLumiPerBX_oldSBR_cor[LBKey], "allLumiPerBX_oldSBR_cor_"+LBKey)
    newfile.WriteTObject(hist_info[LBKey], "hist_info_"+LBKey)
    newfile.WriteTObject(lumimask[LBKey], "lumimask_"+LBKey)







