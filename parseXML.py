import xml.etree.ElementTree as ET
import binascii
tree = ET.parse('LUT-dataRun3_HLTNew_HCAL.xml')

crates = [22, 29, 32]

root = tree.getroot()

average_perboard=True
if average_perboard:
    board_LUT_list={}
    for channel in root :
        crate, slot, fiber, fiberChan, iEta, iPhi, depth, lutType = None, None, None, None, None, None, None, None
        goodLine = False
        for prop in channel :
            if prop.tag == 'Parameter' :
                if prop.attrib['name'] == 'CRATE' : crate = int(prop.text)
                if prop.attrib['name'] == 'SLOT' :   slot = int(prop.text)
                if prop.attrib['name'] == 'FIBER' :  fiber = int(prop.text)
                if prop.attrib['name'] == 'FIBERCHAN' :  fiberChan = int(prop.text)
                if prop.attrib['name'] == 'IETA' :  iEta = int(prop.text)
                if prop.attrib['name'] == 'IPHI' :  iPhi = int(prop.text)
                if prop.attrib['name'] == 'DEPTH' : depth = int(prop.text)
                if prop.attrib['name'] == 'LUT_TYPE' : lutType = int(prop.text)
            elif (prop.tag == 'Data') and (lutType == 1) and depth and (crate in crates) :
                if not (( (fiber==0 or fiber==3 or fiber==12 or fiber==15) and (fiberChan==2 or fiberChan==3)) or ( (fiber==2 or fiber==5 or fiber==14 or fiber==17) and (fiberChan==0 or fiberChan==1))): continue
                outLine = "{0:d} {1:d} {2:d} {3:d} {4:d} {5:d} {6:d}".format(crate, slot, fiber, fiberChan, iEta, iPhi, depth)
                vals = prop.text.split()
                if 'crate_%s_slot_%s'%(str(crate),str(slot)) not in board_LUT_list.keys():
                    board_LUT_list['crate_%s_slot_%s'%(str(crate),str(slot))]={}
                    board_LUT_list['crate_%s_slot_%s'%(str(crate),str(slot))]['entry'] = 1
                    board_LUT_list['crate_%s_slot_%s'%(str(crate),str(slot))]['lut'] = []
                    for i in range(len(vals)) :
                        hex2decimal2=int('0x%s'%(vals[i]),16) & 0x7ff
                        board_LUT_list['crate_%s_slot_%s'%(str(crate),str(slot))]['lut'].append(hex2decimal2)
                else:
                    board_LUT_list['crate_%s_slot_%s'%(str(crate),str(slot))]['entry'] += 1
                    for i in range(len(vals)) :
                        hex2decimal2=int('0x%s'%(vals[i]),16) & 0x7ff
                        board_LUT_list['crate_%s_slot_%s'%(str(crate),str(slot))]['lut'][i] += hex2decimal2
    for iboard in board_LUT_list.keys():
        for i in range(len(board_LUT_list[iboard]['lut'])):
            board_LUT_list[iboard]['lut'][i]/=board_LUT_list[iboard]['entry']


outLines = []
for channel in root :
    crate, slot, fiber, fiberChan, iEta, iPhi, depth, lutType = None, None, None, None, None, None, None, None
    goodLine = False 
    for prop in channel :
        if prop.tag == 'Parameter' :
            if prop.attrib['name'] == 'CRATE' : crate = int(prop.text)
            if prop.attrib['name'] == 'SLOT' :   slot = int(prop.text)
            if prop.attrib['name'] == 'FIBER' :  fiber = int(prop.text)
            if prop.attrib['name'] == 'FIBERCHAN' :  fiberChan = int(prop.text)  
            if prop.attrib['name'] == 'IETA' :  iEta = int(prop.text)
            if prop.attrib['name'] == 'IPHI' :  iPhi = int(prop.text)
            if prop.attrib['name'] == 'DEPTH' : depth = int(prop.text)
            if prop.attrib['name'] == 'LUT_TYPE' : lutType = int(prop.text)
        elif (prop.tag == 'Data') and (lutType == 1) and depth and (crate in crates) :
            goodLine = True 
            if not (( (fiber==0 or fiber==3 or fiber==12 or fiber==15) and (fiberChan==2 or fiberChan==3)) or ( (fiber==2 or fiber==5 or fiber==14 or fiber==17) and (fiberChan==0 or fiberChan==1))): goodLine = False  # only selected channels 
            outLine = "{0:d} {1:d} {2:d} {3:d} {4:d} {5:d} {6:d}".format(crate, slot, fiber, fiberChan, iEta, iPhi, depth)
            vals = prop.text.split()
            for i in range(len(vals)) : 
                if average_perboard:
                    outLine += " {0:d}".format(board_LUT_list['crate_%s_slot_%s'%(str(crate),str(slot))]['lut'][i])
                else:
                    hex2decimal2=int('0x%s'%(vals[i]),16) & 0x7ff
                    outLine += " {0:d}".format(hex2decimal2)
    if goodLine : outLines.append(outLine + "\n")

with open('outLUT_2022.txt','w') as f :
    for outLine in outLines :
        f.write(outLine)

f.close()
    
        

    


