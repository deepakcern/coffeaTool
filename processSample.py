import awkward as ak
from coffea.nanoevents import NanoEventsFactory,BaseSchema
import uproot
from monoHbbProcessor import monoHbbProcessor
import sys, optparse,argparse
import os

usage = "python DataframeToHist.py -F -inDir directoryName -D outputDir "
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", "--inputfile",  dest="inputfile",default="myfiles.root")
parser.add_argument("-F", "--farmout", action="store_true",  dest="farmout")
parser.add_argument("-m", "--multifiles", action="store_true",  dest="multifiles")
parser.add_argument("-y", "--year", dest="year",default="2017")
parser.add_argument("-D", "--outputdir", dest="outputdir",default=".")
parser.add_argument("-tag", "--tag", dest="tag",default="test")
args = parser.parse_args()

infile  = args.inputfile

runOnFiles = args.multifiles 
year       = args.year
print ("year  : ",year)
tag        = args.tag
outputpath = "/eos/cms/store/group/phys_exotica/monoHiggs/monoHbb/coffeaoutput/"+year+"/"+tag


if runOnFiles:
    if '.root' in infile.split('/')[-1]:fname=infile
    if '.txt' in infile.split('/')[-1]:
        fname=open(infile).readline().rstrip()
else:
    fname='/eos/cms/store/group/phys_exotica/monoHiggs/monoHbb/2018_skimmedFiles/merged_bkgrootfiles_splitted/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_0.root'
    #fname='/eos/cms/store/group/phys_exotica/monoHiggs/monoHbb/2018_skimmedFiles/merged_bkgrootfiles_splitted/MET-Run2018D-PromptReco-v2_0.root'
    #fname='/eos/cms/store/group/phys_exotica/monoHiggs/monoHbb/skimmedFiles/TTTest.root'
    outputpath = '.'
    #fname='/eos/cms/store/group/phys_exotica/monoHiggs/monoHbb/2018_skimmedFiles/v12.09_NoJERdata/MET-Run2018D-PromptReco-v2_146.root'

isData = False
if 'Run' in fname.split('/')[-1] or 'MET' in fname.split('/')[-1] or 'EGamma' in fname.split('/')[-1] or 'SingleElectron' in fname.split('/')[-1]: isData = True
else:isData=False

if 'TTT' in fname.split('/')[-1]:dataset="TT"
else:dataset="other"

if not os.path.isdir(outputpath):os.system('mkdir -p '+outputpath)

'''
----------------------------------
INITIALIZE COFFEA AND GET AWKWARD ARRAY OF EVENTS
---------------------------------
'''
print ("inputfile : ",fname)
f = uproot.open(fname)
totalweightedevents = f["h_total_mcweight"].values()[1]
totalevents         = f["h_total"].values()[1]
filename            = fname.split('/')[-1]


events = NanoEventsFactory.from_root(
    fname,
    schemaclass=BaseSchema,    
    treepath="outTree",
    metadata={"dataset":dataset,"isData":isData,"outputpath":outputpath,"filename":filename,"totalweightevents":float(totalweightedevents),"totalevents":int(totalevents),"year":year}
).events()


'''
----------------------------------
CALL THE PROCESSOR TO RUN OVER THE EVENTS AND SAME OUTPUT
---------------------------------
'''
p = monoHbbProcessor() 
out = p.process(events)
#print (out)
print ('\n'+"================= Process completed =============="+'\n')
print ("OUTPUT PATH :",outputpath)
