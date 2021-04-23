import os

#filepath = "/eos/user/d/daebi/cosmic_MC/" #MC cosmic
#filepath = "/eos/cms/store/express/Commissioning2020/ExpressCosmics/FEVT/Express-v1/000/337/973/00000/" #MWGR4 express dataset
#filepath = "/eos/cms/store/express/Commissioning2020/ExpressCosmics/FEVT/Express-v1/000/338/714/" #MWGR5 express dataset
filepath = "/eos/cms/store/express/Commissioning2021/ExpressCosmics/FEVT/Express-v1/000/341/169/" #MWGR3 express dataset
pwd = "/afs/cern.ch/work/d/daebi/analyser/CMSSW_11_3_0_pre5/src/GEMCSCBendingAnalyzer/MuonAnalyser/condor/test_condor/"

n_files_per_python = 20

#Whether or not to use a misalignment
misalign = False
misalign_db = "gemAl.db"
misalign_tag = "GEM"
misalign_APR_tag = "test"

#Set up the current folder to be the working space
os.system("mkdir -p output")
os.system("mkdir -p error")
os.system("mkdir -p log")
os.system("mkdir -p python")
os.system("mkdir -p scripts")
os.system("mkdir -p rootfiles")

sub_all = open("submit_all.sh", "write")
sub_all.write("#!/bin/bash\n")
sub_all.write("cd "+pwd+"\n")
sub_all.write("eval `scramv1 runtime -sh`\n")

filelist = os.listdir(filepath)
print filelist

filecounter = 0
files_in_python = 0
job_filelist = []

for subdir in filelist:

  #subdir is the first folder of dataset, "0000/" and "0001/" usually
  filepath_tmp = filepath+subdir+"/"
  filelist_tmp = os.listdir(filepath_tmp)

  for name in filelist_tmp:
    print "path = ", filepath_tmp, " name = ", name
    if (".root" not in name):
      continue
    name = name[:-5] #Take filename without '.root'
    print name

    if files_in_python == 0:
      ana_name = "python/analyser_"+name+".py"
      ana_script = open(ana_name, "write")

      c_name = "scripts/condor_"+name+".sh"
      c_script = open(c_name, "write")

    job_filelist.append('file:{filepath}{name}.root'.format(name = name, filepath = filepath_tmp))
    files_in_python += 1    

    if (subdir == filelist[-1] and name == filelist_tmp[-1][:-5]):
      print "WE DONE!"
    if files_in_python == n_files_per_python or (subdir == filelist[-1] and name == filelist_tmp[-1][:-5]):
      files_in_python = 0
      print job_filelist

      ana_script.write("import FWCore.ParameterSet.Config as cms \n")
      ana_script.write("from Configuration.Eras.Era_Run3_cff import Run3\n")
      ana_script.write("process = cms.Process('analyser',Run3)\n")
      ana_script.write('process.load("FWCore.MessageService.MessageLogger_cfi")\n')
      #ana_script.write("process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')\n")
      ana_script.write("process.load('Configuration.StandardSequences.MagneticField_0T_cff')\n")
      ana_script.write("process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')\n")
      ana_script.write("process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')\n")
      ana_script.write("process.load('Configuration.StandardSequences.SimIdeal_cff')\n")
      ana_script.write("process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')\n")
      ana_script.write("process.load('Configuration.StandardSequences.GeometryRecoDB_cff')\n")
      ana_script.write("from Configuration.AlCa.GlobalTag import GlobalTag\n")
      ana_script.write("process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')\n")



      if misalign:
        #This is misalignment part
        ana_script.write("process.GlobalTag.toGet = cms.VPSet(\n")
        ana_script.write("    cms.PSet(\n")
        ana_script.write("            connect = cms.string('sqlite_file:{misalign_db}'),\n".format(misalign_db = pwd+misalign_db))
        ana_script.write("            record = cms.string('GEMAlignmentRcd'),\n")
        ana_script.write("            tag = cms.string('{misalign_tag}')\n".format(misalign_tag = misalign_tag))
        ana_script.write("    ),\n")
        ana_script.write("    cms.PSet(\n")
        ana_script.write("            connect = cms.string('sqlite_file:{misalign_db}'),\n".format(misalign_db = pwd+misalign_db))
        ana_script.write("            record = cms.string('GEMAlignmentErrorExtendedRcd'),\n")
        ana_script.write("            tag = cms.string('{misalign_APR_tag}')\n".format(misalign_APR_tag = misalign_APR_tag))
        ana_script.write("    ),\n")
        ana_script.write("    cms.PSet(record=cms.string('GlobalPositionRcd'), tag = cms.string('IdealGeometry'))\n")
        ana_script.write(")\n")
        ana_script.write("process.GEMGeometryESModule.applyAlignment = cms.bool(True)\n")
        #End of misalignment part

      ana_script.write("process.MessageLogger.cerr.FwkReport.reportEvery = 5000\n")
      ana_script.write("from FWCore.ParameterSet.VarParsing import VarParsing\n")
      ana_script.write("options = VarParsing('analysis')\n")
      ana_script.write("options.register ('nEvents',\n")
      ana_script.write("                        -1, #Max number of events\n")
      ana_script.write("                        VarParsing.multiplicity.singleton,\n")
      ana_script.write("                        VarParsing.varType.int,\n")
      ana_script.write('                        "Number of events")\n')
      ana_script.write("options.parseArguments()\n")
      ana_script.write("process.maxEvents = cms.untracked.PSet(\n")
      ana_script.write("  input = cms.untracked.int32(options.nEvents)\n")
      ana_script.write(")\n")
      ana_script.write("process.maxEvents.input = cms.untracked.int32(-1)\n")
      ana_script.write('process.source = cms.Source("PoolSource",\n')
      ana_script.write("                                fileNames = cms.untracked.vstring({input_filelist}),\n".format(input_filelist = job_filelist))
      ana_script.write("                                inputCommands = cms.untracked.vstring(\n")
      ana_script.write('                        "keep *",\n')
      ana_script.write('                        "drop TotemTimingDigiedmDetSetVector_totemTimingRawToDigi_TotemTiming_reRECO",\n')
      ana_script.write('                        "drop TotemTimingRecHitedmDetSetVector_totemTimingRecHits__reRECO"\n')
      ana_script.write("                        )\n")
      ana_script.write("                                )\n")
      #ana_script.write("process.source.fileNames.append('file:{filepath}{name}.root')\n".format(name = name, filepath = filepath_tmp))
      ana_script.write('process.options = cms.untracked.PSet(\n')
      ana_script.write("                        SkipEvent = cms.untracked.vstring('ProductNotFound')\n")
      ana_script.write("                        )\n")
      ana_script.write('process.TFileService = cms.Service("TFileService", fileName = cms.string("{pwd}rootfiles/out_ana_{name}.root"))\n'.format(name = name, pwd = pwd))
      ana_script.write("process.analyser = cms.EDAnalyzer('analyser',\n")
      ana_script.write('        process.MuonServiceProxy,\n')
      ana_script.write('        gemRecHits = cms.InputTag("gemRecHits"),\n')
      ana_script.write('        gemSimHits = cms.InputTag("g4SimHits", "MuonGEMHits"),\n')
      ana_script.write('        muons = cms.InputTag("muons"),\n')
      ana_script.write('        vertexCollection = cms.InputTag("offlinePrimaryVerticies"),\n')
      ana_script.write('        tracker_prop = cms.bool(False),\n')
      ana_script.write('        CSC_prop = cms.bool(True),\n')
      ana_script.write('        debug = cms.bool(False)\n')
      ana_script.write(')\n')
      ana_script.write('process.p = cms.EndPath(process.analyser)\n')


      c_script.write("#!/bin/bash\n")
      c_script.write("date\n")
      c_script.write("cd "+pwd+"\n")
      c_script.write("eval `scramv1 runtime -sh`\n")
      c_script.write("cmsRun "+ana_name+"\n")
      c_script.write("echo 'done'\n")
      c_script.write("date\n")

      condor_sub = open("scripts/submit_"+name+".sh", "write")
      condor_sub.write("""universe                = vanilla
executable              = {fname}
arguments               = no
output                  = {pwd}/output/out_{name}.$(ClusterId).$(ProcId).out
error                   = {pwd}/error/err_{name}.$(ClusterId).$(ProcId).err
log                     = {pwd}/log/log_{name}.$(ClusterId).log
request_memory          = 4000M
+JobFlavour             = "workday"
queue""".format(fname = c_name, pwd = pwd, name = name))

      sub_all.write("echo File {filecounter}\n".format(filecounter = filecounter))
      sub_all.write("condor_submit scripts/submit_"+name+".sh\n")
    filecounter += 1

    os.system("chmod 755 "+c_name)

os.system("chmod 755 submit_all.sh")


sub_all = open("submit_all_condor.sh", "write")
sub_all.write("""executable              = $(filename)
output                  = output/job_$(filename).out
error                   = error/job_$(filename).err
log                     = lob/job_$(filename).log
universe                = vanilla
queue filename matching files scripts/condor*.sh""")
