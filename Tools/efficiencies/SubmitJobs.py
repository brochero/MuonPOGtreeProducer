#!/usr/bin/env python
import os, time, socket, sys

BaseCMSDir     = os.environ["CMSSW_BASE"]+"/src"
OutputName = str(sys.argv[1])
tqu = "1nd"
#tqu = "8nh"
ConfigDir = "config/config"
Configs    = [
    #"DY",
    # "QCD_pT15to20",
    # "QCD_pT20to30",
    # "QCD_pT30to50",
    # "QCD_pT50to80",
    # "QCD_pT80to120",
    # "QCD_pT120to170",
    # "QCD_pT170to300",
    "QCD_pT300to470",
    ## "QCD_pT470to600",
    #"QCD_pT600to800",
    #"QCD_pT800to1000",
    ## "QCD_pT1000toInf",
    ]

RunFileName    = "ToSubmit.sh"

for cfg in Configs: 
    fout = open(RunFileName, "w")
    print>>fout, 'export CMSSW_PROJECT_SRC="' + BaseCMSDir + '"'
    print>>fout, """
export EXE_FILE="MuonPOGtreeProducer/Tools/efficiencies/ISOEfficiency.run"
export TOP="$PWD"
cd $CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
cd $TOP
"""
    print>>fout, '$CMSSW_PROJECT_SRC/$EXE_FILE $CMSSW_PROJECT_SRC/MuonPOGtreeProducer/Tools/efficiencies/' + ConfigDir + cfg + '.ini $CMSSW_PROJECT_SRC/MuonPOGtreeProducer/Tools/efficiencies/MuonResults ' + cfg + '_' + OutputName
    
    fout = None
    os.chmod(RunFileName,0744)
    command = 'bsub -R "pool>300000" -q ' + tqu + ' -J ' + cfg + ' < ' + RunFileName

    print 'Submitting job with command: '
    print str(command)
    os.system( command )

    
os.system( "rm " + RunFileName )
