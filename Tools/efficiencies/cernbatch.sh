# Lxplus Batch Job Script

# Run:$ bsub -R "pool>30000" -q 1nd -J MuonISO < cernbatch.sho

config="config/configQCD_pT470to600.ini"
outputname="TestAll2"

export CMSSW_PROJECT_SRC="MuonIsolation-902X/CMSSW_9_0_2_patch1/src"
export EXE_FILE="MuonPOGtreeProducer/Tools/efficiencies/ISOEfficiency.run"
export TOP="$PWD"

echo /afs/cern.ch/work/b/brochero/$CMSSW_PROJECT_SRC
cd /afs/cern.ch/work/b/brochero/$CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
cd $TOP
/afs/cern.ch/work/b/brochero/$CMSSW_PROJECT_SRC/$EXE_FILE /afs/cern.ch/work/b/brochero/$CMSSW_PROJECT_SRC/MuonPOGtreeProducer/Tools/efficiencies/${config} /afs/cern.ch/work/b/brochero/$CMSSW_PROJECT_SRC/MuonPOGtreeProducer/Tools/efficiencies/MuonResults ${outputname}
