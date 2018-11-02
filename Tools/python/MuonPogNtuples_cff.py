import FWCore.ParameterSet.Config as cms


def appendMuonPogNtuple(process, runOnMC, processTag="HLT", ntupleFileName="MuonPogTree.root") :

    process.load("MuonPOGtreeProducer.Tools.MuonPogTreeProducer_cfi")
    process.load("CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi")

    if processTag != "HLT" :
        print "[MuonPogNtuples]: Customising process tag for TriggerResults / Summary to :", processTag
        process.MuonPogTree.TrigResultsTag = "TriggerResults::"+processTag
        process.MuonPogTree.TrigSummaryTag = "hltTriggerSummaryAOD::"+processTag

    if runOnMC :
        # PUPPI Isolation
        process.PUPPIMuonRelIso = cms.EDProducer('PuppiLeptonIsolation',
                                                 srcLepton = cms.string( 'slimmedMuons' ),
                                                 dR = cms.double( 0.4 ),
                                                 mixFraction = cms.double( 0.5 ), 
                                                 configuration = cms.string( "#detail#" )
                                                 )

        #process.load("MuonPOGtreeProducer.Tools.PrunedGenParticles_cfi")
        #process.muonPogNtuple = cms.Sequence(process.prunedGenParticles + process.MuonPogTree)
        #process.muonPogNtuple = cms.Sequence(process.MuonPogTree)        
        # process.out = cms.OutputModule("PoolOutputModule",
        #                                fileName = cms.untracked.string('myOutputFile.root'),
        #                                outputCommands = cms.untracked.vstring(['keep *'])
        #                                )
        # process.muonPogNtuple = cms.Sequence(process.PUPPIMuonRelIso+process.MuonPogTree+process.out)
        process.muonPogNtuple = cms.Sequence(process.PUPPIMuonRelIso+process.MuonPogTree)

    else :
        process.muonPogNtuple = cms.Sequence(process.MuonPogTree)
        process.MuonPogTree.PileUpInfoTag = cms.untracked.InputTag("none")
        process.MuonPogTree.GenInfoTag = cms.untracked.InputTag("none")
        process.MuonPogTree.GenTag = cms.untracked.InputTag("none")
        
    process.TFileService = cms.Service('TFileService',
        fileName = cms.string(ntupleFileName)
    )

    if hasattr(process,"reconstruction_step") :
        print "[MuonPogNtuples]: Appending goodOfflinePrimaryVertices to RECO step"
        process.AOutput.replace(process.reconstruction_step, process.reconstruction_step + process.goodOfflinePrimaryVertices)
    else :
        print "[MuonPogNtuples]: Creating FastFilter path to host goodOfflinePrimaryVertices"
        process.FastFilters = cms.Path(process.goodOfflinePrimaryVertices)
    
    if hasattr(process,"AOutput") :
        print "[MuonPogNtuples]: EndPath AOutput found, appending ntuples"
        process.AOutput.replace(process.hltOutputA, process.hltOutputA + process.muonPogNtuple)
    else :
        print "[MuonPogNtuples]: EndPath AOuptput not found, creating it for ntuple sequence"
        process.AOutput = cms.EndPath(process.muonPogNtuple)

def customiseHlt(process, pathCut = "all", filterCut = "all") :
    if hasattr(process,"MuonPogTree") :
        print "[MuonPogNtuples]: skimming HLT format using:\n" \
            + "\tpaths : " + pathCut + "\n" \
            + "\tfilters : " + filterCut 
            
        process.MuonPogTree.TrigPathCut = pathCut
        process.MuonPogTree.TrigFilterCut = filterCut
    else : 
        print "[MuonPogNtuples]: muonPogTree not found, check your cfg!"

def customiseMuonCuts(process, minMuPt = 0., minNMu = 0) :
    if hasattr(process,"MuonPogTree") :
        print "[MuonPogNtuples]: skimming ntuple saving only muons that are: " \
            + "# STA || TRK || GLB muons >= " + str(minNMu) + " with muon pT > " + str(minMuPt) 
            
        process.MuonPogTree.MinMuPtCut = cms.untracked.double(minMuPt)
        process.MuonPogTree.MinNMuCut  = cms.untracked.int32(minNMu)

        if hasattr(process,"prunedGenParticles") :
            print "[MuonPogNtuples]: applying pT cut to GEN particles as well"
            
            process.prunedGenParticles.select = cms.vstring("drop *"
                                                            , "++keep pdgId =  13 && pt >" + str(minMuPt) 
                                                            , "++keep pdgId = -13 && pt >" + str(minMuPt)
                                                            )

    else :
        print "[MuonPogNtuples]: muonPogTree not found, check your cfg!"

