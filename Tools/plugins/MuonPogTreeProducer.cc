//////////////////////////////////////
// Ntuplizer that fills muon_pog trees
//////////////////////////////////////

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"

#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JetCollection.h"

// MiniIsolation
#include "PhysicsTools/PatUtils/interface/MiniIsolation.h"

#include "MuonPOGtreeProducer/Tools/src/MuonPogTree.h"
#include "TTree.h"

#include <algorithm>
#include <iostream>

class MuonPogTreeProducer : public edm::EDAnalyzer 
{
public:

  MuonPogTreeProducer(const edm::ParameterSet &);
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();
  
private:
  
  void fillGenInfo(const edm::Handle<std::vector<PileupSummaryInfo> > &,
		   const  edm::Handle<GenEventInfoProduct> &);

  void fillGenParticles(const edm::Handle<reco::GenParticleCollection> &);

  void fillHlt(const edm::Handle<edm::TriggerResults> &, 
	       const edm::Handle<trigger::TriggerEvent> &,
	       const edm::TriggerNames &);
  
  void fillPV(const edm::Handle<std::vector<reco::Vertex> > &);
  
  
  Int_t fillMuons(const edm::Handle<edm::View<pat::Muon> > &,
		  const edm::Handle<std::vector<reco::Vertex> > &,
		  const edm::Handle<reco::BeamSpot> &,
		  const edm::Handle<std::vector<pat::PackedCandidate>> &,
		  const edm::Handle<edm::ValueMap<double>>&,
		  const edm::Handle<edm::ValueMap<double>>&,
		  const edm::Handle<edm::ValueMap<double>>&,
		  const edm::Handle<edm::ValueMap<double>>&,
		  const edm::Handle<edm::ValueMap<double>>&,
		  const edm::Handle<edm::ValueMap<double>>&);

  void fillJets(const edm::Handle<edm::View<pat::Jet> > &);

  void fillGenJets(const edm::Handle<edm::View<reco::GenJet> > &);

  void fillL1(const edm::Handle<l1t::MuonBxCollection> &);

    
  // returns false in case the match is for a RPC chamber
  bool getMuonChamberId(DetId & id, muon_pog::MuonDetType & det, Int_t & r, Int_t & phi, Int_t & eta) const ;
  
  edm::EDGetTokenT<edm::TriggerResults> trigResultsToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> trigSummaryToken_;

  std::string trigFilterCut_;
  std::string trigPathCut_;

  edm::EDGetTokenT<edm::View<pat::Muon> > muonToken_;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetsToken_;
  edm::EDGetTokenT<edm::View<reco::GenJet> > genjetsToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > primaryVertexToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;

  edm::EDGetTokenT<reco::PFMETCollection> pfMetToken_;
  edm::EDGetTokenT<reco::PFMETCollection> pfChMetToken_;
  edm::EDGetTokenT<reco::CaloMETCollection> caloMetToken_;

  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileUpInfoToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;

  edm::EDGetTokenT<LumiScalersCollection> scalersToken_;
    
  edm::EDGetTokenT<l1t::MuonBxCollection> l1Token_;

  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > pfCandsToken_;

  edm::EDGetTokenT<edm::ValueMap<double>> PUPPIToken_,PUPPICHToken_,PUPPINHToken_,PUPPIPHToken_;
  edm::EDGetTokenT<edm::ValueMap<double>> PUPPILepToken_;
  edm::EDGetTokenT<edm::ValueMap<double>> PUPPINoLepToken_;

  Float_t m_minMuPtCut;
  Int_t m_minNMuCut;

  muon_pog::Event event_;
  muon_pog::EventId eventId_;
  std::map<std::string,TTree*> tree_;
  
};


MuonPogTreeProducer::MuonPogTreeProducer( const edm::ParameterSet & cfg )
{

  // Input collections
  edm::InputTag tag = cfg.getUntrackedParameter<edm::InputTag>("TrigResultsTag", edm::InputTag("TriggerResults::HLT"));
  if (tag.label() != "none") trigResultsToken_ = consumes<edm::TriggerResults>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("TrigSummaryTag", edm::InputTag("hltTriggerSummaryAOD::HLT")); 
  if (tag.label() != "none") trigSummaryToken_ =consumes<trigger::TriggerEvent>(tag);

  trigFilterCut_ = cfg.getUntrackedParameter<std::string>("TrigFilterCut", std::string("all"));
  trigPathCut_ = cfg.getUntrackedParameter<std::string>("TrigPathCut", std::string("all"));

  tag = cfg.getUntrackedParameter<edm::InputTag>("MuonTag", edm::InputTag("muons"));
  if (tag.label() != "none") muonToken_ = consumes<edm::View<pat::Muon> >(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("JetTag", edm::InputTag("ak4PFJetsCHS"));
  //tag = cfg.getUntrackedParameter<edm::InputTag>("JetTag", edm::InputTag("slimmedJets"));
  if (tag.label() != "none") jetsToken_ = consumes<edm::View<pat::Jet> >(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("GenJetTag", edm::InputTag("ak4GenJets"));
  if (tag.label() != "none") genjetsToken_ = consumes<edm::View<reco::GenJet> >(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("PrimaryVertexTag", edm::InputTag("offlinePrimaryVertices"));
  if (tag.label() != "none") primaryVertexToken_ = consumes<std::vector<reco::Vertex> >(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("BeamSpotTag", edm::InputTag("offlineBeamSpot"));
  if (tag.label() != "none") beamSpotToken_ = consumes<reco::BeamSpot>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("PFMetTag", edm::InputTag("pfMet"));
  if (tag.label() != "none") pfMetToken_ = consumes<reco::PFMETCollection>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("PFChMetTag", edm::InputTag("pfChMet"));
  if (tag.label() != "none") pfChMetToken_ = consumes<reco::PFMETCollection>(tag);
 
  tag = cfg.getUntrackedParameter<edm::InputTag>("CaloMetTag", edm::InputTag("caloMet"));
  if (tag.label() != "none") caloMetToken_ = consumes<reco::CaloMETCollection>(tag); 

  tag = cfg.getUntrackedParameter<edm::InputTag>("GenTag", edm::InputTag("prunedGenParticles"));
  if (tag.label() != "none") genToken_ = consumes<reco::GenParticleCollection>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("PileUpInfoTag", edm::InputTag("pileupInfo"));
  if (tag.label() != "none") pileUpInfoToken_ = consumes<std::vector<PileupSummaryInfo> >(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("GenInfoTag", edm::InputTag("generator"));
  if (tag.label() != "none") genInfoToken_ = consumes<GenEventInfoProduct>(tag);  

  tag = cfg.getUntrackedParameter<edm::InputTag>("ScalersTag", edm::InputTag("scalersRawToDigi"));
  if (tag.label() != "none") scalersToken_ = consumes<LumiScalersCollection>(tag);
    
  tag = cfg.getUntrackedParameter<edm::InputTag>("l1MuonsTag", edm::InputTag("gmtStage2Digis:Muon:"));
  if (tag.label() != "none") l1Token_ = consumes<l1t::MuonBxCollection>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("pfCandsTag", edm::InputTag("packedPFCandidates"));
  if (tag.label() != "none") pfCandsToken_ = consumes<std::vector<pat::PackedCandidate>>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("PUPPIMuonIsoTag", edm::InputTag("PUPPIMuonRelIso:PuppiCombined"));
  if (tag.label() != "none") PUPPIToken_ = consumes<edm::ValueMap<double>>(tag);
  tag = cfg.getUntrackedParameter<edm::InputTag>("PUPPIMuonIsoCHTag", edm::InputTag("PUPPIMuonRelIso:PuppiCombinedCH"));
  if (tag.label() != "none") PUPPICHToken_ = consumes<edm::ValueMap<double>>(tag);
  tag = cfg.getUntrackedParameter<edm::InputTag>("PUPPIMuonIsoNHTag", edm::InputTag("PUPPIMuonRelIso:PuppiCombinedNH"));
  if (tag.label() != "none") PUPPINHToken_ = consumes<edm::ValueMap<double>>(tag);
  tag = cfg.getUntrackedParameter<edm::InputTag>("PUPPIMuonIsoPHTag", edm::InputTag("PUPPIMuonRelIso:PuppiCombinedPH"));
  if (tag.label() != "none") PUPPIPHToken_ = consumes<edm::ValueMap<double>>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("PUPPIMuonIsoLepTag", edm::InputTag("PUPPIMuonRelIso:PuppiWithLepton"));
  if (tag.label() != "none") PUPPILepToken_ = consumes<edm::ValueMap<double>>(tag);
  tag = cfg.getUntrackedParameter<edm::InputTag>("PUPPIMuonIsoNoLepTag", edm::InputTag("PUPPIMuonRelIso:PuppiWithoutLepton"));
  if (tag.label() != "none") PUPPINoLepToken_ = consumes<edm::ValueMap<double>>(tag);

  m_minMuPtCut = cfg.getUntrackedParameter<double>("MinMuPtCut", 0.);
  m_minNMuCut  = cfg.getUntrackedParameter<int>("MinNMuCut",  0.);

}


void MuonPogTreeProducer::beginJob() 
{
  
  edm::Service<TFileService> fs;
  tree_["muPogTree"] = fs->make<TTree>("MUONPOGTREE","Muon POG Tree");

  int splitBranches = 2;
  tree_["muPogTree"]->Branch("event",&event_,64000,splitBranches);
  tree_["muPogTree"]->Branch("eventId",&eventId_,64000,splitBranches);

}


void MuonPogTreeProducer::beginRun(const edm::Run & run, const edm::EventSetup & config )
{
  
}


void MuonPogTreeProducer::endJob() 
{

}


void MuonPogTreeProducer::analyze (const edm::Event & ev, const edm::EventSetup &)
{

  // Clearing branch variables
  // and setting default values
  event_.hlt.triggers.clear();
  event_.hlt.objects.clear();
  event_.l1muons.clear();

  event_.genParticles.clear();
  event_.genInfos.clear();
  event_.muons.clear();
  event_.jets.clear();
  event_.genjets.clear();
  
  event_.mets.pfMet   = -999; 
  event_.mets.pfChMet = -999; 
  event_.mets.caloMet = -999; 

  for (unsigned int ix=0; ix<3; ++ix) {
    event_.primaryVertex[ix] = 0.;
    for (unsigned int iy=0; iy<3; ++iy) {
      event_.cov_primaryVertex[ix][iy] = 0.;
    }
  }
  event_.nVtx = -1;


  // Fill general information
  // run, luminosity block, event
  event_.runNumber = ev.id().run();
  event_.luminosityBlockNumber = ev.id().luminosityBlock();
  event_.eventNumber = ev.id().event();

  eventId_.runNumber = ev.id().run();
  eventId_.luminosityBlockNumber = ev.id().luminosityBlock();
  eventId_.eventNumber = ev.id().event();
    
  // Fill GEN pile up information
  if (!ev.isRealData()) 
    {
      if (!pileUpInfoToken_.isUninitialized() &&
	  !genInfoToken_.isUninitialized()) 
	{
	  edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
	  edm::Handle<GenEventInfoProduct> genInfo;

	  if (ev.getByToken(pileUpInfoToken_, puInfo) &&
	      ev.getByToken(genInfoToken_, genInfo) ) 
	    fillGenInfo(puInfo,genInfo);
	  else 
	    edm::LogError("") << "[MuonPogTreeProducer]: Pile-Up Info collection does not exist !!!";
	}      
    }
  

  // Fill GEN particles information
  if (!ev.isRealData()) 
    {
      if (!genToken_.isUninitialized() ) 
	{ 
	  edm::Handle<reco::GenParticleCollection> genParticles;
	  if (ev.getByToken(genToken_, genParticles)) 
	    fillGenParticles(genParticles);
	  else 
	    edm::LogError("") << ">>> GEN collection does not exist !!!";
	}
    }

  if (ev.isRealData()) 
    {

      event_.bxId  = ev.bunchCrossing();
      event_.orbit = ev.orbitNumber();

      if (!scalersToken_.isUninitialized()) 
        { 
          edm::Handle<LumiScalersCollection> lumiScalers;
          if (ev.getByToken(scalersToken_, lumiScalers) && 
              lumiScalers->size() > 0 ) 
            event_.instLumi  = lumiScalers->begin()->instantLumi();
          else 
            edm::LogError("") << ">>> Scaler collection does not exist !!!";
        }
    }

  // Fill trigger information
  if (!trigResultsToken_.isUninitialized() &&
      !trigSummaryToken_.isUninitialized()) 
    {
      
      edm::Handle<edm::TriggerResults> triggerResults;
      edm::Handle<trigger::TriggerEvent> triggerEvent;
      
      if (ev.getByToken(trigResultsToken_, triggerResults) &&
	  ev.getByToken(trigSummaryToken_, triggerEvent)) 
	fillHlt(triggerResults, triggerEvent,ev.triggerNames(*triggerResults));
      else 
	edm::LogError("") << "[MuonPogTreeProducer]: Trigger collections do not exist !!!";
    }
  
  
  // Fill vertex information
  edm::Handle<std::vector<reco::Vertex> > vertexes;

  if(!primaryVertexToken_.isUninitialized()) 
    {
      if (ev.getByToken(primaryVertexToken_, vertexes))
	fillPV(vertexes);
      else 
	edm::LogError("") << "[MuonPogTreeProducer]: Vertex collection does not exist !!!";
    }

  // Get beam spot for muons
  edm::Handle<reco::BeamSpot> beamSpot;
  if (!beamSpotToken_.isUninitialized() ) 
    { 
      if (!ev.getByToken(beamSpotToken_, beamSpot)) 
	edm::LogError("") << "[MuonPogTreeProducer]: Beam spot collection not found !!!";
    }

  // Fill (raw) MET information: PF, PF charged, Calo    
  edm::Handle<reco::PFMETCollection> pfMet; 
  if(!pfMetToken_.isUninitialized()) 
    { 
      if (!ev.getByToken(pfMetToken_, pfMet)) 
	edm::LogError("") << "[MuonPogTreeProducer] PFMet collection does not exist !!!"; 
      else { 
	const reco::PFMET &iPfMet = (*pfMet)[0]; 
	event_.mets.pfMet = iPfMet.et(); 
      } 
    } 

  edm::Handle<reco::PFMETCollection> pfChMet; 
  if(!pfChMetToken_.isUninitialized()) 
    { 
      if (!ev.getByToken(pfChMetToken_, pfChMet)) 
	edm::LogError("") << "[MuonPogTreeProducer] PFChMet collection does not exist !!!"; 
      else { 
	const reco::PFMET &iPfChMet = (*pfChMet)[0]; 
	event_.mets.pfChMet = iPfChMet.et(); 
      } 
    } 

  edm::Handle<reco::CaloMETCollection> caloMet; 
  if(!caloMetToken_.isUninitialized()) 
    { 
      if (!ev.getByToken(caloMetToken_, caloMet)) 
	edm::LogError("") << "[MuonPogTreeProducer] CaloMet collection does not exist !!!"; 
      else { 
	const reco::CaloMET &iCaloMet = (*caloMet)[0]; 
	event_.mets.caloMet = iCaloMet.et(); 
      } 
    } 

  // std::cout << "-----------------------------------------------------------" << std::endl;
  // std::cout << "-----------------------------------------------------------" << std::endl;
  // Get muons  
  edm::Handle<edm::View<pat::Muon> > muons;
  if (!muonToken_.isUninitialized() ) 
    { 
      if (!ev.getByToken(muonToken_, muons)) 
	edm::LogError("") << "[MuonPogTreeProducer] Muon collection does not exist !!!";
    }

  // Get PF Candidates  
  edm::Handle< std::vector<pat::PackedCandidate> > pfCands;
  if (!pfCandsToken_.isUninitialized() ) 
    { 
      if (!ev.getByToken(pfCandsToken_, pfCands)) 
	edm::LogError("") << "[MuonPogTreeProducer] PF Candidates collection does not exist !!!";
    }

  // Get PUPPI Isolation  
  edm::Handle< edm::ValueMap<double> > PUPPICom, PUPPICHCom,PUPPINHCom,PUPPIPHCom;
  edm::Handle< edm::ValueMap<double> > PUPPILep;
  edm::Handle< edm::ValueMap<double> > PUPPINoLep;
  if (!PUPPIToken_.isUninitialized() || !PUPPICHToken_.isUninitialized() || !PUPPINHToken_.isUninitialized() || !PUPPIPHToken_.isUninitialized() || 
      !PUPPILepToken_.isUninitialized() || !PUPPINoLepToken_.isUninitialized() ) 
    { 
      if (!ev.getByToken(PUPPIToken_,PUPPICom) || !ev.getByToken(PUPPICHToken_,PUPPICHCom) || !ev.getByToken(PUPPINHToken_,PUPPINHCom) || !ev.getByToken(PUPPIPHToken_,PUPPIPHCom) || 
	  !ev.getByToken(PUPPILepToken_,PUPPILep) || !ev.getByToken(PUPPINoLepToken_,PUPPINoLep)) 
	edm::LogError("") << "[MuonPogTreeProducer] PUPPI isolation collection does not exist !!!";
    }
  

  Int_t nGoodMuons = 0;
  eventId_.maxPTs.clear();
  // Fill muon information
  if (muons.isValid() && vertexes.isValid() && beamSpot.isValid() && pfCands.isValid() &&
      PUPPICom.isValid() && PUPPICHCom.isValid() && PUPPINHCom.isValid() && PUPPIPHCom.isValid() && 
      PUPPILep.isValid() && PUPPINoLep.isValid()) 
    {
      nGoodMuons = fillMuons(muons,vertexes,beamSpot,pfCands,
			     PUPPICom,PUPPICHCom,PUPPINHCom,PUPPIPHCom,
			     PUPPILep,PUPPINoLep);
    }
  eventId_.nMuons = nGoodMuons;


  // Get genjets  
  edm::Handle<edm::View<reco::GenJet> > genjets;
  if (!genjetsToken_.isUninitialized() ) 
    { 
      if (!ev.getByToken(genjetsToken_, genjets)) 
	edm::LogError("") << "[MuonPogTreeProducer] GenJet collection does not exist !!!";
    }
  
  // Fill jet information if there are good muons in the event
  if (genjets.isValid() && nGoodMuons >= m_minNMuCut) 
    {
      fillGenJets(genjets);
    }
    
  // Get jets  
  edm::Handle<edm::View<pat::Jet> > jets;
  if (!jetsToken_.isUninitialized() ) 
    { 
      if (!ev.getByToken(jetsToken_, jets)) 
  	edm::LogError("") << "[MuonPogTreeProducer] PFJet collection does not exist !!!";
    }
  
  // Fill jet information if there are good muons in the event
  if (jets.isValid() && nGoodMuons >= m_minNMuCut) 
    {
      // std::cout << "Filling Jets" << std::endl;
      fillJets(jets);
    }
    
  //Fill L1 informations
  edm::Handle<l1t::MuonBxCollection> l1s;
  if (!l1Token_.isUninitialized() )
    {
        if (!ev.getByToken(l1Token_, l1s))
	  edm::LogError("") << "[MuonPogTreeProducer] L1 muon bx collection does not exist !!!";
        else {
	  //fillL1(l1s);
        }
    }
    
  if (nGoodMuons >= m_minNMuCut)
  tree_["muPogTree"]->Fill();

}

void MuonPogTreeProducer::fillGenInfo(const edm::Handle<std::vector<PileupSummaryInfo> > & puInfo,
				      const edm::Handle<GenEventInfoProduct> & gen)
{

  muon_pog::GenInfo genInfo;
  
  genInfo.trueNumberOfInteractions     = -1.;
  genInfo.actualNumberOfInteractions   = -1.;
  genInfo.genWeight = gen->weight() ;
	
  std::vector<PileupSummaryInfo>::const_iterator puInfoIt  = puInfo->begin();
  std::vector<PileupSummaryInfo>::const_iterator puInfoEnd = puInfo->end();

  for(; puInfoIt != puInfoEnd; ++puInfoIt) 
    {
      int bx = puInfoIt->getBunchCrossing();
	  
      if(bx == 0) 
	{ 
	  genInfo.trueNumberOfInteractions   = puInfoIt->getTrueNumInteractions();
	  genInfo.actualNumberOfInteractions = puInfoIt->getPU_NumInteractions();
	  continue;
	}
    }
  
  event_.genInfos.push_back(genInfo);
  
}


void MuonPogTreeProducer::fillGenParticles(const edm::Handle<reco::GenParticleCollection> & genParticles)
{
  
  unsigned int gensize = genParticles->size();
  
  // Do not record the initial protons
  for (unsigned int i=0; i<gensize; ++i) 
    {

      const reco::GenParticle& part = genParticles->at(i);
    
      muon_pog::GenParticle gensel;
      gensel.pdgId = part.pdgId();
      gensel.status = part.status();
      gensel.energy = part.energy();
      gensel.pt = part.pt();
      gensel.eta = part.eta();
      gensel.phi = part.phi();
      gensel.vx = part.vx();
      gensel.vy = part.vy();
      gensel.vz = part.vz();

      // Full set of GenFlags
      gensel.flags.clear();
      reco::GenStatusFlags statusflags = part.statusFlags();
      if (statusflags.flags_.size() == 15)
	for (unsigned int flag = 0; flag < statusflags.flags_.size(); ++flag)
	  gensel.flags.push_back(statusflags.flags_[flag]);      
      
      gensel.mothers.clear();
      unsigned int nMothers = part.numberOfMothers();

      for (unsigned int iMother=0; iMother<nMothers; ++iMother) 
	{
	  gensel.mothers.push_back(part.motherRef(iMother)->pdgId());
	}

      // Protect agains bug in genParticles (missing mother => first proton)
      if (i>=2 && nMothers==0) gensel.mothers.push_back(0);
      
      event_.genParticles.push_back(gensel);
    }
  
}


void MuonPogTreeProducer::fillHlt(const edm::Handle<edm::TriggerResults> & triggerResults, 
				  const edm::Handle<trigger::TriggerEvent> & triggerEvent,
				  const edm::TriggerNames & triggerNames)
{    

  for (unsigned int iTrig=0; iTrig<triggerNames.size(); ++iTrig) 
    {
      
      if (triggerResults->accept(iTrig)) 
	{
	  std::string pathName = triggerNames.triggerName(iTrig);
	  if (trigPathCut_ == "all" || pathName.find(trigPathCut_) != std::string::npos)
	    event_.hlt.triggers.push_back(pathName);
	}
    }
      
  const trigger::size_type nFilters(triggerEvent->sizeFilters());

  for (trigger::size_type iFilter=0; iFilter!=nFilters; ++iFilter) 
    {
	
      std::string filterTag = triggerEvent->filterTag(iFilter).encode();
      
      if (trigFilterCut_ == "all" || filterTag.find(trigFilterCut_) != std::string::npos)
	{

	  trigger::Keys objectKeys = triggerEvent->filterKeys(iFilter);
	  const trigger::TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());
	
	  for (trigger::size_type iKey=0; iKey<objectKeys.size(); ++iKey) 
	    {  
	      trigger::size_type objKey = objectKeys.at(iKey);
	      const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);
	      
	      muon_pog::HLTObject hltObj;
	      
	      float trigObjPt = triggerObj.pt();
	      float trigObjEta = triggerObj.eta();
	      float trigObjPhi = triggerObj.phi();
	      
	      hltObj.filterTag = filterTag;
	      
	      hltObj.pt  = trigObjPt;
	      hltObj.eta = trigObjEta;
	      hltObj.phi = trigObjPhi;
	      
	      event_.hlt.objects.push_back(hltObj);
	      
	    }
	}
    }

}

void MuonPogTreeProducer::fillL1(const edm::Handle<l1t::MuonBxCollection> & l1MuonBxColl)
{

  for (int ibx = l1MuonBxColl->getFirstBX(); ibx <= l1MuonBxColl->getLastBX(); ++ibx) 
    {
      for (auto l1MuIt = l1MuonBxColl->begin(ibx); l1MuIt != l1MuonBxColl->end(ibx); ++l1MuIt)
	{

	  muon_pog::L1Muon l1part;
	  l1part.pt = l1MuIt->pt();
	  l1part.eta = l1MuIt->eta();
	  l1part.phi = l1MuIt->phi();
	  l1part.charge = l1MuIt->hwChargeValid() ? l1MuIt->charge() : 0;
	  
	  l1part.quality = l1MuIt->hwQual();
	  l1part.bx = ibx;
	  
	  l1part.tfIndex = l1MuIt->tfMuonIndex();
	  
	  event_.l1muons.push_back(l1part);
	  
	}
    }
}



void MuonPogTreeProducer::fillPV(const edm::Handle<std::vector<reco::Vertex> > & vertexes)
{
      
  int nVtx = 0;

  std::vector<reco::Vertex>::const_iterator vertexIt  = vertexes->begin();
  std::vector<reco::Vertex>::const_iterator vertexEnd = vertexes->end();

  for (; vertexIt != vertexEnd; ++vertexIt) 
    {

      const reco::Vertex& vertex = *vertexIt;

      if (!vertex.isValid()) continue;
      ++nVtx;

      if (vertexIt == vertexes->begin()) 
	{
	  event_.primaryVertex[0] = vertex.x();
	  event_.primaryVertex[1] = vertex.y();
	  event_.primaryVertex[2] = vertex.z();

	  for (unsigned int ix=0; ix<3; ++ix) 
	    {
	      for (unsigned int iy=0; iy<3; ++iy) 
		{
		  event_.cov_primaryVertex[ix][iy] = vertex.covariance(ix,iy);
		}
	    }
	}
    }
  
  event_.nVtx = nVtx;
  
}


Int_t MuonPogTreeProducer::fillMuons(const edm::Handle<edm::View<pat::Muon> > & muons,
				     const edm::Handle<std::vector<reco::Vertex> > & vertexes,
				     const edm::Handle<reco::BeamSpot> & beamSpot,
				     const edm::Handle<std::vector<pat::PackedCandidate>> & pfCands,
				     const edm::Handle<edm::ValueMap<double>>& puppiCom,
				     const edm::Handle<edm::ValueMap<double>>& puppiCHCom,
				     const edm::Handle<edm::ValueMap<double>>& puppiNHCom,
				     const edm::Handle<edm::ValueMap<double>>& puppiPHCom,
				     const edm::Handle<edm::ValueMap<double>>& puppiLep,
				     const edm::Handle<edm::ValueMap<double>>& puppiNoLep)
{
  // edm::View<reco::Muon>::const_iterator muonIt  = muons->begin();
  // edm::View<reco::Muon>::const_iterator muonEnd = muons->end();
  
  // for (; muonIt != muonEnd; ++muonIt) 
  //   {
      
  for (unsigned int im=0; im<muons->size(); im++)
    //const reco::Muon& mu = (*muonIt);
    {
      const reco::Muon& mu = muons->at(im);

      bool isGlobal      = mu.isGlobalMuon();
      bool isTracker     = mu.isTrackerMuon();
      bool isTrackerArb  = muon::isGoodMuon(mu, muon::TrackerMuonArbitrated); 
      bool isRPC         = mu.isRPCMuon();
      bool isStandAlone  = mu.isStandAloneMuon();
      bool isPF          = mu.isPFMuon();

      bool hasInnerTrack = !mu.innerTrack().isNull();
      
      muon_pog::Muon ntupleMu;
      
      ntupleMu.pt     = mu.pt();
      ntupleMu.eta    = mu.eta();
      ntupleMu.phi    = mu.phi();
      ntupleMu.charge = mu.charge();

      ntupleMu.isSoft    = 0;	  
      ntupleMu.isTight   = 0;	  
      ntupleMu.isHighPt  = 0;
      ntupleMu.isLoose   = muon::isLooseMuon(mu)  ? 1 : 0;	  
      ntupleMu.isMedium  = muon::isMediumMuon(mu) ? 1 : 0;	  

      // Detector Based Isolation
      reco::MuonIsolation detIso03 = mu.isolationR03();

      ntupleMu.trackerIso = detIso03.sumPt;
      ntupleMu.EMCalIso   = detIso03.emEt;
      ntupleMu.HCalIso    = detIso03.hadEt;
      
      // PF Isolation
      reco::MuonPFIsolation pfIso04 = mu.pfIsolationR04();
      reco::MuonPFIsolation pfIso03 = mu.pfIsolationR03();

      ntupleMu.chargedHadronIso   = pfIso04.sumChargedHadronPt;
      ntupleMu.chargedHadronIsoPU = pfIso04.sumPUPt; 
      ntupleMu.neutralHadronIso   = pfIso04.sumNeutralHadronEt;
      ntupleMu.photonIso          = pfIso04.sumPhotonEt;

      ntupleMu.isoPflow04 = (pfIso04.sumChargedHadronPt+ 
       			     std::max(0.,pfIso04.sumPhotonEt+pfIso04.sumNeutralHadronEt - 0.5*pfIso04.sumPUPt)) / mu.pt();
    
      ntupleMu.isoPflow03 = (pfIso03.sumChargedHadronPt+ 
       			     std::max(0.,pfIso03.sumPhotonEt+pfIso03.sumNeutralHadronEt - 0.5*pfIso03.sumPUPt)) / mu.pt();

      // NEW Selectors
      ntupleMu.Sel_CutBasedIdLoose        = mu.passed(reco::Muon::CutBasedIdLoose);
      ntupleMu.Sel_CutBasedIdMedium       = mu.passed(reco::Muon::CutBasedIdMedium);
      ntupleMu.Sel_CutBasedIdMediumPrompt = mu.passed(reco::Muon::CutBasedIdMediumPrompt);
      ntupleMu.Sel_CutBasedIdTight        = mu.passed(reco::Muon::CutBasedIdTight);
      ntupleMu.Sel_CutBasedIdGlobalHighPt = mu.passed(reco::Muon::CutBasedIdGlobalHighPt);
      ntupleMu.Sel_CutBasedIdTrkHighPt    = mu.passed(reco::Muon::CutBasedIdTrkHighPt);
      ntupleMu.Sel_SoftCutBasedId         = mu.passed(reco::Muon::SoftCutBasedId);
      ntupleMu.Sel_SoftMvaId              = mu.passed(reco::Muon::SoftMvaId);
      ntupleMu.Sel_MvaLoose               = mu.passed(reco::Muon::MvaLoose);
      ntupleMu.Sel_MvaMedium              = mu.passed(reco::Muon::MvaMedium);
      ntupleMu.Sel_MvaTight               = mu.passed(reco::Muon::MvaTight);

      // NEW Isolations
      // PUPPI
      edm::Ptr<reco::Muon> muPtr = muons->ptrAt( im );
      ntupleMu.PUPPIIso   = (*puppiCom)[muPtr];
      ntupleMu.PUPPIIsoCH = (*puppiCHCom)[muPtr];
      ntupleMu.PUPPIIsoNH = (*puppiNHCom)[muPtr];
      ntupleMu.PUPPIIsoPH = (*puppiPHCom)[muPtr];

      ntupleMu.PUPPILepIso   = (*puppiLep)[muPtr];
      ntupleMu.PUPPINoLepIso = (*puppiNoLep)[muPtr];
      // -- MiniIsolation
      const pat::Muon& pmu = muons->at(im);
      ntupleMu.MiniIsoCH = pmu.miniPFIsolation().chargedHadronIso();
      ntupleMu.MiniIsoNH = pmu.miniPFIsolation().neutralHadronIso();
      ntupleMu.MiniIsoPH = pmu.miniPFIsolation().photonIso();
      ntupleMu.MiniIsoPU = pmu.miniPFIsolation().puChargedHadronIso();
      
      // New Gen (SimHit) Matching
      ntupleMu.IsMatchedPrimaryMuon = pmu.simType() == reco::MuonSimType::MatchedPrimaryMuon;
      ntupleMu.IsNotMatched         = pmu.simType() == reco::MuonSimType::NotMatched;
      ntupleMu.IsMatchedElectron    = pmu.simType() == reco::MuonSimType::MatchedElectron;
      ntupleMu.IsMatchedMuonHF      = pmu.simType() == reco::MuonSimType::MatchedMuonFromHeavyFlavour;
      ntupleMu.IsMatchedMuonLF      = pmu.simType() == reco::MuonSimType::MatchedMuonFromLightFlavour;
      ntupleMu.IsGhostPrimaryMuon   = pmu.simType() == reco::MuonSimType::GhostPrimaryMuon;
      ntupleMu.IsGhostMuonHF        = pmu.simType() == reco::MuonSimType::GhostMuonFromHeavyFlavour;
      ntupleMu.IsGhostMuonLF        = pmu.simType() == reco::MuonSimType::GhostMuonFromLightFlavour;
      // Additional Info
      ntupleMu.SimpdgId             = pmu.simPdgId();
      ntupleMu.SimmotherPdgId       = pmu.simMotherPdgId();

      ntupleMu.isGlobal     = isGlobal ? 1 : 0;	
      ntupleMu.isTracker    = isTracker ? 1 : 0;	
      ntupleMu.isTrackerArb = isTrackerArb ? 1 : 0;	
      ntupleMu.isRPC        = isRPC ? 1 : 0;
      ntupleMu.isStandAlone = isStandAlone ? 1 : 0;
      ntupleMu.isPF         = isPF ? 1 : 0;

      ntupleMu.glbNormChi2              = isGlobal      ? mu.globalTrack()->normalizedChi2() : -999; 
      ntupleMu.trkNormChi2	        = hasInnerTrack ? mu.innerTrack()->normalizedChi2()  : -999; 
      ntupleMu.trkMuonMatchedStations   = isTracker     ? mu.numberOfMatchedStations()       : -999; 
      ntupleMu.glbMuonValidHits	        = isGlobal      ? mu.globalTrack()->hitPattern().numberOfValidMuonHits()       : -999; 
      ntupleMu.trkPixelValidHits	= hasInnerTrack ? mu.innerTrack()->hitPattern().numberOfValidPixelHits()       : -999; 
      ntupleMu.trkPixelLayersWithMeas   = hasInnerTrack ? mu.innerTrack()->hitPattern().pixelLayersWithMeasurement()   : -999; 
      ntupleMu.trkTrackerLayersWithMeas = hasInnerTrack ? mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() : -999; 

      ntupleMu.bestMuPtErr              = mu.muonBestTrack()->ptError(); 

      ntupleMu.trkValidHitFrac = hasInnerTrack           ? mu.innerTrack()->validFraction()       : -999; 
      ntupleMu.trkStaChi2      = isGlobal                ? mu.combinedQuality().chi2LocalPosition : -999; 
      ntupleMu.trkKink         = isGlobal                ? mu.combinedQuality().trkKink           : -999; 
      ntupleMu.muSegmComp      = (isGlobal || isTracker) ? muon::segmentCompatibility(mu)         : -999; 

      ntupleMu.isTrkMuOST               = muon::isGoodMuon(mu, muon::TMOneStationTight) ? 1 : 0; 
      ntupleMu.isTrkHP                  = hasInnerTrack && mu.innerTrack()->quality(reco::TrackBase::highPurity) ? 1 : 0; 

      if ( mu.isMatchesValid() && ntupleMu.isTrackerArb )
      	{
      	  for ( reco::MuonChamberMatch match : mu.matches() )
      	    {
      	      muon_pog::ChambMatch ntupleMatch;
	      
      	      if ( getMuonChamberId(match.id,
      				    ntupleMatch.type,ntupleMatch.r,
      				    ntupleMatch.phi,ntupleMatch.eta)
      		   )
      		{
	      
      		  ntupleMatch.x = mu.trackX(match.station(),match.detector());
      		  ntupleMatch.y = mu.trackY(match.station(),match.detector());
      		  ntupleMatch.dXdZ = mu.trackDxDz(match.station(),match.detector());
      		  ntupleMatch.dYdZ = mu.trackDyDz(match.station(),match.detector());

      		  ntupleMatch.errxTk = mu.trackXErr(match.station(),match.detector());
      		  ntupleMatch.erryTk = mu.trackYErr(match.station(),match.detector());

      		  ntupleMatch.errDxDzTk = mu.trackDxDzErr(match.station(),match.detector());
      		  ntupleMatch.errDyDzTk = mu.trackDyDzErr(match.station(),match.detector());
	      
      		  ntupleMatch.dx = mu.dX(match.station(),match.detector());
      		  ntupleMatch.dy = mu.dY(match.station(),match.detector());
      		  ntupleMatch.dDxDz = mu.dDxDz(match.station(),match.detector());
      		  ntupleMatch.dDyDz = mu.dDxDz(match.station(),match.detector());
		  
      		  ntupleMatch.errxSeg = mu.segmentXErr(match.station(),match.detector());
      		  ntupleMatch.errySeg = mu.segmentYErr(match.station(),match.detector());
      		  ntupleMatch.errDxDzSeg = mu.segmentDxDzErr(match.station(),match.detector());
      		  ntupleMatch.errDyDzSeg = mu.segmentDyDzErr(match.station(),match.detector());
		  
      		  ntupleMu.matches.push_back(ntupleMatch);
      		}
      	    }
      	}
      
      ntupleMu.dxyBest  = -999; 
      ntupleMu.dzBest   = -999; 
      ntupleMu.dxyInner = -999; 
      ntupleMu.dzInner  = -999; 

      double dxybs = isGlobal ? mu.globalTrack()->dxy(beamSpot->position()) :
      	hasInnerTrack ? mu.innerTrack()->dxy(beamSpot->position()) : -1000;
      double dzbs  = isGlobal ? mu.globalTrack()->dz(beamSpot->position()) :
      	hasInnerTrack ? mu.innerTrack()->dz(beamSpot->position()) : -1000;

      double dxy = -1000.;
      double dz  = -1000.;

      if (vertexes->size() > 0)
      	{
      	  const reco::Vertex & vertex = vertexes->at(0);

      	  dxy = isGlobal ? mu.globalTrack()->dxy(vertex.position()) :
      	    hasInnerTrack ? mu.innerTrack()->dxy(vertex.position()) : -1000;
      	  dz = isGlobal ? mu.globalTrack()->dz(vertex.position()) :
      	    hasInnerTrack ? mu.innerTrack()->dz(vertex.position()) : -1000;
 
      	  ntupleMu.dxyBest  = mu.muonBestTrack()->dxy(vertex.position()); 
      	  ntupleMu.dzBest   = mu.muonBestTrack()->dz(vertex.position()); 
      	  if(hasInnerTrack) { 
      	    ntupleMu.dxyInner = mu.innerTrack()->dxy(vertex.position()); 
      	    ntupleMu.dzInner  = mu.innerTrack()->dz(vertex.position()); 
      	  } 

      	  ntupleMu.isSoft    = muon::isSoftMuon(mu,vertex)   ? 1 : 0;	  
      	  ntupleMu.isTight   = muon::isTightMuon(mu,vertex)  ? 1 : 0;	  
      	  ntupleMu.isHighPt  = muon::isHighPtMuon(mu,vertex) ? 1 : 0;

      	}

      ntupleMu.dxy    = dxy;
      ntupleMu.dz     = dz;
      ntupleMu.edxy   = isGlobal ? mu.globalTrack()->dxyError() : hasInnerTrack ? mu.innerTrack()->dxyError() : -1000;
      ntupleMu.edz    = isGlobal ? mu.globalTrack()->dzError()  : hasInnerTrack ? mu.innerTrack()->dzError() : -1000;

      ntupleMu.dxybs  = dxybs;
      ntupleMu.dzbs   = dzbs;

      if(mu.isTimeValid()) { 
      	ntupleMu.muonTimeDof = mu.time().nDof; 
      	ntupleMu.muonTime    = mu.time().timeAtIpInOut; 
      	ntupleMu.muonTimeErr = mu.time().timeAtIpInOutErr; 
      } 
      else { 
      	ntupleMu.muonTimeDof = -999; 
      	ntupleMu.muonTime    = -999; 
      	ntupleMu.muonTimeErr = -999; 
      } 

      if(mu.rpcTime().nDof > 0) { 
      	ntupleMu.muonRpcTimeDof = mu.rpcTime().nDof; 
      	ntupleMu.muonRpcTime    = mu.rpcTime().timeAtIpInOut; 
      	ntupleMu.muonRpcTimeErr = mu.rpcTime().timeAtIpInOutErr; 
      } 
      else { 
      	ntupleMu.muonRpcTimeDof = -999; 
      	ntupleMu.muonRpcTime    = -999; 
      	ntupleMu.muonRpcTimeErr = -999; 
      } 

      if ( m_minMuPtCut < 0 ||
	   ( (isTracker || isGlobal || isStandAlone) && mu.pt() > m_minMuPtCut )
	   )
	{
	  event_.muons.push_back(ntupleMu);
        
	}

    }

  return event_.muons.size();

}

void MuonPogTreeProducer::fillJets(const edm::Handle<edm::View<pat::Jet> > & jets)
{
  
  // std::cout << "Number of Jets: " << jets->size() << std::endl;
  // std::cout << "Number of Muons: " << event_.muons.size() << std::endl;

  edm::View<pat::Jet>::const_iterator jetIt  = jets->begin();
  edm::View<pat::Jet>::const_iterator jetEnd = jets->end();
  
  for (; jetIt != jetEnd; ++jetIt) 
    {
      const pat::Jet& pfjet = (*jetIt);
      muon_pog::PF04Jets ntupleJet;
      
      if (pfjet.pt()<20) continue;

      ntupleJet.pt  = pfjet.pt();
      ntupleJet.eta = pfjet.eta();
      ntupleJet.phi = pfjet.phi();
      ntupleJet.Ent = pfjet.et();
      ntupleJet.En  = pfjet.energy();

      float NHF    = pfjet.neutralHadronEnergyFraction();
      float NEMF   = pfjet.neutralEmEnergyFraction();
      float CHF    = pfjet.chargedHadronEnergyFraction();
      float MUF    = pfjet.muonEnergyFraction();
      float CEMF   = pfjet.chargedEmEnergyFraction();
      int NumConst = pfjet.chargedMultiplicity()+pfjet.neutralMultiplicity();
      int NumNeutralParticle =pfjet.neutralMultiplicity();
      int CHM      = pfjet.chargedMultiplicity();

      // Definition from 
      bool looseJetID = false;
      looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && 
	( (std::abs(ntupleJet.eta)<=2.4 && CHF>0.0 && CHM>0.0 && CEMF<0.99) || std::abs(ntupleJet.eta)>2.4 ) 
	&& std::abs(ntupleJet.eta)<=2.7;
      
      bool tightJetID = false;
      tightJetID  = (NHF<0.90 && NEMF<0.90 && NumConst>1) && 
	((std::abs(ntupleJet.eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || std::abs(ntupleJet.eta)>2.4) 
	&& std::abs(ntupleJet.eta)<=2.7;
      
      if ( std::abs(ntupleJet.eta) > 3.0 ) {
      	looseJetID = (NEMF<0.90 && NumNeutralParticle>10 && std::abs(ntupleJet.eta)>3.0 );
      	tightJetID = (NEMF<0.90 && NumNeutralParticle>10 && std::abs(ntupleJet.eta)>3.0 );
      }
      else if (std::abs(ntupleJet.eta) > 2.7){
      	looseJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2 && std::abs(ntupleJet.eta)>2.7 && std::abs(ntupleJet.eta)<=3.0 );
      	tightJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2 && std::abs(ntupleJet.eta)>2.7 && std::abs(ntupleJet.eta)<=3.0 );
      }      

      ntupleJet.isTightID = tightJetID;
      
      ntupleJet.MuonEn = MUF;

      // Muon cleaning  
      bool cleanJet = true; 
      for (unsigned int imuon = 0; imuon < event_.muons.size(); imuon++){
	muon_pog::Muon SelMuon = event_.muons.at(imuon);
	if (cleanJet) cleanJet = deltaR(SelMuon.eta, SelMuon.phi, ntupleJet.eta, ntupleJet.phi) > 0.4;
      }
      
      if(looseJetID && cleanJet){
	event_.jets.push_back(ntupleJet);
	//std::cout << "pt = " << ntupleJet.pt << " ; eta = " << ntupleJet.eta << " ; phi = " << ntupleJet.phi << " ; E = " << ntupleJet.En << " ; ntIsTight = " << ntupleJet.isTightID << std::endl;
      }
      
    }

}


void MuonPogTreeProducer::fillGenJets(const edm::Handle<edm::View<reco::GenJet> > & genjets)
{
  
  // std::cout << "Number of Gen-Jets: " << genjets->size() << std::endl;

  edm::View<reco::GenJet>::const_iterator genjetIt  = genjets->begin();
  edm::View<reco::GenJet>::const_iterator genjetEnd = genjets->end();
  
  for (; genjetIt != genjetEnd; ++genjetIt) 
    {
      const reco::GenJet& genjet = (*genjetIt);
      muon_pog::Gen04Jets ntupleGenJet;
      
      ntupleGenJet.pt  = genjet.pt();
      ntupleGenJet.eta = genjet.eta();
      ntupleGenJet.phi = genjet.phi();
      ntupleGenJet.Ent = genjet.et();
      ntupleGenJet.En  = genjet.energy();

      // Muon cleaning  
      bool cleanGenJet = true; 
      for (unsigned int imuon = 0; imuon < event_.muons.size(); imuon++){
	muon_pog::Muon SelMuon = event_.muons.at(imuon);
	if (cleanGenJet) cleanGenJet = deltaR(SelMuon.eta, SelMuon.phi, ntupleGenJet.eta, ntupleGenJet.phi) > 0.4;
      }

      if (ntupleGenJet.pt>20 && cleanGenJet){
	event_.genjets.push_back(ntupleGenJet);
	// std::cout << "Genpt = " << ntupleGenJet.pt << " ; Geneta = " << ntupleGenJet.eta << " ; Genphi = " << ntupleGenJet.phi  << " ; GenE = " << ntupleGenJet.En << std::endl;
      }
    }

}

bool MuonPogTreeProducer::getMuonChamberId(DetId & id, muon_pog::MuonDetType & det,
					   Int_t & r, Int_t & phi, Int_t & eta) const
{

  if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::DT)
    {
      DTChamberId dtId(id.rawId());  
  
      det = muon_pog::MuonDetType::DT;
      r   = dtId.station();
      phi = dtId.sector();
      eta = dtId.wheel();

      return true;
    }

  if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::CSC)
    {
      CSCDetId cscId(id.rawId());
    
      det = muon_pog::MuonDetType::CSC;
      r   = cscId.station() * cscId.zendcap();
      phi = cscId.chamber();
      eta = cscId.ring();

      return true;
    }

  return false;
      
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(MuonPogTreeProducer);
