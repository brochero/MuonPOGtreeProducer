#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include "TLegend.h"
#include "TRegexp.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"

#include "../src/MuonPogTree.h"
#include "../src/Utils.h"
#include "tdrstyle.C"

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream> 
#include <vector>
#include <regex>
#include <map>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/exceptions.hpp>
#include <boost/lexical_cast.hpp>

#include <TError.h>

// Helper classes defintion *****
// 1. SampleConfig : configuration class containing sample information
// 2. TagAndProbeConfig : configuration class containing TnP cuts information
// 3. Plotter : class containing the plot definition and defining the plot filling 
//              for a given sample <= CB modify this to add new variables
// ******************************


namespace muon_pog {
 
  class SampleConfig {

  public :

    // config parameters (public for direct access)

    TString fileName;  
    TString sampleName;  
    Float_t QCDWeight;
    Float_t nEvents;
    Bool_t applyReweighting;
    std::vector<int> runs;
        
    SampleConfig() {};
    
#ifndef __MAKECINT__ // CB CINT doesn't like boost :'-(    
    SampleConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~SampleConfig() {};

  private:
    std::vector<int> toArray(const std::string & entries); 
    
  };

  class TagAndProbeConfig {

  public :
    
    // config parameters (public for direct access)
    
    Float_t     gen_DrCut;
    
    std::string muon_trackType; // applies to both, tag and probe     
  
    Float_t probe_minPt;
    Float_t probe_isoLooseWP, probe_isoTightWP;      
   
    TagAndProbeConfig() {};
    
#ifndef __MAKECINT__ // CB CINT doesn't like boost :'-(    
    TagAndProbeConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~TagAndProbeConfig() {};
    
  private:
    std::vector<TString> toArray(const std::string & entries); 
  
  };


  class Plotter {

  public :
    
    TString IDName[7] = {"GLB",  "TRK", "LOOSE", "MEDIUM", "TIGHT", "SOFT", "HIGHpT"};
    enum BinType   {GLB=1, TRK, LOOSE, TIGHT, MEDIUM, SOFT, HIGHpT};
    TString region[4] = {"Full", "Barrel", "Endcap", "Overlap"};
    enum BinRegion {Full=1, Barrel, Endcap, Overlap};
    TString isowp[3]  = {"ISONoCut", "ISOLoose", "ISOTight"};  
    enum BinIsoWP  {ISONoCut=1, IsoLoose, IsoTight};
    TString PUr [4] = {"FullPU", "LowPU", "MediumPU", "HighPU"};  
    enum BinPUr{FullPU=1, LowPU, MediumPU, HighPU};

    
    Plotter(muon_pog::TagAndProbeConfig tnpConfig, muon_pog::SampleConfig & sampleConfig) :
      m_tnpConfig(tnpConfig) , m_sampleConfig(sampleConfig) {};
    ~Plotter() {};
    
    void book(TFile *outFile);
    void fill(const std::vector<muon_pog::Muon> & muons, const muon_pog::HLT & hlt, const muon_pog::Event & ev, float weight);
    void fillGen(const std::vector<muon_pog::GenParticle> & genpars, const muon_pog::Event & ev);

    std::map<TString, TH1D *>        m_plots;
    std::map<TString, TH2D *>        m_2Dplots;
    std::map<TString, TProfile *>    m_prof;
    std::map<TString, TEfficiency *> m_effs;

    std::map<TString, TH2D *> m_2Dyields;

    TagAndProbeConfig m_tnpConfig;
    SampleConfig      m_sampleConfig;
        
  };
  
}

// Helper classes defintion *****
// 1. parseConfig : parse the full cfg file
// 1. comparisonPlot : make a plot overlayng data and MC for a given plot
// ******************************

namespace muon_pog {
  void parseConfig(const std::string configFile, TagAndProbeConfig & tpConfig,
		   std::vector<SampleConfig> & sampleConfigs);
  
  void comparisonPlots(std::vector<Plotter> & plotters,
		       TFile *outFile, TString &  outputDir);

  void print_progress(int TreeEntries, Long64_t ievt);
  void copyPhp(const TString &  outputDir);
  Bool_t IsMediumHIP(const muon_pog::Muon & muon);
  Bool_t IsSoftHIP  (const muon_pog::Muon & muon);
  // void setTProfY(TProfile &prof1, TProfile &prof2);
}



// The main program******** *****
// 1. Get configuration file and produces configurations
// 2. Create Plotters and loop on the event to fill them
// 3. Writes results in cnfigurable outuput file
// ******************************

int main(int argc, char* argv[]){
  using namespace muon_pog;

  gErrorIgnoreLevel=kError;


  if (argc != 4) 
    {
      std::cout << "Usage : "
		<< argv[0] << " PATH_TO_CONFIG_FILE PATH_TO_OUTPUT_DIR OUTPUT_FILE_NAME\n";
      exit(100);
    }

  std::string configFile(argv[1]);
  
  std::cout << "[" << argv[0] << "] Using config file " << configFile << std::endl;

  // Output directory
  TString dirName    = argv[2];
  TString outputName = argv[3];
  system("mkdir -p " + dirName);
  TFile* outputFile = TFile::Open(dirName + "/histos_" + outputName  + ".root","RECREATE"); // CB find a better name for output file  

  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 

  // Initialize Root application
  TRint* app = new TRint("CMS Root Application", &argc, argv);

  setTDRStyle();
 
  TagAndProbeConfig tnpConfig;
  std::vector<SampleConfig> sampleConfigs;
  
  parseConfig(configFile,tnpConfig,sampleConfigs);

  std::vector<Plotter> plotters;

  for (auto sampleConfig : sampleConfigs)
    {

      Plotter plotter(tnpConfig, sampleConfig);
      plotter.book(outputFile);
      
      plotters.push_back(plotter);
    }
 
  for (auto plotter : plotters)
    {


      TString fileName = plotter.m_sampleConfig.fileName;
      std::cout << "[" << argv[0] << "] Processing file "
		<< fileName.Data() << std::endl;  
  
      // Initialize pointers to summary and full event structure

      muon_pog::Event*   ev   = new muon_pog::Event();

      TTree* tree;
      TBranch* evBranch;

      // Open file, get tree, set branches

      TFile* inputFile = TFile::Open(fileName,"READONLY");
      tree = (TTree*)inputFile->Get("MUONPOGTREE");
      if (!tree) inputFile->GetObject("MuonPogTree/MUONPOGTREE",tree);

      evBranch = tree->GetBranch("event");
      evBranch->SetAddress(&ev);

      // Watch number of entries
      int nEntries = plotter.m_sampleConfig.nEvents > 0 ? plotter.m_sampleConfig.nEvents : tree->GetEntriesFast();

      std::cout << "[" << argv[0] << "] Number of entries = " << nEntries << std::endl;

      int nFilteredEvents = 0;

      float weight = 1.;
      // Add cross section weight -> To take into account QCD filter efficiencies
      if(plotter.m_sampleConfig.applyReweighting==true)
	weight *= (plotter.m_sampleConfig.QCDWeight);
      
      std::cout << "Weight per event = " << weight << std::endl;	

      for (Long64_t iEvent=0; iEvent<nEntries; ++iEvent) {
	if (tree->LoadTree(iEvent)<0) break;
	
	print_progress(nEntries, iEvent);
	
	evBranch->GetEntry(iEvent);

	plotter.fill(ev->muons, ev->hlt, (*ev), weight);
	
	plotter.fillGen(ev->genParticles, (*ev));
      }
      
      std::cout << "[==================================================] 100% " << std::endl;
      
      delete ev;
      inputFile->Close();
      std::cout << std::endl;
	   
    }
  
  muon_pog::comparisonPlots(plotters,outputFile,dirName);
  muon_pog::copyPhp(dirName);
  
  outputFile->Write();
  
  if (!gROOT->IsBatch()) app->Run();

  return 0;

}


muon_pog::TagAndProbeConfig::TagAndProbeConfig(boost::property_tree::ptree::value_type & vt)
{

  try
    {

      gen_DrCut      = vt.second.get<Float_t>("gen_DrCut");
     
      muon_trackType = vt.second.get<std::string>("muon_trackType");

      probe_minPt  = vt.second.get<Float_t>("probe_minPt");
      probe_isoLooseWP = vt.second.get<Float_t>("probe_isoLooseWP");
      probe_isoTightWP = vt.second.get<Float_t>("probe_isoTightWP");

    }

  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[TagAndProbeConfig] Can't get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }

}

muon_pog::SampleConfig::SampleConfig(boost::property_tree::ptree::value_type & vt)
{
 
 try
    {
      fileName     = TString(vt.second.get<std::string>("fileName").c_str());
      sampleName   = TString(vt.first.c_str());
      QCDWeight = vt.second.get<Float_t>("QCDWeight");
      nEvents = vt.second.get<Float_t>("nEvents"); //CB do we really need this? can't we take nEvents from the file itself?
      applyReweighting = vt.second.get<Bool_t>("applyReweighting");
      runs = toArray(vt.second.get<std::string>("runs"));
    }
  
  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[TagAndProbeConfig] Can't get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }
  
}

std::vector<int> muon_pog::SampleConfig::toArray(const std::string& entries)
{
  
  std::vector<int> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(atoi(item.c_str()));
  return result;

}


std::vector<TString> muon_pog::TagAndProbeConfig::toArray(const std::string& entries)
{
  
  std::vector<TString> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(TString(item));
  return result;

}


void muon_pog::Plotter::book(TFile *outFile)
{

  const double etaBins[9] = {-2.4,-1.8,-1.2,-0.6,0.0,
   			      0.6, 1.2, 1.8, 2.4};
  const double pTBins[7] = {5,15,25,35,50,70,100};
  const double VtxBins[10] = {0,10,20,24,28,32,36,40,50,60};
	  
  TString sampleTag = m_sampleConfig.sampleName;
  
  outFile->cd("/");
  outFile->mkdir(sampleTag);
  outFile->cd(sampleTag);
  
  outFile->mkdir(sampleTag + "/Yields");
 
  TH1::SetDefaultSumw2(kTRUE);
  
  // -- Yields
  outFile->cd(sampleTag + "/Yields/");
  m_2Dyields["RecoMuon"] = new TH2D("Yields_RecoMuon" ,"Yields " + sampleTag + " RecoMuons ",10,0,10,4,0,4);
  m_2Dyields["RecoMuon"]->SetOption("COLTEXT"); 
  for (unsigned nid=0; nid<7; nid++) m_2Dyields["RecoMuon"]->GetXaxis()->SetBinLabel(nid+1,IDName[nid]);
  m_2Dyields["RecoMuon"]->GetXaxis()->SetBinLabel(8, "Total Pass Cut");
  m_2Dyields["RecoMuon"]->GetXaxis()->SetBinLabel(9, "Total");
  m_2Dyields["RecoMuon"]->GetXaxis()->SetBinLabel(10, "Total Events (weighted)");
  for (unsigned nre=0; nre<4; nre++) m_2Dyields["RecoMuon"]->GetYaxis()->SetBinLabel(nre+1,region[nre]);
  
  m_2Dyields["GenMuon"] = new TH2D("Yields_GenMuon" ,"Yields GenMuons " + sampleTag,4,0,4,4,0,4);
  m_2Dyields["GenMuon"]->SetOption("COLTEXT"); 
  m_2Dyields["GenMuon"]->GetXaxis()->SetBinLabel(1, "Total Events");
  m_2Dyields["GenMuon"]->GetXaxis()->SetBinLabel(2, "Gen Muon Z");
  m_2Dyields["GenMuon"]->GetXaxis()->SetBinLabel(3, "Gen Muon Others");
  m_2Dyields["GenMuon"]->GetXaxis()->SetBinLabel(4, "Gen Muon Total");
  for (unsigned nre=0; nre<4; nre++) m_2Dyields["GenMuon"]->GetYaxis()->SetBinLabel(nre+1,region[nre]);
  
  
  for (unsigned nid=0; nid<7; nid++){
    outFile->mkdir(sampleTag + "/KinIso_variables/" + IDName[nid]);
    for (unsigned nre=0; nre<4; nre++){
      outFile->mkdir(sampleTag + "/KinIso_variables/" + IDName[nid] + "/" + region[nre]);
      for (unsigned npu=0; npu<4; npu++){
	outFile->mkdir(sampleTag + "/KinIso_variables/" + IDName[nid] + "/" + region[nre] + "/" + PUr[npu]);
	outFile->cd   (sampleTag + "/KinIso_variables/" + IDName[nid] + "/" + region[nre] + "/" + PUr[npu]);
      
	TString IDRegName  = IDName[nid] + "_" + region[nre] + "_" + PUr[npu];
	TString IDRegTitle = "ID:" + IDName[nid] + " - REGION:" + region[nre] + " - PU:" + PUr[npu];
	
	m_plots[IDName[nid]+region[nre]+PUr[npu]+"RelIso"]       = new TH1D ("RelIso_" + IDRegName,     IDRegTitle + " Relative Isolation ; Relative Iso.; # entries", 100, 0., 1.);
	m_plots[IDName[nid]+region[nre]+PUr[npu]+"RelIsoNoDB"]   = new TH1D ("RelIsoNoDB_" + IDRegName, IDRegTitle + " Relative Isolation w/o #Delta#beta ; Relative Iso. (No #Delta#beta); # entries", 100, 0., 2.);

	m_prof[IDName[nid]+region[nre]+PUr[npu]+"VtxVsRelIso"]   = new TProfile("VtxVsRelIso_" + IDRegName, IDRegTitle + " Vtx Vs Relative Isolation ; N. Vtx; Relative Iso.", 80, 0., 80., 0., 1.);
	m_prof[IDName[nid]+region[nre]+PUr[npu]+"PtVsRelIso"]    = new TProfile("PtVsRelIso_"  + IDRegName, IDRegTitle + " Pt Vs Relative Isolation ; p_{T} [GeV]; Relative Iso.", 100,0., 100, 0., 1.);
	m_prof[IDName[nid]+region[nre]+PUr[npu]+"EtaVsRelIso"]   = new TProfile("EtaVsRelIso_" + IDRegName, IDRegTitle + " Eta Vs Relative Isolation ; #eta; Relative Iso.", 48, -2.4, 2.4, 0., 1.);

	m_prof[IDName[nid]+region[nre]+PUr[npu]+"VtxVsRelIsoNoDB"]   = new TProfile("VtxVsRelIsoNoDB_" + IDRegName, IDRegTitle + " Vtx Vs Relative Isolation  w/o #Delta#beta ; N. Vtx; Relative Iso. (No #Delta#beta)", 80, 0., 80., 0., 1.);
	m_prof[IDName[nid]+region[nre]+PUr[npu]+"PtVsRelIsoNoDB"]    = new TProfile("PtVsRelIsoNoDB_"  + IDRegName, IDRegTitle + " Pt Vs Relative Isolation  w/o #Delta#beta ; p_{T} [GeV]; Relative Iso. (No #Delta#beta)", 100,0., 100, 0., 1.);
	m_prof[IDName[nid]+region[nre]+PUr[npu]+"EtaVsRelIsoNoDB"]   = new TProfile("EtaVsRelIsoNoDB_" + IDRegName, IDRegTitle + " Eta Vs Relative Isolation  w/o #Delta#beta ; #eta; Relative Iso. (No #Delta#beta)", 48, -2.4, 2.4, 0., 1.);
	
	m_plots[IDName[nid]+region[nre]+PUr[npu]+"RelIsoCut"] = new TH1D ("RelIsoCut_" + IDRegName, IDRegTitle + " Relative Isolation Cut ; Relative Iso. Cut; # entries", 101, 0., 1.01);
	
	for (unsigned niso=0; niso<=2; niso++){
	  
	  IDRegName  = IDName[nid] + "_" + region[nre] + "_" + PUr[npu] + "_" + isowp[niso];
	  IDRegTitle = "ID:" + IDName[nid] + " - REGION:" + region[nre] + " - PU:" + PUr[npu] + " - ISO:" + isowp[niso];      
	  
	  m_plots[IDName[nid]+region[nre]+PUr[npu]+isowp[niso]+"Vtx"] = new TH1D ("Vtx_" + IDRegName,  IDRegTitle + " Vtx ; N. Vtx; # entries", 80, 0., 80.);
	  m_plots[IDName[nid]+region[nre]+PUr[npu]+isowp[niso]+"Pt"]  = new TH1D ("Pt_" + IDRegName,  IDRegTitle + " Pt ; p_{T} [GeV]; # entries", 100,0., 100);
	  m_plots[IDName[nid]+region[nre]+PUr[npu]+isowp[niso]+"Eta"] = new TH1D ("Eta_" + IDRegName,  IDRegTitle + " Eta ; #eta; # entries", 48, -2.4, 2.4);
	  
	  m_plots[IDName[nid]+region[nre]+PUr[npu]+isowp[niso]+"ChHadIso"] = new TH1D ("ChHadIso_" + IDRegName,  IDRegTitle + " Charged Hadron Isolation ; Charged Had. Iso.; # entries", 40, 0., 10);
	  m_plots[IDName[nid]+region[nre]+PUr[npu]+isowp[niso]+"NeHadIso"] = new TH1D ("NeHadIso_" + IDRegName,  IDRegTitle + " Neutral Hadron Isolation ; Neutral Had. Iso.; # entries", 40, 0., 10);
	  m_plots[IDName[nid]+region[nre]+PUr[npu]+isowp[niso]+"PhIso"]    = new TH1D ("PhIso_" + IDRegName,     IDRegTitle + " Photon Isolation ; Photon Iso.; # entries", 40, 0., 60);
	  m_plots[IDName[nid]+region[nre]+PUr[npu]+isowp[niso]+"PUIso"]    = new TH1D ("PUIso_" + IDRegName,     IDRegTitle + " PU Isolation ; PU Iso.; # entries", 40, 0., 60);
	  
	  m_plots[IDName[nid]+region[nre]+PUr[npu]+isowp[niso]+"TotalIso"]     = new TH1D ("TotalIso_" + IDRegName,  IDRegTitle + " Total Isolation ; Total Iso.; # entries", 100, 0., 200);
	  m_plots[IDName[nid]+region[nre]+PUr[npu]+isowp[niso]+"TotalIsoNoDB"] = new TH1D ("TotalIsoNoDB_" + IDRegName,  IDRegTitle + " Total Isolation w/o #Delta#beta ; Total Iso. (No #Delta#beta); # entries", 100, 0., 200);

	  m_prof[IDName[nid]+region[nre]+PUr[npu]+isowp[niso]+"VtxVsTotalIso"]   = new TProfile("VtxVsTotalIso_" + IDRegName, IDRegTitle + " Vtx Vs Total Isolation ; N. Vtx; Total Iso.", 80, 0., 80., 0., 200.);
	  m_prof[IDName[nid]+region[nre]+PUr[npu]+isowp[niso]+"PtVsTotalIso"]    = new TProfile("PtVsTotalIso_"  + IDRegName, IDRegTitle + " Pt Vs Total Isolation ; p_{T} [GeV]; Total Iso.", 100,0., 100, 0., 200.);
	  m_prof[IDName[nid]+region[nre]+PUr[npu]+isowp[niso]+"EtaVsTotalIso"]   = new TProfile("EtaVsTotalIso_" + IDRegName, IDRegTitle + " Eta Vs Total Isolation ; #eta; Total Iso.", 48, -2.4, 2.4, 0., 200.);
	  
	  m_prof[IDName[nid]+region[nre]+PUr[npu]+isowp[niso]+"VtxVsTotalIsoNoDB"]   = new TProfile("VtxVsTotalIsoNoDB_" + IDRegName, IDRegTitle + " Vtx Vs Total Isolation  w/o #Delta#beta ; N. Vtx; Total Iso. (No #Delta#beta)", 80, 0., 80., 0., 200.);
	  m_prof[IDName[nid]+region[nre]+PUr[npu]+isowp[niso]+"PtVsTotalIsoNoDB"]    = new TProfile("PtVsTotalIsoNoDB_"  + IDRegName, IDRegTitle + " Pt Vs Total Isolation  w/o #Delta#beta ; p_{T} [GeV]; Total Iso. (No #Delta#beta)", 100,0., 100, 0., 200.);
	  m_prof[IDName[nid]+region[nre]+PUr[npu]+isowp[niso]+"EtaVsTotalIsoNoDB"]   = new TProfile("EtaVsTotalIsoNoDB_" + IDRegName, IDRegTitle + " Eta Vs Total Isolation  w/o #Delta#beta ; #eta; Total Iso. (No #Delta#beta)", 48, -2.4, 2.4, 0., 200.);


	}// for(npu)
      }// for(niso)
    } // for(nre)
  } // for(nid)
  
} // void::book 

void muon_pog::Plotter::fillGen(const std::vector<muon_pog::GenParticle> & genpars, const muon_pog::Event & ev){
  
  float etaBarrel = 1.2;  
  float pTeta = 20.0;
  
  // Number of Events 
  m_2Dyields["GenMuon"]->Fill(0.,0.);
  
  for (auto & genpar : genpars){

    bool genZmuon     = false;
    bool genOthermuon = false;

    if(abs(genpar.pdgId) == 13 &&
       genpar.pt          > m_tnpConfig.probe_minPt &&
       fabs(genpar.eta)   < 2.4 ){

      bool IsMuonpTCorr = false;
      // Check this with flags!!!!
      for (auto motherId : genpar.mothers){
	if(abs(motherId) == 13) IsMuonpTCorr =true;
      }      
      if (IsMuonpTCorr) continue;

      // Z_id = 23
      genZmuon = hasMother(genpar, 23);
      if(!genZmuon) genOthermuon = true;
    }
    if(genZmuon || genOthermuon){
      
      for (unsigned int nre = 0; nre<4; nre++){
	
	bool etaRegion = false;
	if (nre == 0) etaRegion = true;
	if (nre == 1 && std::abs(genpar.eta) <  etaBarrel) etaRegion = true;
	if (nre == 2 && std::abs(genpar.eta) >= etaBarrel) etaRegion = true;
	if (nre == 3 && (std::abs(genpar.eta) > 0.9 && std::abs(genpar.eta) < 1.3)) etaRegion = true;
	
	if (!etaRegion)     continue;

	if (genZmuon)     m_2Dyields["GenMuon"]->Fill(1,nre); // From Z
	if (genOthermuon) m_2Dyields["GenMuon"]->Fill(2,nre); 
	m_2Dyields["GenMuon"]->Fill(3,nre); // Number of GenMuons
	
      }// for(nre)
    }//if(GenMuon)
    
  }// for(genpar)  
}

void muon_pog::Plotter::fill(const std::vector<muon_pog::Muon> & muons,
			     const muon_pog::HLT & hlt, const muon_pog::Event & ev, float weight)
{

  // Total Number of Events 
  m_2Dyields["RecoMuon"]->Fill(9., 0., weight);
  
  float etaBarrel = 1.2;  
  float pTeta = 20.0;

  for (auto & muon : muons){

    // Total Number of Muons
    m_2Dyields["RecoMuon"]->Fill(8., 0., weight);
    
    // General Probe Muons	      
    if(muon.pt > m_tnpConfig.probe_minPt && 
       fabs(muon.eta) < 2.4){
      
      // Number of Muons Passing (pT,eta) cuts
      m_2Dyields["RecoMuon"]->Fill(7., 0., weight);
      
      bool IsSIGN = false;
      
      // Check first if muon comes from Z
      // Check "hasGenMatch" and "hasNoGenMatch"
      IsSIGN = muon_pog::hasGenMatch(muon, ev.genParticles, m_tnpConfig.gen_DrCut, 23);

      bool FillMuon = false;
      TString sampleTag = m_sampleConfig.sampleName;
      if(sampleTag.Contains("SIGNAL") && IsSIGN) FillMuon = true;
      if(sampleTag.Contains("BACKGROUND") && !IsSIGN) FillMuon = true;

      if (!FillMuon) continue;
      
      bool IsMuonID[7];
      for (int bmid = 0; bmid<7; bmid++) IsMuonID[bmid] = false;
      IsMuonID[0] = muon.isGlobal ; 
      IsMuonID[1] = muon.isTracker ; 
      IsMuonID[2] = muon.isLoose ; 
      IsMuonID[3] = muon.isMedium ; 
      IsMuonID[4] = muon.isTight ; 
      IsMuonID[5] = muon.isSoft ; 
      IsMuonID[6] = muon.isHighPt ; 

      // Muon Track
      TLorentzVector muTk(muon_pog::muonTk(muon, m_tnpConfig.muon_trackType));

      bool IsEtaRegion[4];
      for (int bnre = 0; bnre<4; bnre++) IsEtaRegion[bnre] = false;
      IsEtaRegion[0] = true;
      IsEtaRegion[1] = (std::abs(muTk.Eta()) < etaBarrel);
      IsEtaRegion[2] = (std::abs(muTk.Eta()) >= etaBarrel);
      IsEtaRegion[3] = ((std::abs(muTk.Eta()) > 0.9 && std::abs(muTk.Eta()) < 1.3));

      bool IsPURegime[4];
      for (int bnpu = 0; bnpu<4; bnpu++) IsPURegime[bnpu] = false;
      IsPURegime[0] = true;
      IsPURegime[1] = (ev.nVtx < 25);
      IsPURegime[2] = (ev.nVtx >= 25 && ev.nVtx < 45);
      IsPURegime[3] = (ev.nVtx >= 45);
      
      for (unsigned int mid = 0; mid<7; mid++){

	if (!IsMuonID[mid])    continue;

	for (unsigned int nre = 0; nre<4; nre++){

	  if (!IsEtaRegion[nre]) continue;
	  
	  m_2Dyields["RecoMuon"]->Fill(mid,nre, weight);	  
	  
	  for (unsigned int npu = 0; npu<4; npu++){
	    
	    if (!IsPURegime[npu])  continue;
	    
	    // Isolation (R=0.4)
	    float ChHadIso = muon.chargedHadronIso;
	    float NeHadIso = muon.neutralHadronIso;
	    float PhIso    = muon.photonIso;
	    float ChPUIso  = muon.chargedHadronIsoPU;
	    
	    float TotalIso = (ChHadIso + std::max(0., NeHadIso + PhIso - 0.5*ChPUIso));
	    float RelIso   = TotalIso/muon.pt; // muon.isoPflow04;	  
	    
	    float TotalIsoNoDB = ChHadIso + NeHadIso + PhIso;
	    float RelIsoNoDB   = TotalIsoNoDB/muon.pt; // muon.isoPflow04;	  
	    
	    m_plots[IDName[mid]+region[nre]+PUr[npu]+"RelIso"]->Fill(RelIso, weight);
	    m_plots[IDName[mid]+region[nre]+PUr[npu]+"RelIsoNoDB"]->Fill(RelIsoNoDB, weight);

	    m_prof[IDName[mid]+region[nre]+PUr[npu]+"VtxVsRelIso"]->Fill(ev.nVtx, RelIso, weight);
	    m_prof[IDName[mid]+region[nre]+PUr[npu]+"PtVsRelIso"] ->Fill(muTk.Pt(), RelIso, weight);
	    if(muTk.Pt() > pTeta) m_prof[IDName[mid]+region[nre]+PUr[npu]+"EtaVsRelIso"]->Fill(muTk.Eta(), RelIso, weight);

	    m_prof[IDName[mid]+region[nre]+PUr[npu]+"VtxVsRelIsoNoDB"]->Fill(ev.nVtx, RelIsoNoDB, weight);
	    m_prof[IDName[mid]+region[nre]+PUr[npu]+"PtVsRelIsoNoDB"] ->Fill(muTk.Pt(), RelIsoNoDB, weight);
	    if(muTk.Pt() > pTeta) m_prof[IDName[mid]+region[nre]+PUr[npu]+"EtaVsRelIsoNoDB"]->Fill(muTk.Eta(), RelIsoNoDB, weight);
	    
	    int vRelIsoCut = 1 + std::round(RelIso*100.);
	    for (int ibin=vRelIsoCut; ibin<=100; ibin++) m_plots[IDName[mid]+region[nre]+PUr[npu]+"RelIsoCut"]->AddBinContent(ibin, weight);
	    m_plots[IDName[mid]+region[nre]+PUr[npu]+"RelIsoCut"]->AddBinContent(101, weight); // Total Number of Muons 
	    
	    for(int isocat=0;isocat<=2;isocat++){
	      bool IsIso = false;
	      if(isocat == 0) IsIso = true;
	      if(isocat == 1 && RelIso<m_tnpConfig.probe_isoLooseWP) IsIso = true;
	      if(isocat == 2 && RelIso<m_tnpConfig.probe_isoTightWP) IsIso = true;
	      
	      if(!IsIso) continue;
	      
	      m_plots[IDName[mid]+region[nre]+PUr[npu]+isowp[isocat]+"Vtx"]->Fill(ev.nVtx, weight);
	      m_plots[IDName[mid]+region[nre]+PUr[npu]+isowp[isocat]+"Pt"] ->Fill(muTk.Pt(), weight);
	      if(muTk.Pt() > pTeta) m_plots[IDName[mid]+region[nre]+PUr[npu]+isowp[isocat]+"Eta"]->Fill(muTk.Eta(), weight);
	      
	      m_plots[IDName[mid]+region[nre]+PUr[npu]+isowp[isocat]+"ChHadIso"]->Fill(ChHadIso, weight);
	      m_plots[IDName[mid]+region[nre]+PUr[npu]+isowp[isocat]+"NeHadIso"]->Fill(NeHadIso, weight);
	      m_plots[IDName[mid]+region[nre]+PUr[npu]+isowp[isocat]+"PhIso"]   ->Fill(PhIso, weight);
	      m_plots[IDName[mid]+region[nre]+PUr[npu]+isowp[isocat]+"PUIso"]   ->Fill(ChPUIso, weight);
	      
	      m_plots[IDName[mid]+region[nre]+PUr[npu]+isowp[isocat]+"TotalIso"]->Fill(TotalIso, weight);
	      m_plots[IDName[mid]+region[nre]+PUr[npu]+isowp[isocat]+"TotalIsoNoDB"]->Fill(TotalIsoNoDB, weight);

	      m_prof[IDName[mid]+region[nre]+PUr[npu]+isowp[isocat]+"VtxVsTotalIso"]->Fill(ev.nVtx, TotalIso, weight);
	      m_prof[IDName[mid]+region[nre]+PUr[npu]+isowp[isocat]+"PtVsTotalIso"] ->Fill(muTk.Pt(), TotalIso, weight);
	      if(muTk.Pt() > pTeta) m_prof[IDName[mid]+region[nre]+PUr[npu]+isowp[isocat]+"EtaVsTotalIso"]->Fill(muTk.Eta(), TotalIso, weight);
	      
	      m_prof[IDName[mid]+region[nre]+PUr[npu]+isowp[isocat]+"VtxVsTotalIsoNoDB"]->Fill(ev.nVtx, TotalIsoNoDB, weight);
	      m_prof[IDName[mid]+region[nre]+PUr[npu]+isowp[isocat]+"PtVsTotalIsoNoDB"] ->Fill(muTk.Pt(), TotalIsoNoDB, weight);
	      if(muTk.Pt() > pTeta) m_prof[IDName[mid]+region[nre]+PUr[npu]+isowp[isocat]+"EtaVsTotalIsoNoDB"]->Fill(muTk.Eta(), TotalIsoNoDB, weight);
	      
	    } // for(IsoCategory)
	    	    
	  } // for(PURegime)
	
	} // for(etaRegion)
      
      } // for(muonID)
      
    } // if(muon pT eta)
  }// for(Muons)
  
  
} // void::fill


//Functions

void muon_pog::print_progress(int TreeEntries, Long64_t ievt){
  if(TreeEntries < 50) TreeEntries = 50;
  int step = TreeEntries/50;
  if (ievt%(step) == 0){ 
    float progress=(ievt)/(TreeEntries*1.0);
    int barWidth = 50;
    
    std::cout << "[";
    int pos = barWidth * progress;
    
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
    }
    
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
  }  
}


void muon_pog::parseConfig(const std::string configFile, muon_pog::TagAndProbeConfig & tpConfig,
			   std::vector<muon_pog::SampleConfig> & sampleConfigs)
{
  
  boost::property_tree::ptree pt;
  
  try
    {
      boost::property_tree::ini_parser::read_ini(configFile, pt);
    }
  catch (boost::property_tree::ini_parser::ini_parser_error iniParseErr)
    {
      std::cout << "[TagAndProbeConfig] Can't open : " << iniParseErr.filename()
		<< "\n\tin line : " << iniParseErr.line()
		<< "\n\thas error :" << iniParseErr.message()
		<< std::endl;
      throw std::runtime_error("Bad INI parsing");
    }

  for( auto vt : pt )
    {
      if (vt.first.find("TagAndProbe") != std::string::npos)
	tpConfig = muon_pog::TagAndProbeConfig(vt);
      else
	sampleConfigs.push_back(muon_pog::SampleConfig(vt));
    }
}

void muon_pog::comparisonPlots(std::vector<muon_pog::Plotter> & plotters,
 			       TFile *outFile, TString &  outputDir)
{

  outFile->cd("/");
  
}

void muon_pog::copyPhp(const TString &  outputDir)
{
  
  system("cp index.php " + outputDir);
  
  boost::filesystem::directory_iterator dirIt(outputDir.Data());
  boost::filesystem::directory_iterator dirEnd;
  for (;dirIt != dirEnd; ++ dirIt)
    {
      if (boost::filesystem::is_directory(dirIt->status()))
	copyPhp(TString(dirIt->path().string()));
    }

}

//  LocalWords:  IsoTight
