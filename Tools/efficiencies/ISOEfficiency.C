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
std::vector<TString> GetListOfFiles(TString FileName);

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
    
    Float_t gen_DrCut;
      
    Float_t probe_minPt;
    Float_t iso_DeltaBeta;      
    Float_t iso_LooseWP, iso_TightWP;
    Float_t isoPUPPI_LooseWP, isoPUPPI_TightWP, isoPUPPILep_LooseWP, isoPUPPILep_TightWP, isoPUPPINoLep_LooseWP, isoPUPPINoLep_TightWP;
    Float_t Miniiso_LooseWP, Miniiso_TightWP;      
    Float_t IDMedium_CutdXY, IDMedium_CutdZ;
   
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
    
    std::vector<TString> IDName = {"GLB",  "TRK", "LOOSE", "MEDIUM", "TIGHT", "SOFT", "HIGHpT", "MEDIUMPrompt",
				   "MVALOOSE", "MVAMEDIUM", "MVATIGHT"};
    enum BinType   {GLB=0, TRK, LOOSE, MEDIUM, TIGHT, SOFT, HIGHpT, MEDIUMPrompt,
		    MVALOOSE, MVAMEDIUM, MVATIGHT};

    std::vector<TString> region = {"Full", "Barrel", "Endcap", "Overlap"};
    enum BinRegion {Full=0, Barrel, Endcap, Overlap};
    
    std::vector<TString> IsoName  = {"ISOPF","ISOPUPPI","ISOPUPPILep","ISOPUPPINoLep","MiniISO"};
    enum BinIsoName  {ISOPF=0,ISOPUPPI,ISOPUPPILep,ISOPUPPINoLep,MiniISO};
    
    std::vector< std::vector<TString> > isowp  = {{"ISONoCut", "ISOLoose", "ISOTight"},
						  {"PUPPILoose","PUPPITight"},
						  {"PUPPILepLoose","PUPPILepTight"},
						  {"PUPPINoLepLoose","PUPPINoLepTight"},
						  {"MiniISOLoose","MiniISOTight"}};
    
    std::vector<TString> PUr = {"FullPU", "LowPU", "MediumPU", "HighPU"};  
    enum BinPUr{FullPU=0, LowPU, MediumPU, HighPU};

    std::vector<TString> pTr = {"FullpT", "LowpT", "MediumpT", "HighpT"};  
    enum BinpTr{FullpT=0, LowpT, MediumpT, HighpT};
    
    
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


      TString fileListName = plotter.m_sampleConfig.fileName;
      std::cout << "[" << argv[0] << "] Processing file "
		<< fileListName.Data() << std::endl;  

      std::vector<TString> SetOfFiles = GetListOfFiles(fileListName);
  
      std::cout << "Number of files to be processed: " << SetOfFiles.size() << std::endl;
      
      // Initialize pointers to summary and full event structure
      muon_pog::Event*   ev   = new muon_pog::Event();

      // Loop over all files
      for(unsigned int ifile = 0; ifile < SetOfFiles.size(); ifile++){
	
	TFile* inputFile; 
	TTree* tree;
	TBranch* evBranch;
	
	// Open file, get tree, set branches
	inputFile = TFile::Open(SetOfFiles.at(ifile),"READONLY");
	tree = (TTree*)inputFile->Get("MUONPOGTREE");
	if (!tree) inputFile->GetObject("MuonPogTree/MUONPOGTREE",tree);
	
	evBranch = tree->GetBranch("event");
	evBranch->SetAddress(&ev);
	
	// Watch number of entries
	int nEntries;
	if (plotter.m_sampleConfig.nEvents > 0 && 
	    plotter.m_sampleConfig.nEvents < tree->GetEntriesFast()) nEntries = plotter.m_sampleConfig.nEvents;
	else nEntries = tree->GetEntriesFast();
	
	std::cout << "[" << argv[0] << "] Number of entries/sample = " << nEntries << std::endl;
	
	int nFilteredEvents = 0;
	
	float weight = 1.;
	// Add cross section weight -> To take into account QCD filter efficiencies
	if(plotter.m_sampleConfig.applyReweighting==true)
	  weight *= (plotter.m_sampleConfig.QCDWeight/nEntries);
      
	std::cout << "Weight per event = " << weight << std::endl;	
	
	for (Long64_t iEvent=0; iEvent<nEntries; ++iEvent) {
	  if (tree->LoadTree(iEvent)<0) break;
	
	  print_progress(nEntries, iEvent);
	 
	  evBranch->GetEntry(iEvent);

	  plotter.fill(ev->muons, ev->hlt, (*ev), weight);
	
	  plotter.fillGen(ev->genParticles, (*ev));
	}
      
	std::cout << "[==================================================] 100% " << std::endl;
	
	delete tree;
	inputFile->Close();
	delete inputFile;
      } // for(ifiles)
      
      delete ev;
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
     
      iso_DeltaBeta  = vt.second.get<Float_t>("iso_DeltaBeta");

      iso_LooseWP           = vt.second.get<Float_t>("iso_LooseWP");
      iso_TightWP           = vt.second.get<Float_t>("iso_TightWP");
      isoPUPPI_LooseWP      = vt.second.get<Float_t>("isoPUPPI_LooseWP");
      isoPUPPI_TightWP      = vt.second.get<Float_t>("isoPUPPI_TightWP");
      isoPUPPILep_LooseWP   = vt.second.get<Float_t>("isoPUPPILep_LooseWP");
      isoPUPPILep_TightWP   = vt.second.get<Float_t>("isoPUPPILep_TightWP");
      isoPUPPINoLep_LooseWP = vt.second.get<Float_t>("isoPUPPINoLep_LooseWP");
      isoPUPPINoLep_TightWP = vt.second.get<Float_t>("isoPUPPINoLep_TightWP");
      Miniiso_LooseWP       = vt.second.get<Float_t>("Miniiso_LooseWP");
      Miniiso_TightWP       = vt.second.get<Float_t>("Miniiso_TightWP");

      IDMedium_CutdXY = vt.second.get<Float_t>("IDMedium_CutdXY");
      IDMedium_CutdZ  = vt.second.get<Float_t>("IDMedium_CutdZ");

      probe_minPt  = vt.second.get<Float_t>("probe_minPt");

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
	  
  TString sampleTag = m_sampleConfig.sampleName;
  
  outFile->cd("/");
  outFile->mkdir(sampleTag);
  outFile->cd(sampleTag);
  
  outFile->mkdir(sampleTag + "/Yields");
 
  TH1::SetDefaultSumw2(kTRUE);
  
  // -- Yields
  outFile->cd(sampleTag + "/Yields/");
  m_2Dyields["RecoMuon"] = new TH2D("Yields_RecoMuon" ,"Yields " + sampleTag + " RecoMuons ",IDName.size()+3,0,IDName.size()+3, region.size(),0,region.size());
  m_2Dyields["RecoMuon"]->SetOption("COLTEXT"); 
  for (unsigned int nid=0; nid<IDName.size(); nid++) m_2Dyields["RecoMuon"]->GetXaxis()->SetBinLabel(nid+1,IDName[nid]);
  m_2Dyields["RecoMuon"]->GetXaxis()->SetBinLabel(IDName.size()+1, "Total Pass Cut");
  m_2Dyields["RecoMuon"]->GetXaxis()->SetBinLabel(IDName.size()+2, "Total");
  m_2Dyields["RecoMuon"]->GetXaxis()->SetBinLabel(IDName.size()+3, "Total Events (weighted)");
  for (unsigned int nre=0; nre<region.size(); nre++) m_2Dyields["RecoMuon"]->GetYaxis()->SetBinLabel(nre+1,region[nre]);
  
  m_2Dyields["GenMuon"] = new TH2D("Yields_GenMuon" ,"Yields GenMuons " + sampleTag,4,0,4,region.size(),0,region.size());
  m_2Dyields["GenMuon"]->SetOption("COLTEXT"); 
  m_2Dyields["GenMuon"]->GetXaxis()->SetBinLabel(1, "Total Events");
  m_2Dyields["GenMuon"]->GetXaxis()->SetBinLabel(2, "Gen Muon Z");
  m_2Dyields["GenMuon"]->GetXaxis()->SetBinLabel(3, "Gen Muon Others");
  m_2Dyields["GenMuon"]->GetXaxis()->SetBinLabel(4, "Gen Muon Total");
  for (unsigned int nre=0; nre<region.size(); nre++) m_2Dyields["GenMuon"]->GetYaxis()->SetBinLabel(nre+1,region[nre]);

  m_plots["NumberOfGenMuonsFromZ"]  = new TH1D ("NMuons_GenFromZ", " Number of GEN Muons from Z per Evt ; N. Muons; # entries", 4, 0., 4.);
  m_plots["NumberOfGenMuonsOther"]  = new TH1D ("NMuons_GenOther", " Number of GEN Muons (NOT from Z) per Evt ; N. Muons; # entries", 4, 0., 4.);

  m_plots["NMuons"] = new TH1D ("NMuons", " Number of Muons per Evt ; N. Muons; # entries", 4, 0., 4.);

  // Vertices
  for (unsigned int npu=0; npu<4; npu++){
    outFile->mkdir(sampleTag + "/Vtx/" + PUr[npu]);
    outFile->cd   (sampleTag + "/Vtx/" + PUr[npu]);
    // Reco Vertex        
    m_plots[PUr[npu]+"RecoVtx_X"] = new TH1D ("RecoVtx_X_" + PUr[npu], "PU:" + PUr[npu] + " Vertex position (RECO) ; X [cm]; # entries", 200, 0., 0.2);
    m_plots[PUr[npu]+"RecoVtx_Y"] = new TH1D ("RecoVtx_Y_" + PUr[npu], "PU:" + PUr[npu] + " Vertex position (RECO) ; Y [cm]; # entries", 200, 0., 0.2);
    m_plots[PUr[npu]+"RecoVtx_Z"] = new TH1D ("RecoVtx_Z_" + PUr[npu], "PU:" + PUr[npu] + " Vertex position (RECO) ; Z [cm]; # entries", 240, 0., 24);
    m_plots[PUr[npu]+"RecoVtx_XY"] = new TH1D ("RecoVtx_XY_" + PUr[npu], "PU:" + PUr[npu] + " Vertex position (RECO) ; XY [cm]; # entries", 400, 0., 0.4);
    // Gen Vertex
    m_plots[PUr[npu]+"GenVtx_X"] = new TH1D ("GenVtx_X_" + PUr[npu], "PU:" + PUr[npu] + " Vertex position (GEN) ; X [cm]; # entries", 200, 0., 0.2);
    m_plots[PUr[npu]+"GenVtx_Y"] = new TH1D ("GenVtx_Y_" + PUr[npu], "PU:" + PUr[npu] + " Vertex position (GEN) ; Y [cm]; # entries", 200, 0., 0.2);
    m_plots[PUr[npu]+"GenVtx_Z"] = new TH1D ("GenVtx_Z_" + PUr[npu], "PU:" + PUr[npu] + " Vertex position (GEN) ; Z [cm]; # entries", 240, 0., 24);
    m_plots[PUr[npu]+"GenVtx_XY"] = new TH1D ("GenVtx_XY_" + PUr[npu], "PU:" + PUr[npu] + " Vertex position (GEN) ; XY [cm]; # entries", 400, 0., 0.4);
    // Delta Vertex (Reco-Gen)
    m_plots[PUr[npu]+"Vtx_DX"] = new TH1D ("Vtx_DX_" + PUr[npu], "PU:" + PUr[npu] + " Vertex position ; #Delta X [cm]; # entries", 100, 0., 0.01);
    m_plots[PUr[npu]+"Vtx_DY"] = new TH1D ("Vtx_DY_" + PUr[npu], "PU:" + PUr[npu] + " Vertex position ; #Delta Y [cm]; # entries", 100, 0., 0.01);
    m_plots[PUr[npu]+"Vtx_DZ"] = new TH1D ("Vtx_DZ_" + PUr[npu], "PU:" + PUr[npu] + " Vertex position ; #Delta Z [cm]; # entries", 240, 0., 24);
    m_plots[PUr[npu]+"Vtx_DXY"] = new TH1D ("Vtx_DXY_" + PUr[npu], "PU:" + PUr[npu] + " Vertex position ; #Delta XY [cm]; # entries", 100, 0., 0.01);
    // Delta Vertex (Reco-Gen)
    m_2Dplots[PUr[npu]+"Vtx_DXVsDY"] = new TH2D ("Vtx_DXVsDY_" + PUr[npu], "PU:" + PUr[npu] + " Vertex position ; #Delta X [cm]; #Delta Y [cm]", 100, 0., 0.01, 100, 0., 0.01);

  } // for(npu)  
  
  for (unsigned int nid=0; nid<IDName.size(); nid++){
    outFile->mkdir(sampleTag + "/KinIso_variables/" + IDName[nid]);
    for (unsigned int nre=0; nre<region.size(); nre++){
      outFile->mkdir(sampleTag + "/KinIso_variables/" + IDName[nid] + "/" + region[nre]);
      for (unsigned int npu=0; npu<PUr.size(); npu++){
	outFile->mkdir(sampleTag + "/KinIso_variables/" + IDName[nid] + "/" + region[nre] + "/" + PUr[npu]);
	outFile->cd   (sampleTag + "/KinIso_variables/" + IDName[nid] + "/" + region[nre] + "/" + PUr[npu]);
	for (unsigned int npT=0; npT<pTr.size(); npT++){
	  outFile->mkdir(sampleTag + "/KinIso_variables/" + IDName[nid] + "/" + region[nre] + "/" + PUr[npu] + "/" + pTr[npT]);
	  outFile->cd   (sampleTag + "/KinIso_variables/" + IDName[nid] + "/" + region[nre] + "/" + PUr[npu] + "/" + pTr[npT]);
	  
	  TString IDReg      = IDName[nid] + region[nre] + PUr[npu] + pTr[npT];
	  TString IDRegName  = IDName[nid] + "_" + region[nre] + "_" + PUr[npu] + "_" + pTr[npT];
	  TString IDRegTitle = "ID:" + IDName[nid] + " - REGION:" + region[nre] + " - PU:" + PUr[npu] + " - pT:" + pTr[npT];
	  
	  m_prof[IDReg+"dXYVsPt"]   = new TProfile("dXYVsPt_" + IDRegName, IDRegTitle + " #Delta XY Vs Pt ; #Delta XY;  p_{T} [GeV]", 200, 0, 2. , 0., 100.);
	  m_prof[IDReg+"dZVsPt"]    = new TProfile("dZVsPt_"  + IDRegName, IDRegTitle + " #Delta Z Vs Pt ; #Delta Z; p_{T} [GeV]" , 200, 0, 2. , 0., 100.);
	  m_prof[IDReg+"PtVsdXY"]   = new TProfile("PtVsdXY_" + IDRegName, IDRegTitle + " Pt Vs #Delta XY ; p_{T} [GeV]; #Delta XY [cm]", 100,0., 100, 0., 2.);
	  m_prof[IDReg+"PtVsdZ"]    = new TProfile("PtVsdZ_"  + IDRegName, IDRegTitle + " Pt Vs #Delta Z ; p_{T} [GeV]; #Delta Z [cm]", 100,0., 100, 0., 2.);
	  m_prof[IDReg+"PtVsVtx"]   = new TProfile("PtVsVtx_" + IDRegName, IDRegTitle + " Pt Vs #Vtx ; p_{T} [GeV]; N. Vtx", 100,0., 100, 0., 60.);

	  m_2Dplots[IDReg+"PtVspdgID"]   = new TH2D("PtVspdgID_" + IDRegName, IDRegTitle + " Pt Vs PDG ID ; p_{T} [GeV]; Particle ID", 100,0., 100, 1000, -1., 999.);

	  m_plots[IDReg+"NumGenJets"] = new TH1D ("NGenJets_" + IDRegName, IDRegTitle + " Number of GEN Jets ; GEN-Jets Multiplicity; # entries", 7, -0.5, 6.5);
	

	  for (unsigned int nISO=0; nISO<IsoName.size(); nISO++){
	    outFile->mkdir(sampleTag + "/KinIso_variables/" + IDName[nid] + "/" + region[nre] + "/" + PUr[npu] + "/" + pTr[npT] + "/" + IsoName[nISO]);
	    outFile->cd   (sampleTag + "/KinIso_variables/" + IDName[nid] + "/" + region[nre] + "/" + PUr[npu] + "/" + pTr[npT] + "/" + IsoName[nISO]);
	    
	    IDReg      = IDName[nid] + region[nre] + PUr[npu] + pTr[npT]  + IsoName[nISO];
	    IDRegName  = IDName[nid] + "_" + region[nre] + "_" + PUr[npu] + "_" + pTr[npT] + "_" + IsoName[nISO];
	    IDRegTitle = "ID:" + IDName[nid] + " - REGION:" + region[nre] + " - PU:" + PUr[npu] + " - pT:" + pTr[npT] + " - ISO:" + IsoName[nISO];
	    
	    m_plots[IDReg+"ChHadIso"] = new TH1D ("ChHadIso_" + IDRegName,  IDRegTitle + " Charged Hadron Isolation ; Charged Had. Iso.; # entries", 40, 0., 10);
	    m_plots[IDReg+"NeHadIso"] = new TH1D ("NeHadIso_" + IDRegName,  IDRegTitle + " Neutral Hadron Isolation ; Neutral Had. Iso.; # entries", 40, 0., 10);
	    m_plots[IDReg+"PhIso"]    = new TH1D ("PhIso_" + IDRegName,     IDRegTitle + " Photon Isolation ; Photon Iso.; # entries", 40, 0., 60);
	    m_plots[IDReg+"PUIso"]    = new TH1D ("PUIso_" + IDRegName,     IDRegTitle + " PU Isolation ; PU Iso.; # entries", 40, 0., 60);

	    m_prof[IDReg+"VtxVsChHadIso"] = new TProfile("VtxVsChHadIso_" + IDRegName, IDRegTitle + " Vtx Vs Charged Hadron Isolation ; N. Vtx; Charged Had. Iso.", 60, 0., 60., 0., 10.);
	    m_prof[IDReg+"VtxVsNeHadIso"] = new TProfile("VtxVsNeHadIso_" + IDRegName, IDRegTitle + " Vtx Vs Neutral Hadron Isolation ; N. Vtx; Neutral Had. Iso.", 60, 0., 60., 0., 10.); 
	    m_prof[IDReg+"VtxVsPhIso"]    = new TProfile("VtxVsPhIso_"    + IDRegName, IDRegTitle + " Vtx Vs Photon Isolation ; N. Vtx; Photon Iso.", 60, 0., 60., 0., 60.);
	    m_prof[IDReg+"VtxVsPUIso"]    = new TProfile("VtxVsPUIso_"    + IDRegName, IDRegTitle + " Vtx Vs PU Isolation ; N. Vtx; PU Iso.", 60, 0., 60., 0., 60.);
	      
	    m_plots[IDReg+"TotalIso"]     = new TH1D ("TotalIso_" + IDRegName,  IDRegTitle + " Total Isolation ; Total Iso.; # entries", 100, 0., 200);
	    m_prof[IDReg+"VtxVsTotalIso"]   = new TProfile("VtxVsTotalIso_" + IDRegName, IDRegTitle + " Vtx Vs Total Isolation ; N. Vtx; Total Iso.", 60, 0., 60., 0., 200.);
	    m_prof[IDReg+"PtVsTotalIso"]    = new TProfile("PtVsTotalIso_"  + IDRegName, IDRegTitle + " Pt Vs Total Isolation ; p_{T} [GeV]; Total Iso.", 100,0., 100, 0., 200.);
	    m_prof[IDReg+"EtaVsTotalIso"]   = new TProfile("EtaVsTotalIso_" + IDRegName, IDRegTitle + " Eta Vs Total Isolation ; #eta; Total Iso.", 48, -2.4, 2.4, 0., 200.);
	    m_prof[IDReg+"NumGenJetsVsTotalIso"] = new TProfile("NGenJetsVsTotalIso_" + IDRegName, IDRegTitle + " Number of GenJets Vs Total Isolation ; GEN-Jets Multiplicity; Total Iso.", 7, -0.5, 6.5, 0., 200.);

	    m_plots[IDReg+"TotalIsoNoDB"] = new TH1D ("TotalIsoNoDB_" + IDRegName,  IDRegTitle + " Total Isolation w/o #Delta#beta ; Total Iso. (No #Delta#beta); # entries", 100, 0., 200);
	    m_prof[IDReg+"VtxVsTotalIsoNoDB"]   = new TProfile("VtxVsTotalIsoNoDB_" + IDRegName, IDRegTitle + " Vtx Vs Total Isolation  w/o #Delta#beta ; N. Vtx; Total Iso. (No #Delta#beta)", 60, 0., 60., 0., 200.);
	    m_prof[IDReg+"PtVsTotalIsoNoDB"]    = new TProfile("PtVsTotalIsoNoDB_"  + IDRegName, IDRegTitle + " Pt Vs Total Isolation  w/o #Delta#beta ; p_{T} [GeV]; Total Iso. (No #Delta#beta)", 100,0., 100, 0., 200.);
	    m_prof[IDReg+"EtaVsTotalIsoNoDB"]   = new TProfile("EtaVsTotalIsoNoDB_" + IDRegName, IDRegTitle + " Eta Vs Total Isolation  w/o #Delta#beta ; #eta; Total Iso. (No #Delta#beta)", 48, -2.4, 2.4, 0., 200.);
	    m_prof[IDReg+"NumGenJetsVsTotalIsoNoDB"]   = new TProfile("NGenJetsVsTotalIsoNoDB_" + IDRegName, IDRegTitle + " Number of GenJets Vs Total Isolation ; GEN-Jets Multiplicity; Total Iso. (No #Delta#beta)", 7, -0.5, 6.5, 0., 200.);

	    m_plots[IDReg+"RelIso"]       = new TH1D ("RelIso_" + IDRegName,     IDRegTitle + " Relative Isolation ; Relative Iso.; # entries", 100, 0., 1.);
	    m_prof[IDReg+"VtxVsRelIso"]   = new TProfile("VtxVsRelIso_" + IDRegName, IDRegTitle + " Vtx Vs Relative Isolation ; N. Vtx; Relative Iso.", 60, 0., 60., 0., 1.);
	    m_prof[IDReg+"PtVsRelIso"]    = new TProfile("PtVsRelIso_"  + IDRegName, IDRegTitle + " Pt Vs Relative Isolation ; p_{T} [GeV]; Relative Iso.", 100,0., 100, 0., 1.);
	    m_prof[IDReg+"EtaVsRelIso"]   = new TProfile("EtaVsRelIso_" + IDRegName, IDRegTitle + " Eta Vs Relative Isolation ; #eta; Relative Iso.", 48, -2.4, 2.4, 0., 1.);
	    m_prof[IDReg+"dXYVsRelIso"]   = new TProfile("dXYVsRelIso_" + IDRegName, IDRegTitle + " #Delta XY Vs Relative Isolation ; #Delta XY; Relative Iso.", 200, 0, 2. , 0., 1.);
	    m_prof[IDReg+"dZVsRelIso"]    = new TProfile("dZVsRelIso_"  + IDRegName, IDRegTitle + " #Delta Z Vs Relative Isolation ; #Delta Z; Relative Iso.", 200, 0, 2. , 0., 1.);
	    m_prof[IDReg+"NumGenJetsVsRelIso"]   = new TProfile("NGenJetsVsRelIso_" + IDRegName, IDRegTitle + " Number of GenJets Vs Relative Isolation ; GEN-Jets Multiplicity; Relative Iso.", 7, -0.5, 6.5, 0., 1.);

	    m_plots[IDReg+"RelIsoNoDB"]   = new TH1D ("RelIsoNoDB_" + IDRegName, IDRegTitle + " Relative Isolation w/o #Delta#beta ; Relative Iso. (No #Delta#beta); # entries", 100, 0., 2.);
	    m_prof[IDReg+"VtxVsRelIsoNoDB"]   = new TProfile("VtxVsRelIsoNoDB_" + IDRegName, IDRegTitle + " Vtx Vs Relative Isolation  w/o #Delta#beta ; N. Vtx; Relative Iso. (No #Delta#beta)", 60, 0., 60., 0., 1.);
	    m_prof[IDReg+"PtVsRelIsoNoDB"]    = new TProfile("PtVsRelIsoNoDB_"  + IDRegName, IDRegTitle + " Pt Vs Relative Isolation  w/o #Delta#beta ; p_{T} [GeV]; Relative Iso. (No #Delta#beta)", 100,0., 100, 0., 1.);
	    m_prof[IDReg+"EtaVsRelIsoNoDB"]   = new TProfile("EtaVsRelIsoNoDB_" + IDRegName, IDRegTitle + " Eta Vs Relative Isolation  w/o #Delta#beta ; #eta; Relative Iso. (No #Delta#beta)", 48, -2.4, 2.4, 0., 1.);
	    m_prof[IDReg+"dXYVsRelIsoNoDB"]   = new TProfile("dXYVsRelIsoNoDB_" + IDRegName, IDRegTitle + " #Delta XY Vs Relative Isolation w/o #Delta#beta ; #Delta XY; Relative Iso.", 200, 0, 2. , 0., 1.);
	    m_prof[IDReg+"dZVsRelIsoNoDB"]    = new TProfile("dZVsRelIsoNoDB_"  + IDRegName, IDRegTitle + " #Delta Z Vs Relative Isolation w/o #Delta#beta ; #Delta Z; Relative Iso.", 200, 0, 2. , 0., 1.);
	    m_prof[IDReg+"NumGenJetsVsRelIsoNoDB"]   = new TProfile("NGenJetsVsRelIsoNoDB_" + IDRegName, IDRegTitle + " Number of GenJets Vs Relative Isolation w/o #Delta#beta ; GEN-Jets Multiplicity; Relative Iso.", 7, -0.5, 6.5, 0., 1.);
	    
	    
	    m_plots[IDReg+"RelIsoCut"] = new TH1D ("RelIsoCut_" + IDRegName, IDRegTitle + " Relative Isolation Cut ; Relative Iso. Cut; # entries", 101, 0., 1.01);
	    
	    for (unsigned int nISOcut=0; nISOcut<isowp[nISO].size(); nISOcut++){
	      
	      IDReg      = IDName[nid] + region[nre] + PUr[npu] + pTr[npT]  + IsoName[nISO] + isowp[nISO][nISOcut];
	      IDRegName  = IDName[nid] + "_" + region[nre] + "_" + PUr[npu] + "_" + pTr[npT] + "_" + IsoName[nISO] + "_" + isowp[nISO][nISOcut];
	      IDRegTitle = "ID:" + IDName[nid] + " - REGION:" + region[nre] + " - PU:" + PUr[npu] + " - pT:" + pTr[npT] + " - ISO:" + IsoName[nISO] + " - ISOWP:" + isowp[nISO][nISOcut];
	      
	      m_plots[IDReg+"NMuons"] = new TH1D ("NMuons_" + IDRegName, IDRegTitle + " Number of Muons per Evt ; N. Muons; # entries", 4, 0., 4.);
	      
	      m_plots[IDReg+"Vtx"] = new TH1D ("Vtx_" + IDRegName,  IDRegTitle + " Vtx ; N. Vtx; # entries", 60, 0., 60.);
	      m_plots[IDReg+"Pt"]  = new TH1D ("Pt_"  + IDRegName,  IDRegTitle + " Pt ; p_{T} [GeV]; # entries", 100,0., 100);
	      m_plots[IDReg+"Eta"] = new TH1D ("Eta_" + IDRegName,  IDRegTitle + " Eta ; #eta; # entries", 48, -2.4, 2.4);
	      m_plots[IDReg+"dXY"] = new TH1D ("dXY_" + IDRegName,  IDRegTitle + " #Delta XY ; #Delta XY; # entries", 200, 0, 2.);
	      m_plots[IDReg+"dZ"]  = new TH1D ("dZ_"  + IDRegName,  IDRegTitle + " #Delta Z ; #Delta Z; # entries", 200, 0, 2.);
	      
	      m_plots[IDReg+"NumGenJets"] = new TH1D ("NGenJets_" + IDRegName, IDRegTitle + " Number of GEN Jets ; GEN-Jets Multiplicity; # entries", 7, -0.5, 6.5);
	      	      
	    }// for(nISOcut)
	  }// for(nISO)
	}// for(npT)
      }// for(npu)
    } // for(nre)
  } // for(nid)

  std::cout << "All histograms have been created." << std::endl;
  
} // void::book 

void muon_pog::Plotter::fillGen(const std::vector<muon_pog::GenParticle> & genpars, const muon_pog::Event & ev){

  // N. Muons
  int NmuonFromZ = 0;
  int NmuonOther = 0;
  
  float etaBarrel = 1.2;  
  float pTeta = 20.0;
  
  // Number of Events 
  m_2Dyields["GenMuon"]->Fill(0.,0.);

  for (auto & genpar : genpars){

    if(abs(genpar.pdgId) == 13 &&
       genpar.pt          > m_tnpConfig.probe_minPt &&
       fabs(genpar.eta)   < 2.4 ){
      
      bool genZmuon     = false;
      
      // Check this with flags!!!!
      if (hasMother(genpar, 13)) continue;
      
      // Z_id = 23
      genZmuon = hasMother(genpar, 23);
      
      for (unsigned int nre = 0; nre<4; nre++){
	
	bool etaRegion = false;
	if (nre == 0) etaRegion = true;
	if (nre == 1 && std::abs(genpar.eta) <  etaBarrel) etaRegion = true;
	if (nre == 2 && std::abs(genpar.eta) >= etaBarrel) etaRegion = true;
	if (nre == 3 && (std::abs(genpar.eta) > 0.9 && std::abs(genpar.eta) < 1.3)) etaRegion = true;
	
	if (!etaRegion)     continue;
	
	if (genZmuon) m_2Dyields["GenMuon"]->Fill(1,nre); // From Z
	else          m_2Dyields["GenMuon"]->Fill(2,nre); 
	m_2Dyields["GenMuon"]->Fill(3,nre); // Number of GenMuons
	
      }// for(nre)
      
      if (genZmuon) NmuonFromZ++;
      else NmuonOther++;

    }//if(GenMuon)
    
  }// for(genpar)  

  m_plots["NumberOfGenMuonsFromZ"]->Fill(NmuonFromZ);
  m_plots["NumberOfGenMuonsOther"]->Fill(NmuonOther);
  
}
 
void muon_pog::Plotter::fill(const std::vector<muon_pog::Muon> & muons,
			     const muon_pog::HLT & hlt, const muon_pog::Event & ev, float weight)
{

  // Total Number of Events 
  m_2Dyields["RecoMuon"]->Fill(IDName.size()+2.5, 0., weight);

  // Primary Vertices: Reco and from GenParticles
  std::vector<bool> IsPURegime;
  for (int bnpu = 0; bnpu<PUr.size(); bnpu++) IsPURegime.push_back( false );
  IsPURegime[FullPU]   = true;
  IsPURegime[LowPU]    = (ev.nVtx < 10);
  IsPURegime[MediumPU] = (ev.nVtx >= 10 && ev.nVtx < 35);
  IsPURegime[HighPU]   = (ev.nVtx >= 35);

  TVector3 RecoVtx;
  RecoVtx.SetXYZ(ev.primaryVertex[0],ev.primaryVertex[1],ev.primaryVertex[2]);

  TVector3 GenVtx;
  bool FPrPar = true;
  for (auto & genpar : ev.genParticles){
    if(genpar.flags.at(7) == 1 && FPrPar){ // Flag 7: IsHardProcess
      GenVtx.SetXYZ(genpar.vx,genpar.vy,genpar.vz);      
      FPrPar = false;
    } // if(flag == 7)
  } // for(genpar)
  
  for (unsigned int npu = 0; npu<PUr.size(); npu++){
    
    if (!IsPURegime[npu])  continue;

    m_plots[PUr[npu]+"RecoVtx_X"] -> Fill(RecoVtx.X(), weight);
    m_plots[PUr[npu]+"RecoVtx_Y"] -> Fill(RecoVtx.Y(), weight);
    m_plots[PUr[npu]+"RecoVtx_Z"] -> Fill(RecoVtx.Z(), weight);
    m_plots[PUr[npu]+"RecoVtx_XY"] -> Fill(RecoVtx.XYvector().Mod(), weight);
    
    m_plots[PUr[npu]+"GenVtx_X"] -> Fill(GenVtx.X(), weight);
    m_plots[PUr[npu]+"GenVtx_Y"] -> Fill(GenVtx.Y(), weight);
    m_plots[PUr[npu]+"GenVtx_Z"] -> Fill(GenVtx.Z(), weight);
    m_plots[PUr[npu]+"GenVtx_XY"] -> Fill(GenVtx.XYvector().Mod(), weight);
    
    m_plots[PUr[npu]+"Vtx_DX"] -> Fill(std::abs(RecoVtx.X() - GenVtx.X()), weight);
    m_plots[PUr[npu]+"Vtx_DY"] -> Fill(std::abs(RecoVtx.Y() - GenVtx.Y()), weight);
    m_plots[PUr[npu]+"Vtx_DZ"] -> Fill(std::abs(RecoVtx.Z() - GenVtx.Z()), weight);
    m_plots[PUr[npu]+"Vtx_DXY"] -> Fill(std::abs(RecoVtx.XYvector().Mod() - GenVtx.XYvector().Mod()), weight);
    
    m_2Dplots[PUr[npu]+"Vtx_DXVsDY"] -> Fill(std::abs(RecoVtx.X() - GenVtx.X()), std::abs(RecoVtx.Y() - GenVtx.Y()), weight);

  } // for(npu)

  // Number of GEN-Jets
  int NGenJets = ev.genjets.size();

  // Number of muons
  std::map<TString, int>  number;
  // -- Initialization
  number["muons"] = 0;
  for (unsigned int mid = 0; mid<IDName.size(); mid++)
    for (unsigned int nre = 0; nre<region.size(); nre++)
      for (unsigned int npu = 0; npu<PUr.size(); npu++)
  	for (unsigned int npT = 0; npT<pTr.size(); npT++)	  
  	  for(unsigned int iISO=0; iISO < IsoName.size(); iISO++)
  	    for(unsigned int isocat=0; isocat<isowp[iISO].size(); isocat++){
  	      TString iname = IDName[mid]+region[nre]+PUr[npu]+pTr[npT]+IsoName[iISO]+isowp[iISO][isocat]; 
  	      number[iname+"muons"]=0;
  	    }
  
  // Muon General Cuts
  float etaBarrel = 1.2;  
  float pTeta = 20.0;
  
  float DBfactor = m_tnpConfig.iso_DeltaBeta;

  // Isolation Cuts from config file
  std::vector< std::vector<float> > iISOwpFromUser; 
  iISOwpFromUser.push_back({9999, m_tnpConfig.iso_LooseWP,     m_tnpConfig.iso_TightWP});
  iISOwpFromUser.push_back({m_tnpConfig.isoPUPPI_LooseWP,      m_tnpConfig.isoPUPPI_TightWP});
  iISOwpFromUser.push_back({m_tnpConfig.isoPUPPILep_LooseWP,   m_tnpConfig.isoPUPPILep_TightWP});
  iISOwpFromUser.push_back({m_tnpConfig.isoPUPPINoLep_LooseWP, m_tnpConfig.isoPUPPINoLep_TightWP});
  iISOwpFromUser.push_back({m_tnpConfig.Miniiso_LooseWP,       m_tnpConfig.Miniiso_TightWP});
  
  //---------------------------------------
  // PV comparison wrt the HP position  ---
  //---------------------------------------
  // if( std::abs(RecoVtx.XYvector().Mod() - GenVtx.XYvector().Mod()) < 0.001 &&
  //     std::abs(RecoVtx.Z() - GenVtx.Z()) < 5.0 ){

  for (auto & muon : muons){
 
    // Muons/Evt
    number["muons"]++;

    // Total Number of Muons
    m_2Dyields["RecoMuon"]->Fill(IDName.size()+1.5, 0., weight);
    
    // General Probe Muons	      
    if(muon.pt > m_tnpConfig.probe_minPt && 
       fabs(muon.eta) < 2.4){

      // Number of Muons Passing (pT,eta) cuts
      m_2Dyields["RecoMuon"]->Fill(IDName.size()+0.5, 0., weight);
      
      float dXY = std::abs(muon.dxyBest); 
      float dZ =  std::abs(muon.dzBest);

      bool IsSIGN = false;
      
      // ------- OLD ------- 
      // Check first if muon comes from Z
      // Check "hasGenMatch" and "hasNoGenMatch"
      // IsSIGN = muon_pog::hasGenMatch(muon, ev.genParticles, m_tnpConfig.gen_DrCut, 23);

      // ------- NEW ------- 
      // Check the SIMHit info
      IsSIGN = muon.IsMatchedPrimaryMuon;

      bool FillMuon = false;
      TString sampleTag = m_sampleConfig.sampleName;
      // Only prompt muons from signal
      if(sampleTag.Contains("SIGNAL") && IsSIGN) FillMuon = true;
      // All events from QCD samples
      if(sampleTag.Contains("BACKGROUND")) FillMuon = true;

      if (!FillMuon) continue;
      
      // Muon mother 
      // ------- OLD ------- 
      // int GenMatchMotherID = MotherGenMatch(muon, ev.genParticles, m_tnpConfig.gen_DrCut); // OLD Approach
      // ------- NEW ------- 
      int GenMatchMotherID = muon.SimmotherPdgId; // New Approach: SIMHit Matching
      

      std::vector<bool> IsMuonID;
      for (int bmid = 0; bmid<IDName.size(); bmid++) IsMuonID.push_back( false );
      // New Selectors From MiniAOD
      IsMuonID[BinType::GLB]           = muon.isGlobal ; 
      IsMuonID[BinType::TRK]           = muon.isTracker ; 
      IsMuonID[BinType::LOOSE]         = muon.Sel_CutBasedIdLoose ; 
      IsMuonID[BinType::MEDIUM]        = muon.Sel_CutBasedIdMedium ; 
      IsMuonID[BinType::TIGHT]         = muon.Sel_CutBasedIdTight ; 
      IsMuonID[BinType::SOFT]          = muon.Sel_SoftCutBasedId ; 
      IsMuonID[BinType::HIGHpT]        = muon.Sel_CutBasedIdGlobalHighPt ; 
      IsMuonID[BinType::MEDIUMPrompt]  = muon.Sel_CutBasedIdMediumPrompt ; 
      IsMuonID[BinType::MVALOOSE]      = muon.Sel_MvaLoose ; 
      IsMuonID[BinType::MVAMEDIUM]     = muon.Sel_MvaMedium ; 
      IsMuonID[BinType::MVATIGHT]      = muon.Sel_MvaTight ; 

      std::vector<bool> IsEtaRegion;
      for (int bnre = 0; bnre<region.size(); bnre++) IsEtaRegion.push_back( false );
      IsEtaRegion[BinRegion::Full]    = true;
      IsEtaRegion[BinRegion::Barrel]  = (std::abs(muon.eta) < etaBarrel);
      IsEtaRegion[BinRegion::Endcap]  = (std::abs(muon.eta) >= etaBarrel);
      IsEtaRegion[BinRegion::Overlap] = ((std::abs(muon.eta) > 0.9 && std::abs(muon.eta) < 1.3));

      std::vector<bool> IspTRegime;
      for (int bnre = 0; bnre<pTr.size(); bnre++) IspTRegime.push_back( false );
      IspTRegime[BinpTr::FullpT]   = true;
      IspTRegime[BinpTr::LowpT]    = (muon.pt <  20.);
      IspTRegime[BinpTr::MediumpT] = (muon.pt >= 20. && muon.pt< 50.);
      IspTRegime[BinpTr::HighpT]  = (muon.pt >= 50. );

      for (unsigned int mid = 0; mid<IDName.size(); mid++){

	if (!IsMuonID[mid])    continue;

	for (unsigned int nre = 0; nre<region.size(); nre++){
	  
	  if (!IsEtaRegion[nre]) continue;
	  
	  m_2Dyields["RecoMuon"]->Fill(mid, nre, weight);	  
	  
	  for (unsigned int npu = 0; npu<PUr.size(); npu++){
	    
	    if (!IsPURegime[npu])  continue;

	    for (unsigned int npT = 0; npT<pTr.size(); npT++){
	      
	      if (!IspTRegime[npT])  continue;
	      
	      TString HNameRef = IDName[mid]+region[nre]+PUr[npu]+pTr[npT];
	      
	      m_prof[HNameRef+"dXYVsPt"] ->Fill(dXY, muon.pt, weight);
	      m_prof[HNameRef+"dZVsPt"]  ->Fill(dZ,  muon.pt, weight);
	      
	      m_prof[HNameRef+"PtVsdXY"] ->Fill(muon.pt, dXY, weight);
	      m_prof[HNameRef+"PtVsdZ"]  ->Fill(muon.pt, dZ,  weight);
	      
	      m_prof[HNameRef+"PtVsVtx"] ->Fill(muon.pt, ev.nVtx, weight);
	    
	      m_2Dplots[HNameRef+"PtVspdgID"]->Fill(muon.pt, GenMatchMotherID, weight);

	      m_plots[HNameRef+"NumGenJets"] ->Fill(NGenJets, weight);

	      std::vector<float> ChHadIso, NeHadIso, PhIso, ChPUIso;
	      std::vector<float> TotalIso, TotalIsoNoDB, RelIso, RelIsoNoDB;
	      // Isolation (R=0.4)
	      ChHadIso.push_back(muon.chargedHadronIso);
	      NeHadIso.push_back(muon.neutralHadronIso);
	      PhIso.push_back   (muon.photonIso);
	      ChPUIso.push_back (muon.chargedHadronIsoPU);	      	      
	      // PUPPI 
	      ChHadIso.push_back(muon.PUPPIIsoCH);
	      NeHadIso.push_back(muon.PUPPIIsoNH);
	      PhIso.push_back   (muon.PUPPIIsoPH);
	      ChPUIso.push_back (0.0);
	      // PUPPI Lep
	      ChHadIso.push_back(0.0);
	      NeHadIso.push_back(0.0);
	      PhIso.push_back   (0.0);
	      ChPUIso.push_back (0.0);
	      // PUPPI NoLep
	      ChHadIso.push_back(0.0);
	      NeHadIso.push_back(0.0);
	      PhIso.push_back   (0.0);
	      ChPUIso.push_back (0.0);
	      // Mini-Isolation
	      ChHadIso.push_back(muon.MiniIsoCH);
	      NeHadIso.push_back(muon.MiniIsoNH);
	      PhIso.push_back   (muon.MiniIsoPH);
	      ChPUIso.push_back (muon.MiniIsoPU);
	      
	      // --Total Isolations
	      for(unsigned int iISO=0; iISO < IsoName.size(); iISO++){
		float iIso     = ChHadIso[iISO] + std::max(0., NeHadIso[iISO] + PhIso[iISO] - 1.0*DBfactor*ChPUIso[iISO]);
		float iIsoNoDB = ChHadIso[iISO] + NeHadIso[iISO] + PhIso[iISO];

		if(iISO == BinIsoName::ISOPUPPILep){
		  iIso     = muon.PUPPILepIso;
		  iIsoNoDB = muon.PUPPILepIso;
		}
		
		if(iISO == BinIsoName::ISOPUPPINoLep){
		  iIso     = muon.PUPPINoLepIso;
		  iIsoNoDB = muon.PUPPINoLepIso;
		}

		TotalIso.push_back(iIso);
		if(iISO == BinIsoName::ISOPF || iISO == BinIsoName::MiniISO) 
		  RelIso.push_back(iIso/muon.pt); 
		else RelIso.push_back(iIso); 
		
		TotalIsoNoDB.push_back(iIsoNoDB);
		if(iISO == BinIsoName::ISOPF || iISO == BinIsoName::MiniISO) 
		  RelIsoNoDB.push_back(iIsoNoDB/muon.pt);
		else RelIsoNoDB.push_back(iIsoNoDB);

	      }
	      
	      for(unsigned int iISO=0; iISO < IsoName.size(); iISO++){
		
		//std::cout << iISO << " for "<< RelIsoNoDB.size() << " value=" << RelIsoNoDB[iISO] << " relIsoNoDB and " << RelIso.size() << " for RelIso value=" << RelIso[iISO] << std::endl;
		
		HNameRef = IDName[mid]+region[nre]+PUr[npu]+pTr[npT]+IsoName[iISO];

		// Isolation Components
		m_plots[HNameRef+"ChHadIso"]->Fill(ChHadIso[iISO],weight);
		m_plots[HNameRef+"NeHadIso"]->Fill(NeHadIso[iISO],weight);
		m_plots[HNameRef+"PhIso"]   ->Fill(PhIso[iISO],   weight);
		m_plots[HNameRef+"PUIso"]   ->Fill(ChPUIso[iISO], weight);
		// Isolation Components Vs Vtx
		m_prof[HNameRef+"VtxVsChHadIso"]->Fill(ev.nVtx, ChHadIso[iISO], weight);
		m_prof[HNameRef+"VtxVsNeHadIso"]->Fill(ev.nVtx, NeHadIso[iISO], weight);
		m_prof[HNameRef+"VtxVsPhIso"]   ->Fill(ev.nVtx, PhIso[iISO],    weight);
		m_prof[HNameRef+"VtxVsPUIso"]   ->Fill(ev.nVtx, ChPUIso[iISO],  weight);
		// ---------------------------------------------------------------------------------------
		// Total Isolation
		m_plots[HNameRef+"TotalIso"]           ->Fill(TotalIso[iISO],          weight);
		m_prof[HNameRef+"VtxVsTotalIso"]       ->Fill(ev.nVtx, TotalIso[iISO], weight);
		m_prof[HNameRef+"PtVsTotalIso"]        ->Fill(muon.pt, TotalIso[iISO], weight);
		m_prof[HNameRef+"NumGenJetsVsTotalIso"]->Fill(NGenJets,TotalIso[iISO], weight);
		// Total Isolation w/o DB corrections
		m_plots[HNameRef+"TotalIsoNoDB"]           ->Fill(TotalIsoNoDB[iISO],          weight);
		m_prof[HNameRef+"VtxVsTotalIsoNoDB"]       ->Fill(ev.nVtx, TotalIsoNoDB[iISO], weight);
		m_prof[HNameRef+"PtVsTotalIsoNoDB"]        ->Fill(muon.pt, TotalIsoNoDB[iISO], weight);
		m_prof[HNameRef+"NumGenJetsVsTotalIsoNoDB"]->Fill(NGenJets,TotalIsoNoDB[iISO], weight);
		// Relative Isolation
		m_plots[HNameRef+"RelIso"]           ->Fill(RelIso[iISO],          weight);		
		m_prof[HNameRef+"VtxVsRelIso"]       ->Fill(ev.nVtx, RelIso[iISO], weight);
		m_prof[HNameRef+"PtVsRelIso"]        ->Fill(muon.pt, RelIso[iISO], weight);
		m_prof[HNameRef+"NumGenJetsVsRelIso"]->Fill(NGenJets, RelIso[iISO],weight);
		m_prof[HNameRef+"dXYVsRelIso"]       ->Fill(dXY, RelIso[iISO],     weight);
		m_prof[HNameRef+"dZVsRelIso"]        ->Fill(dZ,  RelIso[iISO],     weight);
		// Relative Isolation w/o DB corrections
		m_plots[HNameRef+"RelIsoNoDB"]           ->Fill(RelIsoNoDB[iISO],          weight);
		m_prof[HNameRef+"VtxVsRelIsoNoDB"]       ->Fill(ev.nVtx, RelIsoNoDB[iISO], weight);
		m_prof[HNameRef+"PtVsRelIsoNoDB"]        ->Fill(muon.pt, RelIsoNoDB[iISO], weight);
		m_prof[HNameRef+"NumGenJetsVsRelIsoNoDB"]->Fill(NGenJets, RelIsoNoDB[iISO],weight);
		m_prof[HNameRef+"dXYVsRelIsoNoDB"]       ->Fill(dXY, RelIsoNoDB[iISO],     weight);
		m_prof[HNameRef+"dZVsRelIsoNoDB"]        ->Fill(dZ, RelIsoNoDB[iISO],      weight);

		if(muon.pt > pTeta){
		  m_prof[HNameRef+"EtaVsTotalIso"]    ->Fill(muon.eta, TotalIso[iISO],    weight);
		  m_prof[HNameRef+"EtaVsTotalIsoNoDB"]->Fill(muon.eta, TotalIsoNoDB[iISO],weight);
		  m_prof[HNameRef+"EtaVsRelIso"]      ->Fill(muon.eta, RelIso[iISO],       weight);
		  m_prof[HNameRef+"EtaVsRelIsoNoDB"]  ->Fill(muon.eta, RelIsoNoDB[iISO],  weight);
		}

		int vRelIsoCut = 1 + std::round(RelIso[iISO]*100.);
		for (int ibin=vRelIsoCut; ibin<=100; ibin++) m_plots[HNameRef+"RelIsoCut"]->AddBinContent(ibin, weight);
		m_plots[HNameRef+"RelIsoCut"]->AddBinContent(101, weight); // Total Number of Muons 
		
		// -- Isolation WorkingPoints
		std::vector<float>   iISOwp = iISOwpFromUser[iISO];
		std::vector<TString> nISOwp = isowp[iISO];
		
		for(int isocat=0;isocat < iISOwp.size();isocat++){
		  
		  if(RelIso[iISO] > iISOwp[isocat]) continue;
		  
		  HNameRef = IDName[mid]+region[nre]+PUr[npu]+pTr[npT]+IsoName[iISO]+nISOwp[isocat];
		
		  number[HNameRef+"muons"]++;
		
		  m_plots[HNameRef+"Vtx"]       ->Fill(ev.nVtx, weight);
		  m_plots[HNameRef+"Pt"]        ->Fill(muon.pt, weight);
		  m_plots[HNameRef+"NumGenJets"]->Fill(NGenJets,weight);
		  m_plots[HNameRef+"dXY"]       ->Fill(dXY,     weight);
		  m_plots[HNameRef+"dZ"]        ->Fill(dZ,      weight);
		  
		  if(muon.pt > pTeta) m_plots[HNameRef+"Eta"]              ->Fill(muon.eta, weight);
		
		} // for(IsoCategory)
		
	      } // for(Isolation)

	    } // for(pTRegime)

	  } // for(PURegime)
	
	} // for(etaRegion)
      
      } // for(muonID)
      
    } // if(muon pT eta)
  }// for(Muons)
  
  m_plots["NMuons"]->Fill(number["muons"],weight);
  
  for (unsigned int mid = 0; mid<IDName.size(); mid++)
    for (unsigned int nre = 0; nre<region.size(); nre++)
      for (unsigned int npu = 0; npu<PUr.size(); npu++)
  	for (unsigned int npT = 0; npT<pTr.size(); npT++)	  
  	  for(unsigned int iISO=0; iISO < IsoName.size(); iISO++)
  	    for(unsigned int isocat=0; isocat<isowp[iISO].size(); isocat++){
  	      TString iname = IDName[mid]+region[nre]+PUr[npu]+pTr[npT]+IsoName[iISO]+isowp[iISO][isocat]; 

  	      m_plots[iname+"NMuons"]->Fill(number[iname+"muons"],weight);
  	    }

  
  //}// IF(PV-HP)  
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

std::vector<TString> GetListOfFiles(TString FileName){
  std::vector<TString> ListOfSamples;
  std::ifstream InFile;
  InFile.open(FileName);
  if (!InFile){
    std::cout << "File " << FileName << " not found!" << std::endl;
    std::exit(0);
  }
  else{
    std::string tmpLine;
    while (std::getline(InFile,tmpLine)){
      TString TxtLine = tmpLine;
      if (!TxtLine.Contains("#") && TxtLine.Contains(".root")) ListOfSamples.push_back(TxtLine);
    }// while
  }
  
  return ListOfSamples;
}

//  LocalWords:  IsoTight
