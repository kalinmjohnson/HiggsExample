//
// Original Author:
//         Created:  Fri December 7, 2017 by  N.Z. Jomhari   
//                   with contributions from  A. Geiser
//                                            A. Anuar
//                   and pieces similar to previous analysis examples
// $Id$
// ..
//
// ***************************************************************************
//  Higgs-to-four-lepton analysis example at research level                  *
// ***************************************************************************
//                                                                           *
// Built upon the DEMO setup provided by CMS open data access team,          *
// expanded/upgraded to contain a pedagocigal analysis example for the       *
// approximate reproducttion of the Higgs to four lepton mass spectrum       *
// published in CMS-HIG-12-028                                               *
// Phys.Lett. B716 (2012) 30-61,  arXiv:1207.7235                            *
//                                                                           *
// This research level example is a strongly simplified reimplementation of  *
// parts of the original CMS Higgs to four lepton analysis                   *
//                                                                           * 
// The published reference plot which is being approximated in this example  *
// (in addition to many auxiliary plots) is                                  *
// https://inspirehep.net/record/1124338/files/H4l_mass_v3.png               *
// Other Higgs final states (e.g. Higgs to two photons), which were also     *
// part of the same CMS paper and strongly contributed to the Higgs          *
// discovery, are not covered by this example.                               *
//                                                                           *
// The example addresses users who feel they have at least some minimal      *
// understanding of the content of the CMS paper and of the meaning of this  *
// reference plot, which can be reached via (separate) educational exercises.*
// A Root ntuple containing the kinematic information about the event        *
// candidates, which could be used for educational purposes, is part of the  *
// Root output file.                                                         *
//                                                                           *
// The analysis code provided here recodes the spirit of the original        *
// analysis and (approximately) recodes many of the original cuts on         *
// original data objects, but does not provide the original analysis code    *
// itself. Also, for the sake of simplicity, it skips some of the more       *
// advanced analysis methods of the original paper, and does not use any     *
// corrections beyond those already implicit in the objects used, and no     *
// systematic uncertainties.                                                 * 
// Another reason why the published spectrum is only reproduced very         *
// approximately is that the data sets only partially overlap, and that the  *
// legacy software version and corresponding calibrations differ from those  *
// of the original paper.                                                    *                                                      
// Nevertheless, the example provides a qualitative insight about how the    *
// original result was obtained.                                             *
//                                                                           *
// In addition to the documented core results, the resulting  Root files     *
// also contain many undocumented plots which grew as a side product from    *
// setting up this example and earlier examples.                             *
// And it contains an ntuple with candidate four-vectors as mentioned above. *
// ***************************************************************************

// ***************************************************************************
// Analysis strategy                                                         *
// ***************************************************************************
//
// The analysis strategy is the following: Get the histograms for the 4mu    *
// and 2mu2e final states from the DoubleMu datasets and those for 4e final  *
// state from the DoubleElectron dataset. This avoids double counting due to *
// trigger overlaps.                                                         *
// The code itself is agnostic with respect to the input data set used, and  *
// the appropriate histograms have to selected at the subsequent root        *
// analysis step.                                                            *
// ***************************************************************************

// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <utility>

// user include files, general
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//------ EXTRA HEADER FILES--------------------//
#include "math.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"

// for Root histogramming
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include <TMath.h>

// for tracking information
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

// for vertex information 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// for muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

//for beamspot information
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//for electron informaton
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

// for particle flow information
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"



#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
using namespace std;


// class declaration
class HiggsDemoAnalyzer: public edm::EDAnalyzer {
public:
  explicit HiggsDemoAnalyzer(const edm::ParameterSet&);
  ~HiggsDemoAnalyzer();

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  bool providesGoodLumisection(const edm::Event& iEvent);

  // Declare Root histograms or tree
  // For a description of their content see below

  // TTree *t1;
  // TTree *t2;
  // TTree *t3;
  
  ofstream myfile4mu;
  ofstream myfile4e;
  ofstream myfile2mu2e;
  ofstream myfile;
  TTree *tree;

  // Muons
  const static int max_mu = 100;
  int value_mu_n = 0;
  float value_mu_pt[max_mu];
  float value_mu_eta[max_mu];
  float value_mu_phi[max_mu];
  float value_mu_mass[max_mu];
  int value_mu_charge[max_mu];
  float value_mu_pfreliso[max_mu];
  float value_mu_pfrelisoDB[max_mu];
  float value_mu_SIP3d[max_mu];
  float value_mu_dxy[max_mu];
  float value_mu_dz[max_mu];
  float value_mu_isGM[max_mu];
  const float mu_min_pt = 5;
  
  // Electrons
  const static int max_el = 100;
  int value_el_n = 0;
  float value_el_pt[max_el];
  float value_el_eta[max_el];
  float value_el_phi[max_el];
  float value_el_rawE[max_el];
  int value_el_charge[max_el];
  int value_el_misshits[max_el];
  float value_el_pfreliso[max_mu];
  float value_el_SIP3d[max_mu];
  float value_el_dxy[max_mu];
  float value_el_dz[max_mu];
  const float el_min_pt = 5;

  

  // Declare variables
  int  nGoodGlobalMuon, nGoodRecoMuon;

  double s1, s2, s3, s4, s;
  double dx,dy,dz, rap, pz;

  double mZ12, mZ34, mZ13, mZ24, mZ14, mZ23;
  double mass4mu, pt_4mu, eta_4mu, phi_4mu;
  double px4mu, py4mu, pz4mu, E4mu;
  double pt_mu1, pt_mu2, pt_mu3, pt_mu4;
  double eta_mu1, eta_mu2, eta_mu3, eta_mu4;
  double phi_mu1, phi_mu2, phi_mu3, phi_mu4;
  int cas_mu1, cas_mu2, cas_mu3, cas_mu4;

  double px_mu1, px_mu2, px_mu3, px_mu4;
  double py_mu1, py_mu2, py_mu3, py_mu4;
  double pz_mu1, pz_mu2, pz_mu3, pz_mu4;

  double E_mu1, E_mu2, E_mu3, E_mu4;

  double mZa, mZb;

  double sqm1, mZ, eZ12, eZ34, eZ13, eZ24, eZ14, eZ23; 
  double pxZ12, pxZ34, pxZ13, pxZ24,  pxZ14, pxZ23;
  double pyZ12, pyZ34, pyZ13, pyZ24, pyZ14, pyZ23;
  double pzZ12, pzZ34, pzZ13, pzZ24, pzZ14, pzZ23;
  double pZ12, pZ34, pZ13, pZ24, pZ14, pZ23;
  double pTZ12, pTZ34, pTZ13, pTZ24, pTZ14, pTZ23;

  double dZ12, dZ34, dZ13, dZ24, dZ14, dZ23;
  double dZc1, dZc2, dZc3;
  double eZa, pxZa, pyZa, pzZa, pTZa;
  double eZb, pxZb, pyZb, pzZb, pTZb;

  int nGoodElectron;
  double sqme;
  int misshits;

  double mass4e, pt_4e, eta_4e, phi_4e;
  double px4e, py4e, pz4e, E4e;
  double pt_e1, pt_e2, pt_e3, pt_e4;
  double eta_e1, eta_e2, eta_e3, eta_e4;
  double phi_e1, phi_e2, phi_e3, phi_e4;
  int cas_e1, cas_e2, cas_e3, cas_e4;

  double px_e1, px_e2, px_e3, px_e4;
  double py_e1, py_e2, py_e3, py_e4;
  double pz_e1, pz_e2, pz_e3, pz_e4;

  double E_e1, E_e2, E_e3, E_e4;

  double mass2mu2e, pt_2mu2e, eta_2mu2e, phi_2mu2e;
  double px2mu2e, py2mu2e, pz2mu2e, E2mu2e;
  double pt_2mu1, pt_2mu2, pt_2e1, pt_2e2;
  double eta_2mu1, eta_2mu2, eta_2e1, eta_2e2;
  double phi_2mu1, phi_2mu2, phi_2e1, phi_2e2;
  int cas_2mu1, cas_2mu2, cas_2e1, cas_2e2;

  double px_2mu1, px_2mu2, px_2e1, px_2e2;
  double py_2mu1, py_2mu2, py_2e1, py_2e2;
  double pz_2mu1, pz_2mu2, pz_2e1, pz_2e2;

  double E_2mu1, E_2mu2, E_2e1, E_2e2;

  double goodhit;
  double relPFIso_mu;
  double relPFIso_e;

  double IP3d_mu;
  double ErrIP3d_mu;
  double SIP3d_mu;
  double IP3d_e;
  double ErrIP3d_e;
  double SIP3d_e;

  int nRun,nEvt,nLumi;

  TLorentzVector p4Za, p4Zb, p4H;

  bool muonTree = false;
  bool elecTree = false;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

HiggsDemoAnalyzer::HiggsDemoAnalyzer(const edm::ParameterSet& iConfig) {

  // *****************************************************************
  // This is the main analysis routine
  // The goal is to approximately reproduce the Higgs-to-four-lepton    
  // mass spectrum from HIG-12-028
  // *****************************************************************

  // now do what ever initialization is needed
  edm::Service<TFileService> fs;

  // ************************************
  // book histograms and set axis labels
  // (called once for initialization)
  // ************************************
  tree = fs->make<TTree>("Events", "Events");
  

  // Muons
  tree->Branch("nMuon", &value_mu_n, "nMuon/i");
  tree->Branch("Muon_pt", value_mu_pt, "Muon_pt[nMuon]/F");
  tree->Branch("Muon_eta", value_mu_eta, "Muon_eta[nMuon]/F");
  tree->Branch("Muon_phi", value_mu_phi, "Muon_phi[nMuon]/F");
  tree->Branch("Muon_mass", value_mu_mass, "Muon_mass[nMuon]/F");
  tree->Branch("Muon_charge", value_mu_charge, "Muon_charge[nMuon]/I");
  tree->Branch("Muon_pfRelIso04_dBeta", value_mu_pfrelisoDB, "Muon_pfRelIso04_dBeta[nMuon]/F");
  tree->Branch("Muon_pfRelIso04", value_mu_pfreliso, "Muon_pfRelIso04[nMuon]/F");
  tree->Branch("Muon_SIP3d", value_mu_SIP3d, "value_mu_SIP3d[nMuon]/F");
  tree->Branch("Muon_dxy", value_mu_dxy, "Muon_dxy[nMuon]/F");
  tree->Branch("Muon_dz", value_mu_dz, "Muon_dz[nMuon]/F");
  tree->Branch("Muon_isGM", value_mu_isGM, "Muon_isGM[nMuon]/F");
  // Electrons
  tree->Branch("nElectron", &value_el_n, "nElectron/i");
  tree->Branch("Electron_pt", value_el_pt, "Electron_pt[nElectron]/F");
  tree->Branch("Electron_eta", value_el_eta, "Electron_eta[nElectron]/F");
  tree->Branch("Electron_phi", value_el_phi, "Electron_phi[nElectron]/F");
  tree->Branch("Electron_rawE", value_el_rawE, "Electron_rawE[nElectron]/F");
  tree->Branch("Electron_charge", value_el_charge, "Electron_charge[nElectron]/I");
  tree->Branch("Electron_misshits", value_el_misshits, "Electron_misshits[nElectron]/F");
  tree->Branch("Electron_pfRelIso", value_el_pfreliso, "Electron_pfRelIso[nElectron]/F");
  tree->Branch("Electron_SIP3d", value_el_SIP3d, "value_mu_SIP3d[nElectron]/F");
  tree->Branch("Electron_dxy", value_el_dxy, "Electron_dxy[nElectron]/F");
  tree->Branch("Electron_dz", value_el_dz, "Electron_dz[nElectron]/F");

}

HiggsDemoAnalyzer::~HiggsDemoAnalyzer() {
  //do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------//
void HiggsDemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

// **********************************************
// here each relevant event will get analyzed 
// **********************************************

  value_mu_n = 0;
  value_el_n = 0;

  nRun  = iEvent.run();
  nEvt  = (iEvent.id()).event(); // iEvent: no class named event()
  nLumi = iEvent.luminosityBlock();
  
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif

  // Event is to be analyzed

  edm::LogInfo("Demo");
  /*<< "Starting to analyze \n"
  << "Event number: " << (iEvent.id()).event()
  << ", Run number: " << iEvent.run()
  << ", Lumisection: " << iEvent.luminosityBlock();*/

  //------------------Load (relevant) Event information------------------------//
  // INFO: Getting Data From an Event
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookChapter4#GetData
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMGetDataFromEvent#get_ByLabel
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAodDataTable
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideRecoDataTable

  // INFO: globalMuons
  // NB: note that when using keyword "globalMuons" getByLabel-function returns 
  // reco::TrackCollection
  //edm::Handle<reco::TrackCollection> tracks;
  //iEvent.getByLabel("generalTracks", tracks);

  edm::Handle<reco::TrackCollection> gmuons;
  iEvent.getByLabel("globalMuons", gmuons);

  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel("muons", muons);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);

  edm::Handle<reco::VertexCollection> primvtxHandle;
  iEvent.getByLabel("offlinePrimaryVertices", primvtxHandle);
  math::XYZPoint pv(primvtxHandle->begin()->position());

  edm::Handle<reco::GsfElectronCollection> electrons;
  iEvent.getByLabel("gsfElectrons",electrons);

  reco::BeamSpot beamSpot;
  if ( beamSpotHandle.isValid() )
    {
      beamSpot = *beamSpotHandle;

    } else
    {
      edm::LogInfo("Demo")
	<< "No beam spot available from EventSetup \n";
    }

  reco::VertexCollection primvtx;
  if ( primvtxHandle.isValid() )
    {
      primvtx = *primvtxHandle;

    } else
    {
      edm::LogInfo("Demo")
	<< "No primary vertex available from EventSetup \n";
    }


  // Declare vector that contain a pair of variables u wanna save
  // In this case: size and pT
  std::vector< std::pair<int, double> > vIdPt;
  std::vector< std::pair<int, double> > vIdPtmu;
  std::vector< std::pair<int, double> > vIdPte;

  // Initialize variables
  eZ12 = -9999.; eZ34 = -9999.; eZ13 = -9999.; eZ24 = -9999.; eZ14 = -9999.; eZ23 = -9999.; // select largest, init -
  pxZ12 = -9999.; pxZ34 = -9999.; pxZ13 = -9999.; pxZ24 = -9999.; pxZ14 = -9999.; pxZ23 = -9999.;
  pyZ12 = -9999.; pyZ34 = -9999.; pyZ13 = -9999.; pyZ24 = -9999.; pyZ14 = -9999.; pyZ23 = -9999.;
  pzZ12 = -9999.; pzZ34 = -9999.; pzZ13 = -9999.; pzZ24 = -9999.; pzZ14 = -9999.; pzZ23 = -9999.; 
  pZ12 = -9999.; pZ34 = -9999.; pZ13 = -9999.; pZ24 = -9999.; pZ14 = -9999.; pZ23 = -9999.; 
  pTZ12 = -9999.; pTZ34 = -9999.; pTZ13 = -9999.; pTZ24 = -9999.; pTZ14 = -9999.; pTZ23 = -9999.; 

  mZ12 = -9999.; mZ34 = -9999.; mZ13 = -9999.; mZ24 = -9999.; mZ14 = -9999.; mZ23 = -9999.;

  dZ12 = 9999.; dZ34 = 9999.; dZ13 = 9999.; dZ24 = 9999.; dZ14 = 9999.; dZ23 = 9999.; // select smallest, init +
  dZc1 = 9999.; dZc2 = 9999.; dZc3 = 9999.;

  eZa = -9999.; pxZa = -9999.; pyZa = -9999.; pzZa = -9999.; pTZa = -9999.; mZa = -9999.;
  eZb = -9999.; pxZb = -9999.; pyZb = -9999.; pzZb = -9999.; pTZb = -9999.; mZb = -9999.;
  mass4mu = -9999.; pt_4mu = -9999.; eta_4mu = -9999.; phi_4mu = -9999.;
  px4mu = -9999.; py4mu = -9999.; pz4mu = -9999.; E4mu = -9999.;
  
  pt_mu1 = -9999.; pt_mu2 = -9999.; pt_mu3 = -9999.; pt_mu4 = -9999.;
  eta_mu1 = -9999.; eta_mu2 = -9999.; eta_mu3 = -9999.; eta_mu4 = -9999.;
  phi_mu1 = -9999.; phi_mu2 = -9999.; phi_mu3 = -9999.; phi_mu4 = -9999.;
  cas_mu1 = -999; cas_mu2 = -999; cas_mu3 = -999; cas_mu4 = -999;

  px_mu1 = -9999.; px_mu2 = -9999.; px_mu3 = -9999.; px_mu4 = -9999.;
  py_mu1 = -9999.; py_mu2 = -9999.; py_mu3 = -9999.; py_mu4 = -9999.;
  pz_mu1 = -9999.; pz_mu2 = -9999.; pz_mu3 = -9999.; pz_mu4 = -9999.;

  E_mu1 = -9999.; E_mu2 = -9999.; E_mu3 = -9999.; E_mu4 = -9999.;

  mass4e = -9999.; pt_4e = -9999.; eta_4e = -9999.; phi_4e = -9999.;
  px4e = -9999.; py4e = -9999.; pz4e = -9999.; E4e = -9999.;
  
  pt_e1 = -9999.; pt_e2 = -9999.; pt_e3 = -9999.; pt_e4 = -9999.;
  eta_e1 = -9999.; eta_e2 = -9999.; eta_e3 = -9999.; eta_e4 = -9999.;
  phi_e1 = -9999.; phi_e2 = -9999.; phi_e3 = -9999.; phi_e4 = -9999.;
  cas_e1 = -999; cas_e2 = -999; cas_e3 = -999; cas_e4 = -999;

  px_e1 = -9999.; px_e2 = -9999.; px_e3 = -9999.; px_e4 = -9999.;
  py_e1 = -9999.; py_e2 = -9999.; py_e3 = -9999.; py_e4 = -9999.;
  pz_e1 = -9999.; pz_e2 = -9999.; pz_e3 = -9999.; pz_e4 = -9999.;

  E_e1 = -9999.; E_e2 = -9999.; E_e3 = -9999.; E_e4 = -9999.;
  
  s = -9999.;
  s1 = -9999.; s2 = -9999.; s3 = -9999.; s4 = -9999.;
  dx = 9999.; dy = 9999.; dz = 9999.; rap = -9999.; pz = -9999.;

  mass2mu2e = -9999.; pt_2mu2e = -9999.; eta_2mu2e = -9999.;  phi_2mu2e = -9999.;
  px2mu2e = -9999.; py2mu2e = -9999.; pz2mu2e = -9999.; E2mu2e = -9999.;
  
  pt_2mu1 = -9999.; pt_2mu2 = -9999.; pt_2e1 = -9999.; pt_2e2 = -9999.;
  eta_2mu1 = -9999.; eta_2mu2 = -9999.; eta_2e1 = -9999.; eta_2e2 = -9999.;
  phi_2mu1 = -9999.; phi_2mu2 = -9999.; phi_2e1 = -9999.; phi_2e2 = -9999.;
  
  cas_2mu1 = -999; cas_2mu2 = -999; cas_2e1 = -999; cas_2e2 = -999;

  px_2mu1 = -9999.; px_2mu2 = -9999.; px_2e1 = -9999.; px_2e2 = -9999.;
  py_2mu1 = -9999.; py_2mu2 = -9999.; py_2e1 = -9999.; py_2e2 = -9999.;
  pz_2mu1 = -9999.; pz_2mu2 = -9999.; pz_2e1 = -9999.; pz_2e2 = -9999.;

  E_2mu1 = -9999.; E_2mu2 = -9999.; E_2e1 = -9999.; E_2e2 = -9999.;
  
  goodhit = -9999.;
  relPFIso_mu = -9999.;
  relPFIso_e = -9999.;

  IP3d_mu = -9999;
  ErrIP3d_mu = -9999;
  SIP3d_mu = -9999.;
  IP3d_e = -9999;
  ErrIP3d_e = -9999;
  SIP3d_e = -9999.;

  p4Za.SetPxPyPzE(0., 0., 0., 0.);
  p4Zb.SetPxPyPzE(0., 0., 0., 0.);
  p4H.SetPxPyPzE(0., 0., 0., 0.);

  // constant value squared muon mass, squared electron mass and Z mass
  sqm1 = (0.105658) * (0.105658);
  sqme = (0.0005109989) * (0.0005109989);
  mZ = 91.1876;

  /*h_globalmu_size->Fill(gmuons->size());
  h_recomu_size->Fill(muons->size());
  h_e_size->Fill(electrons->size());*/

  ///////////////////////////////////////////////////////////////////////////////
  /////////////////////// Global Muon Collection Start //////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  //******************************************************************************
  // This is some code similar to the one of the 'dimuon analysis example',      *
  // illustrating the use of global muons, but not directly used for the Higgs   *
  // analysis                                                                    *
  //******************************************************************************

  // Loop over global muons size and select good muons

  for (unsigned t = 0; t < gmuons->size(); t++)
    {
      const reco::Track &iMuon = (*gmuons)[t];
      const reco::HitPattern& p = iMuon.hitPattern();

      /*h_p_gmu->Fill(iMuon.p());
      h_pt_gmu_b4->Fill(iMuon.pt());
      h_eta_gmu_b4->Fill(iMuon.eta());
      h_phi_gmu->Fill(iMuon.phi());

      h_chi2_gmu->Fill(iMuon.chi2());
      h_ndof_gmu->Fill(iMuon.ndof());
      h_normchi2_gmu->Fill(iMuon.normalizedChi2());*/

      // Counter
      int GM_ValidHits = 0;
      int GM_PixelHits = 0;
   
      // loop over the hits of the track
      for (int i = 0; i < p.numberOfHits(); i++) 
	{
	  uint32_t hit = p.getHitPattern(i);

	  // if the hit is valid and in pixel
	  if (p.validHitFilter(hit) && p.pixelHitFilter(hit)) {GM_PixelHits++;}
	  if (p.validHitFilter(hit)) {GM_ValidHits++;}
	}

      /*h_validhits_gmu->Fill(GM_ValidHits);
      h_pixelhits_gmu->Fill(GM_PixelHits);*/

      if (GM_ValidHits >= 12 && GM_PixelHits >= 2 && iMuon.normalizedChi2() < 4.0)
	{
	  // Save a vector that contains 2 variables: gmuon size (.first) and the pT (.second)
	  vIdPt.push_back( std::make_pair(t, iMuon.pt()) );
	}
    }

  // Sort the pT (.second) to decending order (from highest pT to lowest pT)
  std::sort(vIdPt.begin(), vIdPt.end(), [](const std::pair<int, double> &idPt1, const std::pair<int, double> &idPt2) {return (idPt1.second > idPt2.second);});

  ///////////////////////////////////////////////////////////////////////////////
  //////////////////////// Global Muon Collection End ///////////////////////////
  ///////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////////
  ///////////////////////// Reco Muon Collection Start ///////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  
  //******************************************************************************
  // This muon collection is being used here for the Higgs->4 lepton analysis    *
  //******************************************************************************

  // Loop over muons size and select good muons

  for (unsigned u = 0; u < muons->size(); u++)
    {
      const reco::Muon &itMuon = (*muons)[u];

      // math::XYZPoint point(beamSpot.position());
      math::XYZPoint point(primvtx[0].position());
      
      // select global particle flow muons
      // some muons might not have valid track references
      if (itMuon.isPFMuon() && itMuon.isPFIsolationValid() && (itMuon.globalTrack()).isNonnull())
	{


    

	  /*h_p_reco->Fill(itMuon.p());
	  h_pt_reco_b4->Fill(itMuon.pt());
	  h_eta_reco_b4->Fill(itMuon.eta());
	  h_phi_reco->Fill(itMuon.phi());

	  h_chi2_reco->Fill((itMuon.globalTrack())->chi2());
	  h_ndof_reco->Fill((itMuon.globalTrack())->ndof());
	  h_normchi2_reco->Fill((itMuon.globalTrack())->normalizedChi2());*/

	  //======= Use Particle Flow (PF) Muon & PF relative isolation ========//

	  // PF Relative isolation for muons
	  relPFIso_mu = ((itMuon.pfIsolationR04()).sumChargedHadronPt +
			 (itMuon.pfIsolationR04()).sumNeutralHadronEt + 
			 (itMuon.pfIsolationR04()).sumPhotonEt) / itMuon.pt(); 

	  // h_relPFIso_mu->Fill(relPFIso_mu);
	  
	  // Checking hit pattern info
	  const reco::HitPattern& RM_p = (itMuon.globalTrack())->hitPattern();

	  goodhit = RM_p.numberOfValidMuonHits();
	  // h_goodhit->Fill(goodhit);

	  // h_dxy_mu->Fill((itMuon.globalTrack())->dxy(point));

	  IP3d_mu = sqrt((itMuon.globalTrack()->dxy(point) * itMuon.globalTrack()->dxy(point)) + (itMuon.globalTrack()->dz(point) * itMuon.globalTrack()->dz(point)));
	  
	  ErrIP3d_mu = sqrt((itMuon.globalTrack()->d0Error() * itMuon.globalTrack()->d0Error()) + (itMuon.globalTrack()->dzError() * itMuon.globalTrack()->dzError()));
	  
	  SIP3d_mu = IP3d_mu / ErrIP3d_mu;

	  // h_SIP3d_mu_b4->Fill(SIP3d_mu);

	  int RM_ValidHits = 0;
	  int RM_PixelHits = 0;

	  for (int i = 0; i < RM_p.numberOfHits(); i++) {
	    uint32_t hit = RM_p.getHitPattern(i);

	    // If the hit is valid and in pixel
	    if (RM_p.validHitFilter(hit) && RM_p.pixelHitFilter(hit))
	      {RM_PixelHits++;}
	    
	      if (RM_p.validHitFilter(hit)) {RM_ValidHits++;}
	    }

	  // h_validhits_mu->Fill(RM_ValidHits);
	  // h_pixelhits_mu->Fill(RM_PixelHits);

	  if (std::abs(SIP3d_mu) < 4. && std::abs((itMuon.globalTrack())->dxy(point)) < 0.5 && std::abs((itMuon.globalTrack())->dz(point)) < 1. && relPFIso_mu < 0.4)
	    {
	      if (itMuon.pt() > 5. && std::abs(itMuon.eta()) < 2.4)
		{
		  vIdPtmu.push_back( std::make_pair(u, itMuon.pt()) );

      
		}
	    }
	} // end of if (itMuon.isPFMuon().........
    } // end muons loop

  // Sort the pT (.second) to decending order (from highest pT to lowest pT)
  std::sort(vIdPtmu.begin(), vIdPtmu.end(), [] (const std::pair<int, double> &idPtmu1, const std::pair<int, double> &idPtmu2) {return (idPtmu1.second > idPtmu2.second);});

  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////// Muon Collection end /////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  
  ///////////////////////////////////////////////////////////////////////////////
  ////////////////////// Electron Collection Start //////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  for (unsigned te = 0; te < electrons->size(); te++)
    {
      const reco::GsfElectron &iElectron = (*electrons)[te];

      // math::XYZPoint point(beamSpot.position());
       math::XYZPoint point(primvtx[0].position());

       if (iElectron.passingPflowPreselection())
	 {
	   misshits = ((iElectron.gsfTrack())->trackerExpectedHitsInner()).numberOfHits();

	   IP3d_e = sqrt ( (iElectron.gsfTrack()->dxy(point) * iElectron.gsfTrack()->dxy(point)) + (iElectron.gsfTrack()->dz(point) * iElectron.gsfTrack()->dz(point)) );
	   ErrIP3d_e = sqrt ( (iElectron.gsfTrack()->d0Error() * iElectron.gsfTrack()->d0Error()) + (iElectron.gsfTrack()->dzError() * iElectron.gsfTrack()->dzError()) );
	   SIP3d_e = IP3d_e / ErrIP3d_e;

	   /*h_SIP3d_e_b4->Fill(SIP3d_e);

	   h_p_e->Fill(iElectron.p());
	   h_et_e->Fill(iElectron.et());
	   h_pt_e_b4->Fill(iElectron.pt());
	   h_eta_e_b4->Fill(iElectron.eta());
	   h_phi_e->Fill(iElectron.phi());
	   h_sc_eta->Fill((iElectron.superCluster())->eta());
	   h_sc_rawE->Fill(std::abs((iElectron.superCluster())->rawEnergy()));
	   h_misshite->Fill(misshits);*/
      
	   // Relative isolation for electron
	   relPFIso_e = ((iElectron.pfIsolationVariables()).chargedHadronIso +
			 (iElectron.pfIsolationVariables()).neutralHadronIso +
			 (iElectron.pfIsolationVariables()).photonIso) /iElectron.pt();
	
	   /*h_relPFIso_e->Fill(relPFIso_e);

	   h_relPFIso_pt_e->Fill(relPFIso_e, iElectron.pt());

	   h_dxy_e->Fill((iElectron.gsfTrack())->dxy(point));*/

	   // Electron selection
	   if (iElectron.pt() > 7.)
	     {
	       if ((std::abs((iElectron.superCluster())->eta())) < 2.5)
		 {
		   if (misshits <= 1 && std::abs(SIP3d_e) < 4.)
		     {
		       if (std::abs(iElectron.gsfTrack()->dxy(point)) < 0.5 && std::abs(iElectron.gsfTrack()->dz(point)) < 1.)
			 {
			   if (iElectron.isEB())
			     {
			       if (relPFIso_e < 0.4)
				 {
				   vIdPte.push_back( std::make_pair(te, iElectron.pt()) );
				 }
			     }
			   else if (iElectron.isEE())
			     {
			       if (relPFIso_e < 0.4)
				 {
				   vIdPte.push_back( std::make_pair(te, iElectron.pt()) );
				 }
			     }
			 }
		     }
		 }
	     }
	 }
    } // for (unsigned te = 0; te < electrons->size(); te++)
  
  // Sort the pT (.second) to decending order (from highest pT to lowest pT)
  std::sort(vIdPte.begin(), vIdPte.end(), [](const std::pair<int, double> &idPte1, const std::pair<int, double> &idPte2) {return (idPte1.second > idPte2.second);});

  ////////////////////////////////////////////////////////////////////////////////
  ///////////////////////// Electron Collection end //////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  // these index (vIdPt...) is declare as good candidates in integer
  nGoodGlobalMuon = vIdPt.size(); 
  nGoodRecoMuon = vIdPtmu.size(); 
  nGoodElectron = vIdPte.size(); 

  // h_nggmu->Fill(nGoodGlobalMuon);
  // h_ngmu->Fill(nGoodRecoMuon);
  // h_nge->Fill(nGoodElectron);

  ////////////////////////////////////////////////////////////////////////////////
  /////////////////////// All calculation start here /////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  
  //====================== Dimuon using Global Muon start ======================//

  // For the case of nGoodGlobalMuon > 2, the selected two muons are always have the highest pT as we have sorted it before
  if (nGoodGlobalMuon >= 2)
    {
      const reco::Track &gmuon1 = (*gmuons)[vIdPt.at(0).first];
      const reco::Track &gmuon2 = (*gmuons)[vIdPt.at(1).first];
      
      // The sum of 2 charge global muon should be 0 (neutral)
      if (gmuon1.charge() + gmuon2.charge() == 0)
	{
	  for (unsigned i = 0; i < vIdPt.size(); i++)
	    {
	      // These pT and eta are filled after all the cuts
	      // h_pt_gmu_after->Fill(vIdPt.at(i).second); // access directly .second as the second pair is already pT
	      // h_eta_gmu_after->Fill(((*gmuons)[vIdPt.at(i).first]).eta());
	    }

	  s1 = sqrt(((gmuon1.p()) * (gmuon1.p()) + sqm1) * ((gmuon2.p()) * (gmuon2.p()) + sqm1));
	  s2 = gmuon1.px() * gmuon2.px() + gmuon1.py() * gmuon2.py() + gmuon1.pz() * gmuon2.pz();
	  s = sqrt(2.0 * (sqm1 + (s1 - s2)));

	  /*h_m1_gmu->Fill(s);
	  h_m2_gmu->Fill(s);
	  h_m3_gmu->Fill(s);*/

	}
    }

  //======================== Dimuon using Global Muon end ======================//

  //===================== ZTo2Muon using Reco Muon start =======================//

  if (nGoodRecoMuon >= 2)
    {
      const reco::Muon &muon1 = (*muons)[vIdPtmu.at(0).first];
      const reco::Muon &muon2 = (*muons)[vIdPtmu.at(1).first];

      if (muon1.charge() + muon2.charge() == 0)
	{
	  for (unsigned i = 0; i < vIdPtmu.size(); i++)
	    {
	      // These pT and eta are filled after all the cuts
	      // access directly .second as the second pair is already pT
	      // h_pt_after_Zto2mu->Fill(vIdPtmu.at(i).second); 
	      // h_eta_after_Zto2mu->Fill(((*muons)[vIdPtmu.at(i).first]).eta());
	    }

	  s1 = sqrt(((muon1.p()) * (muon1.p()) + sqm1) * ((muon2.p()) * (muon2.p()) + sqm1));
	  s2 = muon1.px() * muon2.px() + muon1.py() * muon2.py() + muon1.pz() * muon2.pz();
	  s = sqrt(2.0 * (sqm1 + (s1 - s2)));

	  // h_mZ_2mu->Fill(s);
	}
    }

  //===================== ZTo2Muon using Reco Muon end =========================//


  //============================ ZTo2Electron start ============================//

  if (nGoodElectron >= 2)
    {
      const reco::GsfElectron &elec1 = (*electrons)[vIdPte.at(0).first];
      const reco::GsfElectron &elec2 = (*electrons)[vIdPte.at(1).first];

      if (elec1.charge() + elec2.charge() == 0)
	{
	  for (unsigned i = 0; i < vIdPte.size(); i++)
	    {
	      // These pT and eta are filled after all the cuts
	      // access directly .second as the second pair is already pT
	      // h_pt_e_after_Zto2e->Fill(vIdPte.at(i).second); 
	      // h_eta_e_after_Zto2e->Fill((((*electrons)[vIdPte.at(i).first]).superCluster())->eta());
	    }

	  s1 = sqrt(((elec1.p()) * (elec1.p()) + sqme) * ((elec2.p()) * (elec2.p()) + sqme));
	  s2 = elec1.px() * elec2.px() + elec1.py() * elec2.py() + elec1.pz() * elec2.pz();
	  s = sqrt(2.0 * (sqme + (s1 - s2)));

	  // h_mZ_2e->Fill(s);
	}
    }

  //========================== ZTo2Electron end ================================//


  //======================== ZZ/ZZ*To4Muon start ===============================//

  // Now, for these goodmuons, pair up and calculate mass
  if (nGoodRecoMuon >= 4)
    {
      const reco::Muon &muon1 = (*muons)[vIdPtmu.at(0).first];
      const reco::Muon &muon2 = (*muons)[vIdPtmu.at(1).first];
      const reco::Muon &muon3 = (*muons)[vIdPtmu.at(2).first];
      const reco::Muon &muon4 = (*muons)[vIdPtmu.at(3).first];

      if (muon1.charge() + muon2.charge() + muon3.charge() + muon4.charge() == 0)
	{
	  // First combination: Combine muon 1234
	  if (muon1.charge() + muon2.charge() == 0) // each lepton pair cas = 0
	    {
	      eZ12 = (sqrt(muon1.p() * muon1.p() + sqm1)) +
		     (sqrt(muon2.p() * muon2.p() + sqm1));
	
	      pxZ12 = muon1.px() + muon2.px();
	      pyZ12 = muon1.py() + muon2.py();
	      pzZ12 = muon1.pz() + muon2.pz();

	      if (muon3.charge() + muon4.charge() == 0)
		{
		  eZ34 = (sqrt(muon3.p() * muon3.p() + sqm1)) +
		         (sqrt(muon4.p() * muon4.p() + sqm1));

		  pxZ34 = muon3.px() + muon4.px();
		  pyZ34 = muon3.py() + muon4.py();
		  pzZ34 = muon3.pz() + muon4.pz();

		  // Calculate p4
		  pZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12) + (pzZ12 * pzZ12));
		  pZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34) + (pzZ34 * pzZ34));

		  pTZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12));
		  pTZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34));

		  mZ12 = sqrt((eZ12 * eZ12) - (pZ12 * pZ12));
		  mZ34 = sqrt((eZ34 * eZ34) - (pZ34 * pZ34));

		  // if (mZ12 > 0.) h_mZ12_4mu->Fill(mZ12);
		  // if (mZ34 > 0.) h_mZ34_4mu->Fill(mZ34);
		}
	    }

	  dZ12 = std::abs( mZ12 - mZ );
	  dZ34 = std::abs( mZ34 - mZ );

	  // take the smallest difference between mass
	  // to use for 4muon combination
	  dZc1 = (dZ12 < dZ34) ? dZ12 : dZ34; 
 
	  // Second combination: Combine muon 1324
	  if (muon1.charge() + muon3.charge() == 0)
	    {
	      eZ13 = (sqrt(muon1.p() * muon1.p() + sqm1)) +
		     (sqrt(muon3.p() * muon3.p() + sqm1));

	      pxZ13 = muon1.px() + muon3.px();
	      pyZ13 = muon1.py() + muon3.py();
	      pzZ13 = muon1.pz() + muon3.pz();

	      if (muon2.charge() + muon4.charge() == 0)
		{
		  eZ24 = (sqrt(muon2.p() * muon2.p() + sqm1)) +
		         (sqrt(muon4.p() * muon4.p() + sqm1));

		  pxZ24 = muon2.px() + muon4.px();
		  pyZ24 = muon2.py() + muon4.py();
		  pzZ24 = muon2.pz() + muon4.pz();

		  // Calculate p4
		  pZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13) + (pzZ13 * pzZ13));
		  pZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24) + (pzZ24 * pzZ24));

		  pTZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13));
		  pTZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24));

		  mZ13 = sqrt((eZ13 * eZ13) - (pZ13 * pZ13));
		  mZ24 = sqrt((eZ24 * eZ24) - (pZ24 * pZ24));

		  // if (mZ13 > 0.) h_mZ13_4mu->Fill(mZ13);
		  // if (mZ24 > 0.) h_mZ24_4mu->Fill(mZ24);
		}
	    }

	  dZ13 = std::abs( mZ13 - mZ );
	  dZ24 = std::abs( mZ24 - mZ );

	  dZc2 = (dZ13 < dZ24) ? dZ13 : dZ24; 

	  // Third combination: Combine muon 1423
	  if (muon1.charge() + muon4.charge() == 0)
	    {
	      eZ14 = (sqrt(muon1.p() * muon1.p() + sqm1)) +
		     (sqrt(muon4.p() * muon4.p() + sqm1));

	      pxZ14 = muon1.px() + muon4.px();
	      pyZ14 = muon1.py() + muon4.py();
	      pzZ14 = muon1.pz() + muon4.pz();

	      if (muon2.charge() + muon3.charge() == 0)
		{
		  eZ23 = sqrt((muon2.p() * muon2.p() + sqm1)) +
		         (sqrt(muon3.p() * muon3.p() + sqm1));
       
		  pxZ23 = muon2.px() + muon3.px();
		  pyZ23 = muon2.py() + muon3.py();
		  pzZ23 = muon2.pz() + muon3.pz();

		  // Calculate p4
		  pZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14) + (pzZ14 * pzZ14));
		  pZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23) + (pzZ23 * pzZ23));

		  pTZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14));
		  pTZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23));

		  mZ14 = sqrt((eZ14 * eZ14) - (pZ14 * pZ14));
		  mZ23 = sqrt((eZ23 * eZ23) - (pZ23 * pZ23));

		  // if (mZ14 > 0.) h_mZ14_4mu->Fill(mZ14);
		  // if (mZ23 > 0.) h_mZ23_4mu->Fill(mZ23);
		}
	    }

	  dZ14 = std::abs( mZ14 - mZ );
	  dZ23 = std::abs( mZ23 - mZ );

	  dZc3 = (dZ14 < dZ23) ? dZ14 : dZ23;

	  bool ptZadaug = false;

	  if (dZc1 < dZc2 && dZc1 < dZc3)
	    {
	      if (dZ12 < dZ34)
		{
		  eZa  = eZ12;     
		  pxZa = pxZ12;
		  pyZa = pyZ12;
		  pzZa = pzZ12;
		  pTZa = pTZ12;
		  mZa  = mZ12;

		  if (muon1.pt() > 20. and muon2.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ34;
		  pxZb = pxZ34;
		  pyZb = pyZ34;
		  pzZb = pzZ34;
		  pTZb = pTZ34;
		  mZb  = mZ34;
		}
	      else
		{
		  eZa  = eZ34;
		  pxZa = pxZ34;
		  pyZa = pyZ34;
		  pzZa = pzZ34;
		  pTZa = pTZ34;
		  mZa  = mZ34;

		  if (muon3.pt() > 20. and muon4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ12;
		  pxZb = pxZ12;
		  pyZb = pyZ12;
		  pzZb = pzZ12;
		  pTZb = pTZ12;
		  mZb  = mZ12;
		}
	    }

	  else if (dZc2 < dZc1 && dZc2 < dZc3)
	    {
	      if (dZ13 < dZ24)
		{
		  eZa  = eZ13;
		  pxZa = pxZ13;
		  pyZa = pyZ13;
		  pzZa = pzZ13;
		  pTZa = pTZ13;
		  mZa  = mZ13;

		  if (muon1.pt() > 20. and muon3.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ24;
		  pxZb = pxZ24;
		  pyZb = pyZ24;
		  pzZb = pzZ24;
		  pTZb = pTZ24;
		  mZb  = mZ24;
		}
	      else
		{
		  eZa  = eZ24;
		  pxZa = pxZ24;
		  pyZa = pyZ24;
		  pzZa = pzZ24;
		  pTZa = pTZ24;
		  mZa  = mZ24;

		  if (muon2.pt() > 20. and muon4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ13;
		  pxZb = pxZ13;
		  pyZb = pyZ13;
		  pzZb = pzZ13;
		  pTZb = pTZ13;
		  mZb  = mZ13;
		}
	    }

	  else if (dZc3 < dZc1 && dZc3 < dZc2)
	    {
	      if (dZ14 < dZ23)
		{
		  eZa  = eZ14;
		  pxZa = pxZ14;
		  pyZa = pyZ14;
		  pzZa = pzZ14;
		  pTZa = pTZ14;
		  mZa  = mZ14;

		  if (muon1.pt() > 20. and muon4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ23;
		  pxZb = pxZ23;
		  pyZb = pyZ23;
		  pzZb = pzZ23;
		  pTZb = pTZ23;
		  mZb  = mZ23;
		}
	      else
		{
		  eZa  = eZ23;
		  pxZa = pxZ23;
		  pyZa = pyZ23;
		  pzZa = pzZ23;
		  pTZa = pTZ23;
		  mZa  = mZ23;

		  if (muon2.pt() > 20. and muon3.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ14;
		  pxZb = pxZ14;
		  pyZb = pyZ14;
		  pzZb = pzZ14;
		  pTZb = pTZ14;
		  mZb  = mZ14;
		}
	    }


    // put it in this section - do a loop or something to get all four objects into the tree - similar to what is happening below

    
        //value_mu_isGM[value_mu_n] = it->isGlobalMuon();
        //auto iso04 = it->pfIsolationR04();
        //value_mu_pfrelisoDB[value_mu_n] = (iso04.sumChargedHadronPt + iso04.sumNeutralHadronEt + iso04.sumPhotonEt - 0.5*iso04.sumPUPt)/it->pt(); // OD version is not DB corrected
        //value_mu_pfreliso[value_mu_n] = (iso04.sumChargedHadronPt + iso04.sumNeutralHadronEt + iso04.sumPhotonEt)/it->pt(); // OD version
        //value_mu_dxy[value_mu_n] = trk->dxy(pv);
        //value_mu_dz[value_mu_n] = trk->dz(pv);
        //float d3 = sqrt( pow(trk->dxy(pv),2) + pow(trk->dz(pv),2) );
        //float d3Err = sqrt( pow(trk->d0Error(),2) + pow(trk->dzError(),2) );
        //value_mu_SIP3d[value_mu_n] = d3/d3Err;


    value_mu_n = 4;
        
    value_mu_pt[0] = muon1.pt();
    value_mu_pt[1] = muon2.pt();
    value_mu_pt[2] = muon3.pt();
    value_mu_pt[3] = muon4.pt();

    value_mu_eta[0] = muon1.eta();
		value_mu_eta[1] = muon2.eta();
		value_mu_eta[2] = muon3.eta();
		value_mu_eta[3] = muon4.eta();

    value_mu_phi[0] = muon1.phi();
		value_mu_phi[1] = muon2.phi();
		value_mu_phi[2] = muon3.phi();
		value_mu_phi[3] = muon4.phi();

		value_mu_charge[0] = muon1.charge();
		value_mu_charge[1] = muon2.charge();
		value_mu_charge[2] = muon3.charge();
		value_mu_charge[3] = muon4.charge();

    value_mu_mass[0] = mass4mu;
		value_mu_mass[1] = mass4mu;
		value_mu_mass[2] = mass4mu;
		value_mu_mass[3] = mass4mu;

    //I'm not sure about this one
    value_mu_isGM[0] = muon1.isGlobalMuon();
		value_mu_isGM[1] = muon2.isGlobalMuon();
		value_mu_isGM[2] = muon3.isGlobalMuon();
		value_mu_isGM[3] = muon4.isGlobalMuon();

    //Or this one
    value_mu_pfrelisoDB[0] = ((muon1.pfIsolationR04()).sumChargedHadronPt + (muon1.pfIsolationR04()).sumNeutralHadronEt + (muon1.pfIsolationR04()).sumPhotonEt - 0.5*(muon1.pfIsolationR04()).sumPUPt) / (muon1.pt());
		value_mu_pfrelisoDB[1] = ((muon2.pfIsolationR04()).sumChargedHadronPt + (muon2.pfIsolationR04()).sumNeutralHadronEt + (muon2.pfIsolationR04()).sumPhotonEt - 0.5*(muon2.pfIsolationR04()).sumPUPt) / (muon2.pt());
		value_mu_pfrelisoDB[2] = ((muon3.pfIsolationR04()).sumChargedHadronPt + (muon3.pfIsolationR04()).sumNeutralHadronEt + (muon3.pfIsolationR04()).sumPhotonEt - 0.5*(muon3.pfIsolationR04()).sumPUPt) / (muon3.pt());
		value_mu_pfrelisoDB[3] = ((muon4.pfIsolationR04()).sumChargedHadronPt + (muon4.pfIsolationR04()).sumNeutralHadronEt + (muon4.pfIsolationR04()).sumPhotonEt - 0.5*(muon4.pfIsolationR04()).sumPUPt) / (muon4.pt());

    //This one as well
    value_mu_pfreliso[0] = ((muon1.pfIsolationR04()).sumChargedHadronPt + (muon1.pfIsolationR04()).sumNeutralHadronEt + (muon1.pfIsolationR04()).sumPhotonEt) / (muon1.pt());
		value_mu_pfreliso[1] = ((muon2.pfIsolationR04()).sumChargedHadronPt + (muon2.pfIsolationR04()).sumNeutralHadronEt + (muon2.pfIsolationR04()).sumPhotonEt) / (muon2.pt());
		value_mu_pfreliso[2] = ((muon3.pfIsolationR04()).sumChargedHadronPt + (muon3.pfIsolationR04()).sumNeutralHadronEt + (muon3.pfIsolationR04()).sumPhotonEt) / (muon3.pt());
		value_mu_pfreliso[3] = ((muon4.pfIsolationR04()).sumChargedHadronPt + (muon4.pfIsolationR04()).sumNeutralHadronEt + (muon4.pfIsolationR04()).sumPhotonEt) / (muon4.pt());

    value_mu_dxy[0] = muon1.globalTrack()->dxy(pv);
		value_mu_dxy[1] = muon2.globalTrack()->dxy(pv);
		value_mu_dxy[2] = muon3.globalTrack()->dxy(pv);
		value_mu_dxy[3] = muon4.globalTrack()->dxy(pv);

    value_mu_dz[0] = muon1.globalTrack()->dz(pv);
		value_mu_dz[1] = muon2.globalTrack()->dz(pv);
		value_mu_dz[2] = muon3.globalTrack()->dz(pv);
		value_mu_dz[3] = muon4.globalTrack()->dz(pv);

    value_mu_SIP3d[0] = (sqrt(pow(muon1.globalTrack()->dxy(pv),2) + pow(muon1.globalTrack()->dz(pv),2)))/(sqrt(pow(muon1.globalTrack()->d0Error(),2) + pow(muon1.globalTrack()->dzError(),2)));
    value_mu_SIP3d[1] = (sqrt(pow(muon2.globalTrack()->dxy(pv),2) + pow(muon2.globalTrack()->dz(pv),2)))/(sqrt(pow(muon2.globalTrack()->d0Error(),2) + pow(muon2.globalTrack()->dzError(),2)));
    value_mu_SIP3d[2] = (sqrt(pow(muon3.globalTrack()->dxy(pv),2) + pow(muon3.globalTrack()->dz(pv),2)))/(sqrt(pow(muon3.globalTrack()->d0Error(),2) + pow(muon3.globalTrack()->dzError(),2)));
    value_mu_SIP3d[3] = (sqrt(pow(muon4.globalTrack()->dxy(pv),2) + pow(muon4.globalTrack()->dz(pv),2)))/(sqrt(pow(muon4.globalTrack()->d0Error(),2) + pow(muon4.globalTrack()->dzError(),2)));

    muonTree = true;



	   if (ptZadaug) {
	    if (mZa > 40. && mZa < 120.) {
	      if (mZb > 12. && mZb < 120.) {

		// h_mZa_4mu->Fill(mZa);
		// h_mZb_4mu->Fill(mZb);

		// 4 vector
		p4Za.SetPxPyPzE(pxZa, pyZa, pzZa, eZa);
		p4Zb.SetPxPyPzE(pxZb, pyZb, pzZb, eZb);

		p4H = p4Za + p4Zb;

		mass4mu = p4H.M();
		pt_4mu = p4H.Pt();
		eta_4mu = p4H.Eta();
		phi_4mu = p4H.Phi();

		px4mu = p4H.Px();
		py4mu = p4H.Py();
		pz4mu = p4H.Pz();
		E4mu = p4H.E();

		pt_mu1 = muon1.pt();
		pt_mu2 = muon2.pt();
		pt_mu3 = muon3.pt();
		pt_mu4 = muon4.pt();

		eta_mu1 = muon1.eta();
		eta_mu2 = muon2.eta();
		eta_mu3 = muon3.eta();
		eta_mu4 = muon4.eta();

		phi_mu1 = muon1.phi();
		phi_mu2 = muon2.phi();
		phi_mu3 = muon3.phi();
		phi_mu4 = muon4.phi();

		cas_mu1 = muon1.charge();
		cas_mu2 = muon2.charge();
		cas_mu3 = muon3.charge();
		cas_mu4 = muon4.charge();

		px_mu1 = muon1.px();
		px_mu2 = muon2.px();
		px_mu3 = muon3.px();
		px_mu4 = muon4.px();

		py_mu1 = muon1.py();
		py_mu2 = muon2.py();
		py_mu3 = muon3.py();
		py_mu4 = muon4.py();
		
		pz_mu1 = muon1.pz();
		pz_mu2 = muon2.pz();
		pz_mu3 = muon3.pz();
		pz_mu4 = muon4.pz();

		E_mu1 = sqrt((muon1.p() * muon1.p()) + sqm1);
		E_mu2 = sqrt((muon2.p() * muon2.p()) + sqm1);
		E_mu3 = sqrt((muon3.p() * muon3.p()) + sqm1);
		E_mu4 = sqrt((muon4.p() * muon4.p()) + sqm1);




		if (mass4mu > 70.)
		  {
		    /*h_m1_m4mu->Fill(mass4mu);
		    h_m2_m4mu->Fill(mass4mu);
		    h_m3_m4mu->Fill(mass4mu);
		    h_m4_m4mu->Fill(mass4mu);*/

        if (E_mu1 != -999) {
          if (myfile4mu.is_open()){
			      myfile4mu << E_mu1 << ", " << E_mu2 << ", " << E_mu3 << ", " << E_mu4 << ", " << px_mu1 << ", " << px_mu2 << ", " << px_mu3 << ", " << px_mu4 << ", " << py_mu1 << ", " << py_mu2 << ", " << py_mu3 << ", " << py_mu4 << ", " << pz_mu1 << ", " << pz_mu2 << ", " << pz_mu3 << ", " << pz_mu4 << "\n";
          }
        }

		    for (unsigned i = 0; i < vIdPtmu.size(); i++)
		      {
			relPFIso_mu = ((((*muons)[vIdPtmu.at(i).first]).pfIsolationR04()).sumChargedHadronPt +
				       (((*muons)[vIdPtmu.at(i).first]).pfIsolationR04()).sumNeutralHadronEt +
				       (((*muons)[vIdPtmu.at(i).first]).pfIsolationR04()).sumPhotonEt) / (((*muons)[vIdPtmu.at(i).first]).pt()); 

			/*h_relPFIso_mu_after->Fill(relPFIso_mu);

			h_pt_after->Fill(vIdPtmu.at(i).second);
			h_eta_after->Fill(((*muons)[vIdPtmu.at(i).first]).eta());*/
		      }
		    // t1->Fill();
		  }
	      } // end of mZb
	    } // end of mZa
	   } // end of ptZadaug
	} // end of total charge
    } // end of nGoodRecoMuon

  //============================ ZZ/ZZ*To4Muon end =============================//

  //============================ ZZ/ZZ*To4e start ==============================//

  // Now, for these goodelectrons, pair up and calculate mass
  if (nGoodElectron >= 4)
    {
	  myfile4e << "thing three" << endl;
      const reco::GsfElectron &elec1 = (*electrons)[vIdPte.at(0).first];
      const reco::GsfElectron &elec2 = (*electrons)[vIdPte.at(1).first];
      const reco::GsfElectron &elec3 = (*electrons)[vIdPte.at(2).first];
      const reco::GsfElectron &elec4 = (*electrons)[vIdPte.at(3).first];

      if (elec1.charge() + elec2.charge() + elec3.charge() + elec4.charge() == 0)
	{
	  // First combination: Combine elec 1234
	  if (elec1.charge() + elec2.charge() == 0)
	    {
	      eZ12 = (sqrt(elec1.p() * elec1.p() + sqme)) + (sqrt(elec2.p() * elec2.p() + sqme));
	
	      pxZ12 = elec1.px() + elec2.px();
	      pyZ12 = elec1.py() + elec2.py();
	      pzZ12 = elec1.pz() + elec2.pz();

	      if (elec3.charge() + elec4.charge() == 0)
		{
		  eZ34 = (sqrt(elec3.p() * elec3.p() + sqme)) + (sqrt(elec4.p() * elec4.p() + sqme));

		  pxZ34 = elec3.px() + elec4.px();
		  pyZ34 = elec3.py() + elec4.py();
		  pzZ34 = elec3.pz() + elec4.pz();

		  // Calculate the momentum and invariant mass of elec 1 and 2, elec 3 and 4
		  pZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12) + (pzZ12 * pzZ12));
		  pZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34) + (pzZ34 * pzZ34));

		  pTZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12));
		  pTZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34));

		  mZ12 = sqrt((eZ12 * eZ12) - (pZ12 * pZ12));
		  mZ34 = sqrt((eZ34 * eZ34) - (pZ34 * pZ34));

		  // if (mZ12 > 0.) h_mZ12_4e->Fill(mZ12);
		  // if (mZ34 > 0.) h_mZ34_4e->Fill(mZ34);
		}
	    }

	  dZ12 = std::abs( mZ12 - mZ );
	  dZ34 = std::abs( mZ34 - mZ );

	  // take the smallest diff between mass to use for 4electron combination
	  dZc1 = (dZ12 < dZ34) ? dZ12 : dZ34; 

	  // Second combination: Combine elec 1324
	  if (elec1.charge() + elec3.charge() == 0)
	    {
	      eZ13 = (sqrt(elec1.p() * elec1.p() + sqme)) + (sqrt(elec3.p() * elec3.p() + sqme));

	      pxZ13 = elec1.px() + elec3.px();
	      pyZ13 = elec1.py() + elec3.py();
	      pzZ13 = elec1.pz() + elec3.pz();

	      if (elec2.charge() + elec4.charge() == 0)
		{

		  eZ24 = (sqrt(elec2.p() * elec2.p() + sqme)) + (sqrt(elec4.p() * elec4.p() + sqme));
		  
		  pxZ24 = elec2.px() + elec4.px();
		  pyZ24 = elec2.py() + elec4.py();
		  pzZ24 = elec2.pz() + elec4.pz();

		  // Calculate the momentum and invariant mass of elec 1 and 3, elec 2 and 4
		  pZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13) + (pzZ13 * pzZ13));
		  pZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24) + (pzZ24 * pzZ24));

		  pTZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13) + (pzZ13 * pzZ13));
		  pTZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24) + (pzZ24 * pzZ24));

		  mZ13 = sqrt((eZ13 * eZ13) - (pZ13 * pZ13));
		  mZ24 = sqrt((eZ24 * eZ24) - (pZ24 * pZ24));

		  // if (mZ13 > 0.) h_mZ13_4e->Fill(mZ13);
		  // if (mZ24 > 0.) h_mZ24_4e->Fill(mZ24);
		}
	    }

	  dZ13 = std::abs( mZ13 - mZ );
	  dZ24 = std::abs( mZ24 - mZ );

	  dZc2 = (dZ13 < dZ24) ? dZ13 : dZ24; 

	  // Third combination: Combine elec 1423
	  if (elec1.charge() + elec4.charge() == 0)
	    {
	      eZ14 = (sqrt(elec1.p() * elec1.p() + sqme)) + (sqrt(elec4.p() * elec4.p() + sqme));

	      pxZ14 = elec1.px() + elec4.px();
	      pyZ14 = elec1.py() + elec4.py();
	      pzZ14 = elec1.pz() + elec4.pz();

	      if (elec2.charge() + elec3.charge() == 0)
		{
		  eZ23 = sqrt((elec2.p() * elec2.p() + sqme)) + (sqrt(elec3.p() * elec3.p() + sqme));
       
		  pxZ23 = elec2.px() + elec3.px();
		  pyZ23 = elec2.py() + elec3.py();
		  pzZ23 = elec2.pz() + elec3.pz();

		  // Calculate the momentum and invariant mass of elec 1 and 4, elec 2 and 3
		  pZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14) + (pzZ14 * pzZ14));
		  pZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23) + (pzZ23 * pzZ23));

		  pTZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14) + (pzZ14 * pzZ14));
		  pTZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23) + (pzZ23 * pzZ23));

		  mZ14 = sqrt((eZ14 * eZ14) - (pZ14 * pZ14));
		  mZ23 = sqrt((eZ23 * eZ23) - (pZ23 * pZ23));

		  // if (mZ14 > 0.) h_mZ14_4e->Fill(mZ14);
		  // if (mZ23 > 0.) h_mZ23_4e->Fill(mZ23);
		}
	    }

	  dZ14 = std::abs( mZ14 - mZ );
	  dZ23 = std::abs( mZ23 - mZ );

	  dZc3 = (dZ14 < dZ23) ? dZ14 : dZ23; 

	  bool ptZadaug = false;

	  // Now whichever have the smallest diff is considered the best comb. 
	  if (dZc1 < dZc2 && dZc1 < dZc3)
	    {
	      if (dZ12 < dZ34)
		{
		  eZa  = eZ12;     
		  pxZa = pxZ12;
		  pyZa = pyZ12;
		  pzZa = pzZ12;
		  pTZa = pTZ12;
		  mZa  = mZ12;

		  if (elec1.pt() > 20. and elec2.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ34;
		  pxZb = pxZ34;
		  pyZb = pyZ34;
		  pzZb = pzZ34;
		  pTZb = pTZ34;
		  mZb  = mZ34;
		}
	      else
		{
		  eZa  = eZ34;  
		  pxZa = pxZ34;
		  pyZa = pyZ34;
		  pzZa = pzZ34;
		  pTZa = pTZ34;
		  mZa  = mZ34;

		  if (elec3.pt() > 20. and elec4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ12;
		  pxZb = pxZ12;
		  pyZb = pyZ12;
		  pzZb = pzZ12;
		  pTZb = pTZ12;
		  mZb  = mZ12;
		}
	    }

	  else if (dZc2 < dZc1 && dZc2 < dZc3)
	    {
	      if (dZ13 < dZ24)
		{
		  eZa  = eZ13;
		  pxZa = pxZ13;
		  pyZa = pyZ13;
		  pzZa = pzZ13;
		  pTZa = pTZ13;
		  mZa  = mZ13;

		  if (elec1.pt() > 20. and elec3.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ24;
		  pxZb = pxZ24;
		  pyZb = pyZ24;
		  pzZb = pzZ24;
		  pTZb = pTZ24;
		  mZb  = mZ24;
		}
	      else
		{
		  eZa  = eZ24;
		  pxZa = pxZ24;
		  pyZa = pyZ24;
		  pzZa = pzZ24;
		  pTZa = pTZ24;
		  mZa  = mZ24;

		  if (elec2.pt() > 20. and elec4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ13;
		  pxZb = pxZ13;
		  pyZb = pyZ13;
		  pzZb = pzZ13;
		  pTZb = pTZ13;
		  mZb  = mZ13;
		}
	    }

	  else if (dZc3 < dZc1 && dZc3 < dZc2)
	    {
	      if (dZ14 < dZ23)
		{
		  eZa  = eZ14;
		  pxZa = pxZ14;
		  pyZa = pyZ14;
		  pzZa = pzZ14;
		  pTZa = pTZ14;
		  mZa  = mZ14;

		  if (elec1.pt() > 20. and elec4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ23;
		  pxZb = pxZ23;
		  pyZb = pyZ23;
		  pzZb = pzZ23;
		  pTZb = pTZ23;
		  mZb  = mZ23;
		}
	      else
		{
		  eZa  = eZ23;
		  pxZa = pxZ23;
		  pyZa = pyZ23;
		  pzZa = pzZ23;
		  pTZa = pTZ23;
		  mZa  = mZ23;

		  if (elec2.pt() > 20. and elec3.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ14;
		  pxZb = pxZ14;
		  pyZb = pyZ14;
		  pzZb = pzZ14;
		  pTZb = pTZ14;
		  mZb  = mZ14;
		}
	    }

          //value_el_pt[value_el_n] = it->pt();
          //value_el_eta[value_el_n] = it->eta(); // USE THEIR SYNTAX Done
          //value_el_phi[value_el_n] = it->phi(); // USE THEIR SYNTAX Done
          //value_el_charge[value_el_n] = it->charge();
          //value_el_rawE[value_el_n] = it->superCluster()->rawEnergy(); // USE THEIR SYNTAX
          //value_el_misshits[value_el_n] = it->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
          //auto iso03 = it->pfIsolationVariables();
          //value_el_pfreliso[value_el_n] = (iso03.chargedHadronIso + iso03.neutralHadronIso + iso03.photonIso)/it->pt();
          //auto trk = it->gsfTrack();
          //value_el_dxy[value_el_n] = trk->dxy(pv);
          //value_el_dz[value_el_n] = trk->dz(pv);
          //float d3 = sqrt( pow(trk->dxy(pv),2) + pow(trk->dz(pv),2) );
          //float d3Err = sqrt( pow(trk->d0Error(),2) + pow(trk->dzError(),2) );
          //value_el_SIP3d[value_el_n] = d3/d3Err;
      
	  value_el_n = 4;
    
	  value_el_pt[0] = elec1.pt();
	  value_el_pt[1] = elec2.pt();
	  value_el_pt[2] = elec3.pt();
	  value_el_pt[3] = elec4.pt();

	  value_el_eta[0] = elec1.eta();
	  value_el_eta[1] = elec2.eta();
	  value_el_eta[2] = elec3.eta();
	  value_el_eta[3] = elec4.eta();

	  value_el_phi[0] = elec1.phi();
	  value_el_phi[1] = elec2.phi();
	  value_el_phi[2] = elec3.phi();
	  value_el_phi[3] = elec4.phi();
      
	  value_el_charge[0] = elec1.charge();
	  value_el_charge[1] = elec2.charge();
	  value_el_charge[2] = elec3.charge();
	  value_el_charge[3] = elec4.charge();

	  value_el_rawE[0] = elec1.superCluster()->rawEnergy();
	  value_el_rawE[1] = elec2.superCluster()->rawEnergy();
	  value_el_rawE[2] = elec3.superCluster()->rawEnergy();
	  value_el_rawE[3] = elec4.superCluster()->rawEnergy();

	  value_el_misshits[0] = ((elec1.gsfTrack())->trackerExpectedHitsInner()).numberOfHits();
	  value_el_misshits[1] = ((elec2.gsfTrack())->trackerExpectedHitsInner()).numberOfHits();
	  value_el_misshits[2] = ((elec3.gsfTrack())->trackerExpectedHitsInner()).numberOfHits();
	  value_el_misshits[3] = ((elec4.gsfTrack())->trackerExpectedHitsInner()).numberOfHits();

	  value_el_pfreliso[0] = ((elec1.pfIsolationVariables()).chargedHadronIso + (elec1.pfIsolationVariables()).neutralHadronIso + (elec1.pfIsolationVariables()).photonIso) / (elec1.pt());
	  value_el_pfreliso[1] = ((elec2.pfIsolationVariables()).chargedHadronIso + (elec2.pfIsolationVariables()).neutralHadronIso + (elec2.pfIsolationVariables()).photonIso) / (elec2.pt());
	  value_el_pfreliso[2] = ((elec3.pfIsolationVariables()).chargedHadronIso + (elec3.pfIsolationVariables()).neutralHadronIso + (elec3.pfIsolationVariables()).photonIso) / (elec3.pt());
	  value_el_pfreliso[3] = ((elec4.pfIsolationVariables()).chargedHadronIso + (elec4.pfIsolationVariables()).neutralHadronIso + (elec4.pfIsolationVariables()).photonIso) / (elec4.pt());

	  value_el_dxy[0] = elec1.gsfTrack()->dxy(pv);
	  value_el_dxy[1] = elec2.gsfTrack()->dxy(pv);
	  value_el_dxy[2] = elec3.gsfTrack()->dxy(pv);
	  value_el_dxy[3] = elec4.gsfTrack()->dxy(pv);

	  value_el_dz[0] = elec1.gsfTrack()->dz(pv);
	  value_el_dz[1] = elec2.gsfTrack()->dz(pv);
	  value_el_dz[2] = elec3.gsfTrack()->dz(pv);
	  value_el_dz[3] = elec4.gsfTrack()->dz(pv);

	  value_el_SIP3d[0] = (sqrt(pow(elec1.gsfTrack()->dxy(pv),2) + pow(elec1.gsfTrack()->dz(pv),2)))/(sqrt(pow(elec1.gsfTrack()->d0Error(),2) + pow(elec1.gsfTrack()->dzError(),2)));
	  value_el_SIP3d[1] = (sqrt(pow(elec2.gsfTrack()->dxy(pv),2) + pow(elec2.gsfTrack()->dz(pv),2)))/(sqrt(pow(elec2.gsfTrack()->d0Error(),2) + pow(elec2.gsfTrack()->dzError(),2)));
	  value_el_SIP3d[2] = (sqrt(pow(elec3.gsfTrack()->dxy(pv),2) + pow(elec3.gsfTrack()->dz(pv),2)))/(sqrt(pow(elec3.gsfTrack()->d0Error(),2) + pow(elec3.gsfTrack()->dzError(),2)));
	  value_el_SIP3d[3] = (sqrt(pow(elec4.gsfTrack()->dxy(pv),2) + pow(elec4.gsfTrack()->dz(pv),2)))/(sqrt(pow(elec4.gsfTrack()->d0Error(),2) + pow(elec4.gsfTrack()->dzError(),2)));

	  elecTree = true;


	  if (ptZadaug) {
	    if (mZa > 40. && mZa < 120.) {
	      if (mZb > 12. && mZb < 120.) {
		// h_mZa_4e->Fill(mZa);
		// h_mZb_4e->Fill(mZb);
		
		// Calculate 4 elec
		p4Za.SetPxPyPzE(pxZa, pyZa, pzZa, eZa);
		p4Zb.SetPxPyPzE(pxZb, pyZb, pzZb, eZb);

		p4H = p4Za + p4Zb;

		mass4e = p4H.M();
		pt_4e = p4H.Pt();
		eta_4e = p4H.Eta();
		phi_4e = p4H.Phi();
		
		px4e = p4H.Px();
		py4e = p4H.Py();
		pz4e = p4H.Pz();
		E4e = p4H.E();		

		pt_e1 = elec1.pt();
		pt_e2 = elec2.pt();
		pt_e3 = elec3.pt();
		pt_e4 = elec4.pt();

		eta_e1 = elec1.eta();
		eta_e2 = elec2.eta();
		eta_e3 = elec3.eta();
		eta_e4 = elec4.eta();

		phi_e1 = elec1.phi();
		phi_e2 = elec2.phi();
		phi_e3 = elec3.phi();
		phi_e4 = elec4.phi();

		cas_e1 = elec1.charge();
		cas_e2 = elec2.charge();
		cas_e3 = elec3.charge();
		cas_e4 = elec4.charge();

		px_e1 = elec1.px();
		px_e2 = elec2.px();
		px_e3 = elec3.px();
		px_e4 = elec4.px();

		py_e1 = elec1.py();
		py_e2 = elec2.py();
		py_e3 = elec3.py();
		py_e4 = elec4.py();
		
		pz_e1 = elec1.pz();
		pz_e2 = elec2.pz();
		pz_e3 = elec3.pz();
		pz_e4 = elec4.pz();

		E_e1 = elec1.energy();
		E_e2 = elec2.energy();
		E_e3 = elec3.energy();
		E_e4 = elec4.energy();

		myfile4e << "thing two" << endl;

		if (mass4e > 70.)
		  {
		    /*h_m1_m4e->Fill(mass4e);
		    h_m2_m4e->Fill(mass4e);
		    h_m3_m4e->Fill(mass4e);
		    h_m4_m4e->Fill(mass4e);*/
        
      		if (E_mu1 != -999) {
        		if (myfile4e.is_open()){
				      myfile4e << "thing one" << endl;
					  myfile4e << E_e1 << ", " << E_e2 << ", " << E_e3 << ", " << E_e4 << ", " << px_e1 << ", " << px_e2 << ", " << px_e3 << ", " << px_e4 << ", " << py_e1 << ", " << py_e2 << ", " << py_e3 << ", " << py_e4 << ", " << pz_e1 << ", " << pz_e2 << ", " << pz_e3 << ", " << pz_e4 << "\n";
          		}
    		}

		    for (unsigned i = 0; i < vIdPte.size(); i++)
		      {
			relPFIso_e = ((((*electrons)[vIdPte.at(i).first]).pfIsolationVariables()).chargedHadronIso +
				      (((*electrons)[vIdPte.at(i).first]).pfIsolationVariables()).neutralHadronIso +
				      (((*electrons)[vIdPte.at(i).first]).pfIsolationVariables()).photonIso) / (((*electrons)[vIdPte.at(i).first]).pt());

			// h_relPFIso_e_after->Fill(relPFIso_e);

			// h_pt_e_after->Fill(vIdPte.at(i).second);
			// h_eta_e_after->Fill((((*electrons)[vIdPte.at(i).first]).superCluster())->eta());
		      }
		    // t2->Fill();
		  }
	      } // end of mZb 
	    } // end of mZa
	  } // end of ptZdaug
	} // end of total charge
    } // end of nGoodElectron

  //=============================== ZZ/ZZ*To4e end =============================//


  //=========================== ZZ/ZZ*To2mu2e start ============================//

  if (nGoodRecoMuon >= 2 && nGoodElectron >= 2)
    {
      const reco::Muon &muon1 = (*muons)[vIdPtmu.at(0).first];
      const reco::Muon &muon2 = (*muons)[vIdPtmu.at(1).first];
      const reco::GsfElectron &elec1 = (*electrons)[vIdPte.at(0).first];
      const reco::GsfElectron &elec2 = (*electrons)[vIdPte.at(1).first];

      if (muon1.charge() + muon2.charge() + elec1.charge() + elec2.charge() == 0)
	{
	  // For case 2mu2e, there is only 1 combination
	  if (muon1.charge() + muon2.charge() == 0)
	    {
	      eZ12 = (sqrt(muon1.p() * muon1.p() + sqm1)) + (sqrt(muon2.p() * muon2.p() + sqm1));
	
	      pxZ12 = muon1.px() + muon2.px();
	      pyZ12 = muon1.py() + muon2.py();
	      pzZ12 = muon1.pz() + muon2.pz();

	      if (elec1.charge() + elec2.charge() == 0)
		{
		  eZ34 = (sqrt(elec1.p() * elec1.p() + sqme)) + (sqrt(elec2.p() * elec2.p() + sqme));

		  pxZ34 = elec1.px() + elec2.px();
		  pyZ34 = elec1.py() + elec2.py();
		  pzZ34 = elec1.pz() + elec2.pz();

		  // Calculate the momentum and invariant mass of muon 12, elec 34
		  pZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12) + (pzZ12 * pzZ12));
		  pZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34) + (pzZ34 * pzZ34));

		  pTZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12));
		  pTZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34));

		  mZ12 = sqrt((eZ12 * eZ12) - (pZ12 * pZ12));
		  mZ34 = sqrt((eZ34 * eZ34) - (pZ34 * pZ34));

		  // if (mZ12 > 0.) h_mZmu_2mu2e->Fill(mZ12);
		  // if (mZ34 > 0.) h_mZe_2mu2e->Fill(mZ34);

		}
	    }

	  dZ12 = std::abs(mZ12 - mZ); // mu
	  dZ34 = std::abs(mZ34 - mZ); // e

	  bool ptZadaug = false;

	  if (dZ12 < dZ34)
	    {
	      eZa  = eZ12;
	      pxZa = pxZ12;
	      pyZa = pyZ12;
	      pzZa = pzZ12;
	      pTZa = pTZ12;
	      mZa  = mZ12;

	      if (muon1.pt() > 20. and muon2.pt() > 10.)
		ptZadaug = true;

	      eZb  = eZ34;
	      pxZb = pxZ34;
	      pyZb = pyZ34;
	      pzZb = pzZ34;
	      pTZb = pTZ34;
	      mZb  = mZ34;
	    }
	  else
	    {
	      eZa  = eZ34;
	      pxZa = pxZ34;
	      pyZa = pyZ34;
	      pzZa = pzZ34;
	      pTZa = pTZ34;
	      mZa  = mZ34;

	      if (elec1.pt() > 20. and elec2.pt() > 10.)
		ptZadaug = true;

	      eZb  = eZ12;
	      pxZb = pxZ12;
	      pyZb = pyZ12;
	      pzZb = pzZ12;
	      pTZb = pTZ12;
	      mZb  = mZ12;
	  }
  
    

    if (elecTree == false && muonTree == false) {

      value_el_n = 2;
      value_mu_n = 2;
      value_el_pt[0] = elec1.pt();
      value_el_pt[1] = elec2.pt();

      value_el_eta[0] = elec1.eta();
		  value_el_eta[1] = elec2.eta();

      value_el_phi[0] = elec1.phi();
		  value_el_phi[1] = elec2.phi();
      
      value_el_charge[0] = elec1.charge();
		  value_el_charge[1] = elec2.charge();

      value_el_rawE[0] = elec1.superCluster()->rawEnergy();
      value_el_rawE[1] = elec2.superCluster()->rawEnergy();

      value_el_misshits[0] = ((elec1.gsfTrack())->trackerExpectedHitsInner()).numberOfHits();
      value_el_misshits[1] = ((elec2.gsfTrack())->trackerExpectedHitsInner()).numberOfHits();

      value_el_pfreliso[0] = ((elec1.pfIsolationVariables()).chargedHadronIso + (elec1.pfIsolationVariables()).neutralHadronIso + (elec1.pfIsolationVariables()).photonIso) / (elec1.pt());
      value_el_pfreliso[1] = ((elec2.pfIsolationVariables()).chargedHadronIso + (elec2.pfIsolationVariables()).neutralHadronIso + (elec2.pfIsolationVariables()).photonIso) / (elec2.pt());
      
      value_el_dxy[0] = elec1.gsfTrack()->dxy(pv);
      value_el_dxy[1] = elec2.gsfTrack()->dxy(pv);
      
      value_el_dz[0] = elec1.gsfTrack()->dz(pv);
      value_el_dz[1] = elec2.gsfTrack()->dz(pv);
      
      value_el_SIP3d[0] = (sqrt(pow(elec1.gsfTrack()->dxy(pv),2) + pow(elec1.gsfTrack()->dz(pv),2)))/(sqrt(pow(elec1.gsfTrack()->d0Error(),2) + pow(elec1.gsfTrack()->dzError(),2)));
      value_el_SIP3d[1] = (sqrt(pow(elec2.gsfTrack()->dxy(pv),2) + pow(elec2.gsfTrack()->dz(pv),2)))/(sqrt(pow(elec2.gsfTrack()->d0Error(),2) + pow(elec2.gsfTrack()->dzError(),2)));
    
      value_mu_pt[0] = muon1.pt();
      value_mu_pt[1] = muon2.pt();

      //myfile << muon1.pt() << ", " << muon2.pt() << ", " << value_mu_pt[0] << ", " << value_mu_pt[1] << endl;

      value_mu_eta[0] = muon1.eta();
		  value_mu_eta[1] = muon2.eta();

      //myfile4mu << muon1.eta() << ", " << muon2.eta() << ", " << value_mu_eta[0] << ", " << value_mu_eta[1] << endl;
		
      value_mu_phi[0] = muon1.phi();
		  value_mu_phi[1] = muon2.phi();
		
		  value_mu_charge[0] = muon1.charge();
		  value_mu_charge[1] = muon2.charge();
		
      value_mu_mass[0] = mass2mu2e;
		  value_mu_mass[1] = mass2mu2e;
		
    //I'm not sure about this one
      value_mu_isGM[0] = muon1.isGlobalMuon();
		  value_mu_isGM[1] = muon2.isGlobalMuon();
		
    //Or this one
      value_mu_pfrelisoDB[0] = ((muon1.pfIsolationR04()).sumChargedHadronPt + (muon1.pfIsolationR04()).sumNeutralHadronEt + (muon1.pfIsolationR04()).sumPhotonEt - 0.5*(muon1.pfIsolationR04()).sumPUPt) / (muon1.pt());
		  value_mu_pfrelisoDB[1] = ((muon2.pfIsolationR04()).sumChargedHadronPt + (muon2.pfIsolationR04()).sumNeutralHadronEt + (muon2.pfIsolationR04()).sumPhotonEt - 0.5*(muon2.pfIsolationR04()).sumPUPt) / (muon2.pt());
		
    //This one as well
      value_mu_pfreliso[0] = ((muon1.pfIsolationR04()).sumChargedHadronPt + (muon1.pfIsolationR04()).sumNeutralHadronEt + (muon1.pfIsolationR04()).sumPhotonEt) / (muon1.pt());
		  value_mu_pfreliso[1] = ((muon2.pfIsolationR04()).sumChargedHadronPt + (muon2.pfIsolationR04()).sumNeutralHadronEt + (muon2.pfIsolationR04()).sumPhotonEt) / (muon2.pt());
		
      value_mu_dxy[0] = muon1.globalTrack()->dxy(pv);
		  value_mu_dxy[1] = muon2.globalTrack()->dxy(pv);
		
      value_mu_dz[0] = muon1.globalTrack()->dz(pv);
		  value_mu_dz[1] = muon2.globalTrack()->dz(pv);
		
      value_mu_SIP3d[0] = (sqrt(pow(muon1.globalTrack()->dxy(pv),2) + pow(muon1.globalTrack()->dz(pv),2)))/(sqrt(pow(muon1.globalTrack()->d0Error(),2) + pow(muon1.globalTrack()->dzError(),2)));
      value_mu_SIP3d[1] = (sqrt(pow(muon2.globalTrack()->dxy(pv),2) + pow(muon2.globalTrack()->dz(pv),2)))/(sqrt(pow(muon2.globalTrack()->d0Error(),2) + pow(muon2.globalTrack()->dzError(),2)));
      //cout << "end of the tree input inside of the loop" << endl;
    
    
    }
    //cout << "end of the tree input outside of the loop" << endl;


	  if (ptZadaug) {
	    if (mZa > 40. && mZa < 120.) {
	      if (mZb > 12. && mZb < 120.) {
		// h_mZa_2mu2e->Fill(mZa);
		// h_mZb_2mu2e->Fill(mZb);

		// Now combine these 2 muons and 2 electrons
 
		// Calculate 4 lepton: 2muons 2electrons

		p4Za.SetPxPyPzE(pxZa, pyZa, pzZa, eZa);
		p4Zb.SetPxPyPzE(pxZb, pyZb, pzZb, eZb);

		p4H = p4Za + p4Zb;

		mass2mu2e = p4H.M();
		pt_2mu2e = p4H.Pt();
		eta_2mu2e = p4H.Eta();
		phi_2mu2e = p4H.Phi();
		
		px2mu2e = p4H.Px();
		py2mu2e = p4H.Py();
		pz2mu2e = p4H.Pz();
		E2mu2e = p4H.E();

		pt_2mu1 = muon1.pt();
		pt_2mu2 = muon2.pt();
		pt_2e1 = elec1.pt();
		pt_2e2 = elec2.pt();

		eta_2mu1 = muon1.eta();
		eta_2mu2 = muon2.eta();
		eta_2e1 = elec1.eta();
		eta_2e2 = elec2.eta();

		phi_2mu1 = muon1.phi();
		phi_2mu2 = muon2.phi();
		phi_2e1 = elec1.phi();
		phi_2e2 = elec2.phi();

		cas_2mu1 = muon1.charge();
		cas_2mu2 = muon2.charge();
		cas_2e1 = elec1.charge();
		cas_2e2 = elec2.charge();

		px_2mu1 = muon1.px();
		px_2mu2 = muon2.px();
		px_2e1 = elec1.px();
		px_2e2 = elec2.px();

		py_2mu1 = muon1.py();
		py_2mu2 = muon2.py();
		py_2e1 = elec1.py();
		py_2e2 = elec2.py();
		
		pz_2mu1 = muon1.pz();
		pz_2mu2 = muon2.pz();
		pz_2e1 = elec1.pz();
		pz_2e2 = elec2.pz();

		E_2mu1 = sqrt((muon1.p() * muon1.p()) + sqm1);
		E_2mu2 = sqrt((muon2.p() * muon2.p()) + sqm1);
		E_2e1 = elec1.energy();
		E_2e2 = elec2.energy();

		if (mass2mu2e > 70.)
		  {
		    /*h_m1_m2mu2e->Fill(mass2mu2e);
		    h_m2_m2mu2e->Fill(mass2mu2e);
		    h_m3_m2mu2e->Fill(mass2mu2e);
		    h_m4_m2mu2e->Fill(mass2mu2e);*/

        cout << "filling the txt file" << endl;
        if (E_2mu1 != -999) {
		      if (myfile2mu2e.is_open()){
			      myfile2mu2e << E_2mu1 << ", " << E_2mu2 << ", " << E_2e1 << ", " << E_2e2 << ", " << px_2mu1 << ", " << px_2mu2 << ", " << px_2e1 << ", " << px_2e2 << ", " << py_2mu1 << ", " << py_2mu2 << ", " << py_2e1 << ", " << py_2e2 << ", " << pz_2mu1 << ", " << pz_2mu2 << ", " << pz_2e1 << ", " << pz_2e2 << "\n";
          }
        }

		    for (unsigned i = 0; i < vIdPtmu.size(); i++)
		      {
			relPFIso_mu = ( (((*muons)[vIdPtmu.at(i).first]).pfIsolationR04()).sumChargedHadronPt +
					(((*muons)[vIdPtmu.at(i).first]).pfIsolationR04()).sumNeutralHadronEt +
					(((*muons)[vIdPtmu.at(i).first]).pfIsolationR04()).sumPhotonEt ) / (((*muons)[vIdPtmu.at(i).first]).pt());

			/*h_relPFIso_2mu_after->Fill(relPFIso_mu);
			h_pt_after_2mu2e->Fill(vIdPtmu.at(i).second);
			h_eta_after_2mu2e->Fill(((*muons)[vIdPtmu.at(i).first]).eta());*/
		      }

		    for (unsigned i = 0; i < vIdPte.size(); i++)
		      {
			relPFIso_e = ((((*electrons)[vIdPte.at(i).first]).pfIsolationVariables()).chargedHadronIso +
				      (((*electrons)[vIdPte.at(i).first]).pfIsolationVariables()).neutralHadronIso +
				      (((*electrons)[vIdPte.at(i).first]).pfIsolationVariables()).photonIso) / (((*electrons)[vIdPte.at(i).first]).pt()); 
	      
			/*h_relPFIso_2e_after->Fill(relPFIso_e);
			h_pt_e_after_2mu2e->Fill(vIdPte.at(i).second);
			h_eta_e_after_2mu2e->Fill(((*electrons)[vIdPte.at(i).first]).eta());*/
		      }
		    // t3->Fill();
		  }
	      }
	    }
	  }
	} // end of total charge
    } // end of nGoodMuonNElectron

  //============================= ZZ/ZZ*To2mu2e end ============================//
  
  // Fill event
  tree->Fill();

} // HiggsDemoAnalyzer::analyze ends


// ------ method called once each job just before starting event loop ---------//

void HiggsDemoAnalyzer::beginJob() {

  // *******************************************************
  // book the ntuple for the surviving 4 lepton candidates *
  // in the mass range 70 < m4l < 181 GeV                  *
  // *******************************************************


    myfile4mu.open ("4mu.txt", std::ios_base::app);
    myfile4e.open ("4e.txt", std::ios_base::app);
    myfile2mu2e.open ("2mu2e.txt", std::ios_base::app);

/*
  t1 = new TTree("tree4mu", "tree4mu");
  t2 = new TTree("tree4e", "tree4e");
  t3 = new TTree("tree2mu2e", "tree2mu2e");

  // tree 4mu
  t1->Branch("nRun", &nRun, "nRun/I");
  t1->Branch("nEvt", &nEvt, "nEvt/I");
  t1->Branch("nLumi", &nLumi, "nLumi/I");
  t1->Branch("mass4mu", &mass4mu, "mass4mu/D");
  t1->Branch("pt_4mu", &pt_4mu, "pt_4mu/D");
  t1->Branch("eta_4mu", &eta_4mu, "eta_4mu/D");
  t1->Branch("phi_4mu", &phi_4mu, "phi_4mu/D");

  t1->Branch("px4mu", &px4mu, "px4mu/D");
  t1->Branch("py4mu", &py4mu, "py4mu/D");
  t1->Branch("pz4mu", &pz4mu, "pz4mu/D");
  t1->Branch("E4mu", &E4mu, "E4mu/D");
  
  t1->Branch("pt_mu1", &pt_mu1, "pt_mu1/D");
  t1->Branch("pt_mu2", &pt_mu2, "pt_mu2/D");
  t1->Branch("pt_mu3", &pt_mu3, "pt_mu3/D");
  t1->Branch("pt_mu4", &pt_mu4, "pt_mu4/D");
  t1->Branch("eta_mu1", &eta_mu1, "eta_mu1/D");
  t1->Branch("eta_mu2", &eta_mu2, "eta_mu2/D");
  t1->Branch("eta_mu3", &eta_mu3, "eta_mu3/D");
  t1->Branch("eta_mu4", &eta_mu4, "eta_mu4/D");

  t1->Branch("phi_mu1", &phi_mu1, "phi_mu1/D");
  t1->Branch("phi_mu2", &phi_mu2, "phi_mu2/D");
  t1->Branch("phi_mu3", &phi_mu3, "phi_mu3/D");
  t1->Branch("phi_mu4", &phi_mu4, "phi_mu4/D");

  t1->Branch("cas_mu1", &cas_mu1, "cas_mu1/I");
  t1->Branch("cas_mu2", &cas_mu2, "cas_mu2/I");
  t1->Branch("cas_mu3", &cas_mu3, "cas_mu3/I");
  t1->Branch("cas_mu4", &cas_mu4, "cas_mu4/I");

  t1->Branch("px_mu1", &px_mu1, "px_mu1/D");
  t1->Branch("px_mu2", &px_mu2, "px_mu2/D");
  t1->Branch("px_mu3", &px_mu3, "px_mu3/D");
  t1->Branch("px_mu4", &px_mu4, "px_mu4/D");

  t1->Branch("py_mu1", &py_mu1, "py_mu1/D");
  t1->Branch("py_mu2", &py_mu2, "py_mu2/D");
  t1->Branch("py_mu3", &py_mu3, "py_mu3/D");
  t1->Branch("py_mu4", &py_mu4, "py_mu4/D");

  t1->Branch("pz_mu1", &pz_mu1, "pz_mu1/D");
  t1->Branch("pz_mu2", &pz_mu2, "pz_mu2/D");
  t1->Branch("pz_mu3", &pz_mu3, "pz_mu3/D");
  t1->Branch("pz_mu4", &pz_mu4, "pz_mu4/D");

  t1->Branch("E_mu1", &E_mu1, "E_mu1/D");
  t1->Branch("E_mu2", &E_mu2, "E_mu2/D");
  t1->Branch("E_mu3", &E_mu3, "E_mu3/D");
  t1->Branch("E_mu4", &E_mu4, "E_mu4/D");

  t1->Branch("mZa", &mZa, "mZa/D");
  t1->Branch("mZb", &mZb, "mZb/D");
  
  // tree 4e
  t2->Branch("nRun", &nRun, "nRun/I");
  t2->Branch("nEvt", &nEvt, "nEvt/I");
  t2->Branch("nLumi", &nLumi, "nLumi/I");
  t2->Branch("mass4e", &mass4e, "mass4e/D");
  t2->Branch("pt_4e", &pt_4e, "pt_4e/D");
  t2->Branch("eta_4e", &eta_4e, "eta_4e/D");
  t2->Branch("phi_4e", &phi_4e, "phi_4e/D");

  t2->Branch("px4e", &px4e, "px4e/D");
  t2->Branch("py4e", &py4e, "py4e/D");
  t2->Branch("pz4e", &pz4e, "pz4e/D");
  t2->Branch("E4e", &E4e, "E4e/D");
  
  t2->Branch("pt_e1", &pt_e1, "pt_e1/D");
  t2->Branch("pt_e2", &pt_e2, "pt_e2/D");
  t2->Branch("pt_e3", &pt_e3, "pt_e3/D");
  t2->Branch("pt_e4", &pt_e4, "pt_e4/D");
  t2->Branch("eta_e1", &eta_e1, "eta_e1/D");
  t2->Branch("eta_e2", &eta_e2, "eta_e2/D");
  t2->Branch("eta_e3", &eta_e3, "eta_e3/D");
  t2->Branch("eta_e4", &eta_e4, "eta_e4/D");

  t2->Branch("phi_e1", &phi_e1, "phi_e1/D");
  t2->Branch("phi_e2", &phi_e2, "phi_e2/D");
  t2->Branch("phi_e3", &phi_e3, "phi_e3/D");
  t2->Branch("phi_e4", &phi_e4, "phi_e4/D");

  t2->Branch("cas_e1", &cas_e1, "cas_e1/I");
  t2->Branch("cas_e2", &cas_e2, "cas_e2/I");
  t2->Branch("cas_e3", &cas_e3, "cas_e3/I");
  t2->Branch("cas_e4", &cas_e4, "cas_e4/I");

  t2->Branch("px_e1", &px_e1, "px_e1/D");
  t2->Branch("px_e2", &px_e2, "px_e2/D");
  t2->Branch("px_e3", &px_e3, "px_e3/D");
  t2->Branch("px_e4", &px_e4, "px_e4/D");

  t2->Branch("py_e1", &py_e1, "py_e1/D");
  t2->Branch("py_e2", &py_e2, "py_e2/D");
  t2->Branch("py_e3", &py_e3, "py_e3/D");
  t2->Branch("py_e4", &py_e4, "py_e4/D");

  t2->Branch("pz_e1", &pz_e1, "pz_e1/D");
  t2->Branch("pz_e2", &pz_e2, "pz_e2/D");
  t2->Branch("pz_e3", &pz_e3, "pz_e3/D");
  t2->Branch("pz_e4", &pz_e4, "pz_e4/D");

  t2->Branch("E_e1", &E_e1, "E_e1/D");
  t2->Branch("E_e2", &E_e2, "E_e2/D");
  t2->Branch("E_e3", &E_e3, "E_e3/D");
  t2->Branch("E_e4", &E_e4, "E_e4/D");

  t2->Branch("mZa", &mZa, "mZa/D");
  t2->Branch("mZb", &mZb, "mZb/D");
  
  // tree 2mu 2e
  t3->Branch("nRun", &nRun, "nRun/I");
  t3->Branch("nEvt", &nEvt, "nEvt/I");
  t3->Branch("nLumi", &nLumi, "nLumi/I");
  t3->Branch("mass2mu2e", &mass2mu2e, "mass2mu2e/D");
  t3->Branch("pt_2mu2e", &pt_2mu2e, "pt_2mu2e/D");
  t3->Branch("eta_2mu2e", &eta_2mu2e, "eta_2mu2e/D");
  t3->Branch("phi_2mu2e", &phi_2mu2e, "phi_2mu2e/D");

  t3->Branch("px2mu2e", &px2mu2e, "px2mu2e/D");
  t3->Branch("py2mu2e", &py2mu2e, "py2mu2e/D");
  t3->Branch("pz2mu2e", &pz2mu2e, "pz2mu2e/D");
  t3->Branch("E2mu2e", &E2mu2e, "E2mu2e/D");
  
  t3->Branch("pt_2mu1", &pt_2mu1, "pt_2mu1/D");
  t3->Branch("pt_2mu2", &pt_2mu2, "pt_2mu2/D");
  t3->Branch("pt_2e1", &pt_2e1, "pt_2e1/D");
  t3->Branch("pt_2e2", &pt_2e2, "pt_2e2/D");
  t3->Branch("eta_2mu1", &eta_2mu1, "eta_2mu1/D");
  t3->Branch("eta_2mu2", &eta_2mu2, "eta_2mu2/D");
  t3->Branch("eta_2e1", &eta_2e1, "eta_2e1/D");
  t3->Branch("eta_2e2", &eta_2e2, "eta_2e2/D");

  t3->Branch("phi_2mu1", &phi_2mu1, "phi_2mu1/D");
  t3->Branch("phi_2mu2", &phi_2mu2, "phi_2mu2/D");
  t3->Branch("phi_2e1", &phi_2e1, "phi_2e1/D");
  t3->Branch("phi_2e2", &phi_2e2, "phi_2e2/D");

  t3->Branch("cas_2mu1", &cas_2mu1, "cas_2mu1/I");
  t3->Branch("cas_2mu2", &cas_2mu2, "cas_2mu2/I");
  t3->Branch("cas_2e1", &cas_2e1, "cas_2e1/I");
  t3->Branch("cas_2e2", &cas_2e2, "cas_2e2/I");

  t3->Branch("px_2mu1", &px_2mu1, "px_2mu1/D");
  t3->Branch("px_2mu2", &px_2mu2, "px_2mu2/D");
  t3->Branch("px_2e1", &px_2e1, "px_2e1/D");
  t3->Branch("px_2e2", &px_2e2, "px_2e2/D");

  t3->Branch("py_2mu1", &py_2mu1, "py_2mu1/D");
  t3->Branch("py_2mu2", &py_2mu2, "py_2mu2/D");
  t3->Branch("py_2e1", &py_2e1, "py_2e1/D");
  t3->Branch("py_2e2", &py_2e2, "py_2e2/D");

  t3->Branch("pz_2mu1", &pz_2mu1, "pz_2mu1/D");
  t3->Branch("pz_2mu2", &pz_2mu2, "pz_2mu2/D");
  t3->Branch("pz_2e1", &pz_2e1, "pz_2e1/D");
  t3->Branch("pz_2e2", &pz_2e2, "pz_2e2/D");

  t3->Branch("E_2mu1", &E_2mu1, "E_2mu1/D");
  t3->Branch("E_2mu2", &E_2mu2, "E_2mu2/D");
  t3->Branch("E_2e1", &E_2e1, "E_2e1/D");
  t3->Branch("E_2e2", &E_2e2, "E_2e2/D");

  t3->Branch("mZa", &mZa, "mZa/D");
  t3->Branch("mZb", &mZb, "mZb/D");
  */
  
}

// ------------ method called once each job just after ending the event loop  ------------
void HiggsDemoAnalyzer::endJob() {
  myfile4mu.close();
  myfile4e.close();
  myfile2mu2e.close();
  myfile.close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiggsDemoAnalyzer);

