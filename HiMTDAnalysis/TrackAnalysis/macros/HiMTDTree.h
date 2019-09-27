//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 23 13:04:52 2019 by ROOT version 6.06/01
// from TTree HiMTDTree/
// found on file: /storage1/users/wl33/MTD/HiMTDTree_numEvent1000.root
//////////////////////////////////////////////////////////

#ifndef HiMTDTree_h
#define HiMTDTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TVector3.h"
#include "vector"
#include "vector"

class HiMTDTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          Event_Run;
   UShort_t        Event_Lumi;
   UInt_t          Event_Bx;
   ULong64_t       Event_Orbit;
   ULong64_t       Event_Number;
   UChar_t         Event_nPV;
   TVector3        *Event_PriVtx_Pos;
   TVector3        *Event_PriVtx_Err;
   Float_t         Event_hiHF;
   Int_t           Event_hiBin;
   UInt_t          Reco_Track_N;
   vector<float>   *Reco_Track_beta_PI;
   vector<float>   *Reco_Track_t0_PI;
   vector<float>   *Reco_Track_t0Err_PI;
   vector<float>   *Reco_Track_beta_PV;
   vector<float>   *Reco_Track_t0_PV;
   vector<float>   *Reco_Track_t0Err_PV;
   vector<float>   *Reco_Track_massSq_PV;
   vector<float>   *Reco_Track_beta_PID;
   vector<float>   *Reco_Track_t0_PID;
   vector<float>   *Reco_Track_t0Err_PID;
   vector<float>   *Reco_Track_massSq_PID;
   vector<float>   *Reco_Track_pathLength;
   vector<float>   *Reco_Track_p;
   vector<float>   *Reco_Track_tMTD;
   vector<float>   *Reco_Track_tMTDErr;
   vector<float>   *Reco_Track_eta;
   vector<float>   *Reco_Track_phi;
   vector<float>   *Reco_Track_charge;
   vector<int>     *Reco_Track_quality;
   vector<float>   *Reco_Track_pt;
   vector<float>   *Reco_Track_ptErr;
   vector<float>   *Reco_Track_dcaDz;
   vector<float>   *Reco_Track_dcaDxy;
   vector<float>   *Reco_Track_normChi2;
   vector<int>     *Reco_Track_nHits;
   vector<float>   *Reco_Track_dEdx;
   vector<float>   *Reco_Track_dEdxErr;
   vector<int>     *Reco_Track_nSatMea;
   vector<int>     *Reco_Track_nMea;
   vector<float>   *Gen_Track_pdgId;
   vector<float>   *Gen_Track_pt;
   vector<float>   *Gen_Track_eta;
   vector<float>   *Gen_Track_phi;
   vector<float>   *Gen_Track_t0;
   vector<float>   *Gen_Track_charge;
   vector<float>   *Gen_Track_mass;
   vector<float>   *Gen_Track_momPdgId;
   UInt_t          Reco_DiTrack_N;
   vector<int>     *Reco_DiTrack_Pion_idx;
   vector<int>     *Reco_DiTrack_Kaon_idx;
   vector<float>   *Reco_DiTrack_pt;
   vector<float>   *Reco_DiTrack_rap;
   vector<float>   *Reco_DiTrack_mass;
   vector<float>   *Reco_DiTrack_vProb;
   vector<float>   *Reco_DiTrack_alpha;
   vector<float>   *Reco_DiTrack_d0Sig;
   vector<float>   *Gen_DiTrack_pdgId;
   vector<float>   *Gen_DiTrack_pt;
   vector<float>   *Gen_DiTrack_eta;
   vector<float>   *Gen_DiTrack_phi;
   vector<float>   *Gen_DiTrack_charge;
   vector<float>   *Gen_DiTrack_mass;
   vector<int>     *Gen_DiTrack_isSwap;
   UInt_t          Reco_TriTrack_N;
   vector<int>     *Reco_TriTrack_Pion_idx;
   vector<int>     *Reco_TriTrack_Kaon_idx;
   vector<int>     *Reco_TriTrack_Proton_idx;
   vector<float>   *Reco_TriTrack_pt;
   vector<float>   *Reco_TriTrack_rap;
   vector<float>   *Reco_TriTrack_mass;
   vector<float>   *Reco_TriTrack_vProb;
   vector<float>   *Reco_TriTrack_alpha;
   vector<float>   *Reco_TriTrack_d0Sig;

   // List of branches
   TBranch        *b_Event_Run;   //!
   TBranch        *b_Event_Lumi;   //!
   TBranch        *b_Event_Bx;   //!
   TBranch        *b_Event_Orbit;   //!
   TBranch        *b_Event_Number;   //!
   TBranch        *b_Event_nPV;   //!
   TBranch        *b_Event_PriVtx_Pos;   //!
   TBranch        *b_Event_PriVtx_Err;   //!
   TBranch        *b_Event_hiHF;   //!
   TBranch        *b_Event_hiBin;   //!
   TBranch        *b_Reco_Track_N;   //!
   TBranch        *b_Reco_Track_beta_PI;   //!
   TBranch        *b_Reco_Track_t0_PI;   //!
   TBranch        *b_Reco_Track_t0Err_PI;   //!
   TBranch        *b_Reco_Track_beta_PV;   //!
   TBranch        *b_Reco_Track_t0_PV;   //!
   TBranch        *b_Reco_Track_t0Err_PV;   //!
   TBranch        *b_Reco_Track_massSq_PV;   //!
   TBranch        *b_Reco_Track_beta_PID;   //!
   TBranch        *b_Reco_Track_t0_PID;   //!
   TBranch        *b_Reco_Track_t0Err_PID;   //!
   TBranch        *b_Reco_Track_massSq_PID;   //!
   TBranch        *b_Reco_Track_pathLength;   //!
   TBranch        *b_Reco_Track_p;   //!
   TBranch        *b_Reco_Track_tMTD;   //!
   TBranch        *b_Reco_Track_tMTDErr;   //!
   TBranch        *b_Reco_Track_eta;   //!
   TBranch        *b_Reco_Track_phi;   //!
   TBranch        *b_Reco_Track_charge;   //!
   TBranch        *b_Reco_Track_quality;   //!
   TBranch        *b_Reco_Track_pt;   //!
   TBranch        *b_Reco_Track_ptErr;   //!
   TBranch        *b_Reco_Track_dcaDz;   //!
   TBranch        *b_Reco_Track_dcaDxy;   //!
   TBranch        *b_Reco_Track_normChi2;   //!
   TBranch        *b_Reco_Track_nHits;   //!
   TBranch        *b_Reco_Track_dEdx;   //!
   TBranch        *b_Reco_Track_dEdxErr;   //!
   TBranch        *b_Reco_Track_nSatMea;   //!
   TBranch        *b_Reco_Track_nMea;   //!
   TBranch        *b_Gen_Track_pdgId;   //!
   TBranch        *b_Gen_Track_pt;   //!
   TBranch        *b_Gen_Track_eta;   //!
   TBranch        *b_Gen_Track_phi;   //!
   TBranch        *b_Gen_Track_t0;   //!
   TBranch        *b_Gen_Track_charge;   //!
   TBranch        *b_Gen_Track_mass;   //!
   TBranch        *b_Gen_Track_momPdgId;   //!
   TBranch        *b_Reco_DiTrack_N;   //!
   TBranch        *b_Reco_DiTrack_Pion_idx;   //!
   TBranch        *b_Reco_DiTrack_Kaon_idx;   //!
   TBranch        *b_Reco_DiTrack_pt;   //!
   TBranch        *b_Reco_DiTrack_rap;   //!
   TBranch        *b_Reco_DiTrack_mass;   //!
   TBranch        *b_Reco_DiTrack_vProb;   //!
   TBranch        *b_Reco_DiTrack_alpha;   //!
   TBranch        *b_Reco_DiTrack_d0Sig;   //!
   TBranch        *b_Gen_DiTrack_pdgId;   //!
   TBranch        *b_Gen_DiTrack_pt;   //!
   TBranch        *b_Gen_DiTrack_eta;   //!
   TBranch        *b_Gen_DiTrack_phi;   //!
   TBranch        *b_Gen_DiTrack_charge;   //!
   TBranch        *b_Gen_DiTrack_mass;   //!
   TBranch        *b_Gen_DiTrack_isSwap;   //!
   TBranch        *b_Reco_TriTrack_N;   //!
   TBranch        *b_Reco_TriTrack_Pion_idx;   //!
   TBranch        *b_Reco_TriTrack_Kaon_idx;   //!
   TBranch        *b_Reco_TriTrack_Proton_idx;   //!
   TBranch        *b_Reco_TriTrack_pt;   //!
   TBranch        *b_Reco_TriTrack_rap;   //!
   TBranch        *b_Reco_TriTrack_mass;   //!
   TBranch        *b_Reco_TriTrack_vProb;   //!
   TBranch        *b_Reco_TriTrack_alpha;   //!
   TBranch        *b_Reco_TriTrack_d0Sig;   //!

   HiMTDTree(TTree *tree=0);
   virtual ~HiMTDTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HiMTDTree_cxx
HiMTDTree::HiMTDTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/storage1/users/wl33/MTD/HiMTDTree_numEvent1000.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/storage1/users/wl33/MTD/HiMTDTree_numEvent1000.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/storage1/users/wl33/MTD/HiMTDTree_numEvent1000.root:/timeAna");
      dir->GetObject("HiMTDTree",tree);

   }
   Init(tree);
}

HiMTDTree::~HiMTDTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HiMTDTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HiMTDTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HiMTDTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Event_PriVtx_Pos = 0;
   Event_PriVtx_Err = 0;
   Reco_Track_beta_PI = 0;
   Reco_Track_t0_PI = 0;
   Reco_Track_t0Err_PI = 0;
   Reco_Track_beta_PV = 0;
   Reco_Track_t0_PV = 0;
   Reco_Track_t0Err_PV = 0;
   Reco_Track_massSq_PV = 0;
   Reco_Track_beta_PID = 0;
   Reco_Track_t0_PID = 0;
   Reco_Track_t0Err_PID = 0;
   Reco_Track_massSq_PID = 0;
   Reco_Track_pathLength = 0;
   Reco_Track_p = 0;
   Reco_Track_tMTD = 0;
   Reco_Track_tMTDErr = 0;
   Reco_Track_eta = 0;
   Reco_Track_phi = 0;
   Reco_Track_charge = 0;
   Reco_Track_quality = 0;
   Reco_Track_pt = 0;
   Reco_Track_ptErr = 0;
   Reco_Track_dcaDz = 0;
   Reco_Track_dcaDxy = 0;
   Reco_Track_normChi2 = 0;
   Reco_Track_nHits = 0;
   Reco_Track_dEdx = 0;
   Reco_Track_dEdxErr = 0;
   Reco_Track_nSatMea = 0;
   Reco_Track_nMea = 0;
   Gen_Track_pdgId = 0;
   Gen_Track_pt = 0;
   Gen_Track_eta = 0;
   Gen_Track_phi = 0;
   Gen_Track_t0 = 0;
   Gen_Track_charge = 0;
   Gen_Track_mass = 0;
   Gen_Track_momPdgId = 0;
   Reco_DiTrack_Pion_idx = 0;
   Reco_DiTrack_Kaon_idx = 0;
   Reco_DiTrack_pt = 0;
   Reco_DiTrack_rap = 0;
   Reco_DiTrack_mass = 0;
   Reco_DiTrack_vProb = 0;
   Reco_DiTrack_alpha = 0;
   Reco_DiTrack_d0Sig = 0;
   Gen_DiTrack_pdgId = 0;
   Gen_DiTrack_pt = 0;
   Gen_DiTrack_eta = 0;
   Gen_DiTrack_phi = 0;
   Gen_DiTrack_charge = 0;
   Gen_DiTrack_mass = 0;
   Gen_DiTrack_isSwap = 0;
   Reco_TriTrack_Pion_idx = 0;
   Reco_TriTrack_Kaon_idx = 0;
   Reco_TriTrack_Proton_idx = 0;
   Reco_TriTrack_pt = 0;
   Reco_TriTrack_rap = 0;
   Reco_TriTrack_mass = 0;
   Reco_TriTrack_vProb = 0;
   Reco_TriTrack_alpha = 0;
   Reco_TriTrack_d0Sig = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event_Run", &Event_Run, &b_Event_Run);
   fChain->SetBranchAddress("Event_Lumi", &Event_Lumi, &b_Event_Lumi);
   fChain->SetBranchAddress("Event_Bx", &Event_Bx, &b_Event_Bx);
   fChain->SetBranchAddress("Event_Orbit", &Event_Orbit, &b_Event_Orbit);
   fChain->SetBranchAddress("Event_Number", &Event_Number, &b_Event_Number);
   fChain->SetBranchAddress("Event_nPV", &Event_nPV, &b_Event_nPV);
   fChain->SetBranchAddress("Event_PriVtx_Pos", &Event_PriVtx_Pos, &b_Event_PriVtx_Pos);
   fChain->SetBranchAddress("Event_PriVtx_Err", &Event_PriVtx_Err, &b_Event_PriVtx_Err);
   fChain->SetBranchAddress("Event_hiHF", &Event_hiHF, &b_Event_hiHF);
   fChain->SetBranchAddress("Event_hiBin", &Event_hiBin, &b_Event_hiBin);
   fChain->SetBranchAddress("Reco_Track_N", &Reco_Track_N, &b_Reco_Track_N);
   fChain->SetBranchAddress("Reco_Track_beta_PI", &Reco_Track_beta_PI, &b_Reco_Track_beta_PI);
   fChain->SetBranchAddress("Reco_Track_t0_PI", &Reco_Track_t0_PI, &b_Reco_Track_t0_PI);
   fChain->SetBranchAddress("Reco_Track_t0Err_PI", &Reco_Track_t0Err_PI, &b_Reco_Track_t0Err_PI);
   fChain->SetBranchAddress("Reco_Track_beta_PV", &Reco_Track_beta_PV, &b_Reco_Track_beta_PV);
   fChain->SetBranchAddress("Reco_Track_t0_PV", &Reco_Track_t0_PV, &b_Reco_Track_t0_PV);
   fChain->SetBranchAddress("Reco_Track_t0Err_PV", &Reco_Track_t0Err_PV, &b_Reco_Track_t0Err_PV);
   fChain->SetBranchAddress("Reco_Track_massSq_PV", &Reco_Track_massSq_PV, &b_Reco_Track_massSq_PV);
   fChain->SetBranchAddress("Reco_Track_beta_PID", &Reco_Track_beta_PID, &b_Reco_Track_beta_PID);
   fChain->SetBranchAddress("Reco_Track_t0_PID", &Reco_Track_t0_PID, &b_Reco_Track_t0_PID);
   fChain->SetBranchAddress("Reco_Track_t0Err_PID", &Reco_Track_t0Err_PID, &b_Reco_Track_t0Err_PID);
   fChain->SetBranchAddress("Reco_Track_massSq_PID", &Reco_Track_massSq_PID, &b_Reco_Track_massSq_PID);
   fChain->SetBranchAddress("Reco_Track_pathLength", &Reco_Track_pathLength, &b_Reco_Track_pathLength);
   fChain->SetBranchAddress("Reco_Track_p", &Reco_Track_p, &b_Reco_Track_p);
   fChain->SetBranchAddress("Reco_Track_tMTD", &Reco_Track_tMTD, &b_Reco_Track_tMTD);
   fChain->SetBranchAddress("Reco_Track_tMTDErr", &Reco_Track_tMTDErr, &b_Reco_Track_tMTDErr);
   fChain->SetBranchAddress("Reco_Track_eta", &Reco_Track_eta, &b_Reco_Track_eta);
   fChain->SetBranchAddress("Reco_Track_phi", &Reco_Track_phi, &b_Reco_Track_phi);
   fChain->SetBranchAddress("Reco_Track_charge", &Reco_Track_charge, &b_Reco_Track_charge);
   fChain->SetBranchAddress("Reco_Track_quality", &Reco_Track_quality, &b_Reco_Track_quality);
   fChain->SetBranchAddress("Reco_Track_pt", &Reco_Track_pt, &b_Reco_Track_pt);
   fChain->SetBranchAddress("Reco_Track_ptErr", &Reco_Track_ptErr, &b_Reco_Track_ptErr);
   fChain->SetBranchAddress("Reco_Track_dcaDz", &Reco_Track_dcaDz, &b_Reco_Track_dcaDz);
   fChain->SetBranchAddress("Reco_Track_dcaDxy", &Reco_Track_dcaDxy, &b_Reco_Track_dcaDxy);
   fChain->SetBranchAddress("Reco_Track_normChi2", &Reco_Track_normChi2, &b_Reco_Track_normChi2);
   fChain->SetBranchAddress("Reco_Track_nHits", &Reco_Track_nHits, &b_Reco_Track_nHits);
   fChain->SetBranchAddress("Reco_Track_dEdx", &Reco_Track_dEdx, &b_Reco_Track_dEdx);
   fChain->SetBranchAddress("Reco_Track_dEdxErr", &Reco_Track_dEdxErr, &b_Reco_Track_dEdxErr);
   fChain->SetBranchAddress("Reco_Track_nSatMea", &Reco_Track_nSatMea, &b_Reco_Track_nSatMea);
   fChain->SetBranchAddress("Reco_Track_nMea", &Reco_Track_nMea, &b_Reco_Track_nMea);
   fChain->SetBranchAddress("Gen_Track_pdgId", &Gen_Track_pdgId, &b_Gen_Track_pdgId);
   fChain->SetBranchAddress("Gen_Track_pt", &Gen_Track_pt, &b_Gen_Track_pt);
   fChain->SetBranchAddress("Gen_Track_eta", &Gen_Track_eta, &b_Gen_Track_eta);
   fChain->SetBranchAddress("Gen_Track_phi", &Gen_Track_phi, &b_Gen_Track_phi);
   fChain->SetBranchAddress("Gen_Track_t0", &Gen_Track_t0, &b_Gen_Track_t0);
   fChain->SetBranchAddress("Gen_Track_charge", &Gen_Track_charge, &b_Gen_Track_charge);
   fChain->SetBranchAddress("Gen_Track_mass", &Gen_Track_mass, &b_Gen_Track_mass);
   fChain->SetBranchAddress("Gen_Track_momPdgId", &Gen_Track_momPdgId, &b_Gen_Track_momPdgId);
   fChain->SetBranchAddress("Reco_DiTrack_N", &Reco_DiTrack_N, &b_Reco_DiTrack_N);
   fChain->SetBranchAddress("Reco_DiTrack_Pion_idx", &Reco_DiTrack_Pion_idx, &b_Reco_DiTrack_Pion_idx);
   fChain->SetBranchAddress("Reco_DiTrack_Kaon_idx", &Reco_DiTrack_Kaon_idx, &b_Reco_DiTrack_Kaon_idx);
   fChain->SetBranchAddress("Reco_DiTrack_pt", &Reco_DiTrack_pt, &b_Reco_DiTrack_pt);
   fChain->SetBranchAddress("Reco_DiTrack_rap", &Reco_DiTrack_rap, &b_Reco_DiTrack_rap);
   fChain->SetBranchAddress("Reco_DiTrack_mass", &Reco_DiTrack_mass, &b_Reco_DiTrack_mass);
   fChain->SetBranchAddress("Reco_DiTrack_vProb", &Reco_DiTrack_vProb, &b_Reco_DiTrack_vProb);
   fChain->SetBranchAddress("Reco_DiTrack_alpha", &Reco_DiTrack_alpha, &b_Reco_DiTrack_alpha);
   fChain->SetBranchAddress("Reco_DiTrack_d0Sig", &Reco_DiTrack_d0Sig, &b_Reco_DiTrack_d0Sig);
   fChain->SetBranchAddress("Gen_DiTrack_pdgId", &Gen_DiTrack_pdgId, &b_Gen_DiTrack_pdgId);
   fChain->SetBranchAddress("Gen_DiTrack_pt", &Gen_DiTrack_pt, &b_Gen_DiTrack_pt);
   fChain->SetBranchAddress("Gen_DiTrack_eta", &Gen_DiTrack_eta, &b_Gen_DiTrack_eta);
   fChain->SetBranchAddress("Gen_DiTrack_phi", &Gen_DiTrack_phi, &b_Gen_DiTrack_phi);
   fChain->SetBranchAddress("Gen_DiTrack_charge", &Gen_DiTrack_charge, &b_Gen_DiTrack_charge);
   fChain->SetBranchAddress("Gen_DiTrack_mass", &Gen_DiTrack_mass, &b_Gen_DiTrack_mass);
   fChain->SetBranchAddress("Gen_DiTrack_isSwap", &Gen_DiTrack_isSwap, &b_Gen_DiTrack_isSwap);
   fChain->SetBranchAddress("Reco_TriTrack_N", &Reco_TriTrack_N, &b_Reco_TriTrack_N);
   fChain->SetBranchAddress("Reco_TriTrack_Pion_idx", &Reco_TriTrack_Pion_idx, &b_Reco_TriTrack_Pion_idx);
   fChain->SetBranchAddress("Reco_TriTrack_Kaon_idx", &Reco_TriTrack_Kaon_idx, &b_Reco_TriTrack_Kaon_idx);
   fChain->SetBranchAddress("Reco_TriTrack_Proton_idx", &Reco_TriTrack_Proton_idx, &b_Reco_TriTrack_Proton_idx);
   fChain->SetBranchAddress("Reco_TriTrack_pt", &Reco_TriTrack_pt, &b_Reco_TriTrack_pt);
   fChain->SetBranchAddress("Reco_TriTrack_rap", &Reco_TriTrack_rap, &b_Reco_TriTrack_rap);
   fChain->SetBranchAddress("Reco_TriTrack_mass", &Reco_TriTrack_mass, &b_Reco_TriTrack_mass);
   fChain->SetBranchAddress("Reco_TriTrack_vProb", &Reco_TriTrack_vProb, &b_Reco_TriTrack_vProb);
   fChain->SetBranchAddress("Reco_TriTrack_alpha", &Reco_TriTrack_alpha, &b_Reco_TriTrack_alpha);
   fChain->SetBranchAddress("Reco_TriTrack_d0Sig", &Reco_TriTrack_d0Sig, &b_Reco_TriTrack_d0Sig);
   Notify();
}

Bool_t HiMTDTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HiMTDTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HiMTDTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HiMTDTree_cxx
