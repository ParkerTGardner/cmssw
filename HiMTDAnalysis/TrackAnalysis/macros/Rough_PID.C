#include "HiMTDTree.C"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraph.h"
//#include <iostream>
#include <vector>
#include "math.h"

double yval(double M, double ETA, double PT) {
    double result;
    result=log(
        (pow(M*M+PT*PT*cosh(ETA)*cosh(ETA),0.5)
        + PT*sinh(ETA))
        / (pow(M*M+PT*PT,0.5))
        );
    return result; 
}


void Rough_PID()
{

TH1::SetDefaultSumw2(kTRUE);
TProfile::SetDefaultSumw2(kTRUE); 


TFile* f1 = new TFile("/storage1/users/wl33/MTD/HiMTDTree_numEvent10000.root");
//Replace with pwd to Tree we want to analyze.

TTree* t= (TTree*) f1->Get("timeAna/HiMTDTree");
HiMTDTree tree(t);

//long int nevents = tree.fChain->GetEntries();
int nevents=30;

double betaE1=1/1.2;
double betaE2=1/0.5;
double mP=.938272;
double mPI=.129570;
double mK=.493677;
float multrms=1;
int betaBin=300;
int pBin=100;
double HLmt=0.3;
float err=1.0;

TH2D* ibetaVp = new TH2D("InvBetaVP", "InvBeta vs. P", pBin, 0, 5, betaBin, betaE1, betaE2);

TH2D* ibetaMTP = new TH2D("InvBetaMTP", "1/B Meas minus Theory  vs. P (Proton)", pBin, 0, 5, betaBin, -HLmt, HLmt);
TH2D* ibetaMTPI = new TH2D("InvBetaMTPI", "1/B Meas minus Theory  vs. P (Pion)", pBin, 0, 5, betaBin, -HLmt, HLmt);
TH2D* ibetaMTK = new TH2D("InvBetaMTPK", "1/B Meas minus Theory  vs. P (Kaon)", pBin, 0, 5, betaBin, -HLmt, HLmt);

for(long int ievent=0; ievent<nevents; ievent++){
    tree.GetEntry(ievent);
    for(int itrk=0; itrk<tree.Reco_Track_p->size(); itrk++){
        bool passCut = false;
        if((*tree.Reco_Track_beta_PV)[itrk] > 0.){ passCut = true;}
        if(passCut){
            ibetaVp->Fill((*tree.Reco_Track_p)[itrk], 1./(*tree.Reco_Track_beta_PV)[itrk]);
            if(fabs((*tree.Gen_Track_pdgId)[itrk]) < 2212 + err  && fabs((*tree.Gen_Track_pdgId)[itrk]) > 2211 - err ){
                ibetaMTP->Fill((*tree.Reco_Track_p)[itrk],  (  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mP,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5))  );
                }
            if(fabs((*tree.Gen_Track_pdgId)[itrk]) < 211 + err && fabs((*tree.Gen_Track_pdgId)[itrk]) > 211 - err ){
                ibetaMTPI->Fill((*tree.Reco_Track_p)[itrk], (  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mPI,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5)) );
                }
            if(fabs((*tree.Gen_Track_pdgId)[itrk]) < 321 + err  && fabs((*tree.Gen_Track_pdgId)[itrk]) > 321 - err ){
                ibetaMTK->Fill((*tree.Reco_Track_p)[itrk],  (  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mK,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5))  );
                }
            //cout << (*tree.Gen_Track_pdgId)[1] << endl;
            }
        }   
    }


TF1 *faP = new TF1("faP","(1+([0]^2)/(x^2))^(0.5)",0,5);
TF1 *faPI = new TF1("faPI","(1+([0]^2)/(x^2))^(0.5)",0,5);
TF1 *faK = new TF1("faK","(1+([0]^2)/(x^2))^(0.5)",0,5);

faP->SetParameter(0,mP);
faPI->SetParameter(0,mPI);
faK->SetParameter(0,mK);

TCanvas* c1 = new TCanvas("c1", "", 900, 800);
c1->cd();

ibetaVp->DrawCopy("COLZ");
faP->Draw("SAME");
faPI->Draw("SAME");
faK->Draw("SAME");

TCanvas* c2 = new TCanvas("c2", "", 1200, 800);
c2->Divide(3);

TLine *lLOW = new TLine(0,-.1,5,-.1);
TLine *lHI = new TLine(0,.1,5,.1);

c2->cd(1);
ibetaMTP->DrawCopy("COLZ");

c2->cd(2);
ibetaMTPI->DrawCopy("COLZ");

c2->cd(3);
ibetaMTK->DrawCopy("COLZ");

cout << ibetaVp->GetEntries() << endl;
TFile o1file("betaVp.root", "recreate");
ibetaVp->Write();
cout << ibetaMTP->GetEntries() << endl;
TFile o2file("MTP.root", "recreate");
ibetaMTP->Write();
cout << ibetaMTPI->GetEntries() << endl;
TFile o3file("MTPI.root", "recreate");
ibetaMTPI->Write();
cout << ibetaMTK->GetEntries() << endl;
TFile o4file("MTK.root", "recreate");
ibetaMTK->Write();


int profBin=10;
int totBin=floor(pBin/profBin);

TCanvas* c6 = new TCanvas("c6", "", 1200, 800);
c6->Divide(3);
c6->cd(1);
auto pxMTK = ibetaMTK->ProfileX("PXK", 1, -1, "S");
pxMTK->Rebin(profBin);
pxMTK->DrawCopy();
c6->cd(2);
auto pxMTP   = ibetaMTP->ProfileX("PXP", 1, -1, "S");
pxMTP->Rebin(profBin);
pxMTP->DrawCopy();
c6->cd(3);
auto pxMTPI = ibetaMTPI->ProfileX("PXPI", 1, -1, "S");
pxMTPI->Rebin(profBin);
pxMTPI->DrawCopy();

std::vector<double> Prms(totBin);
std::vector<double> Krms(totBin);
std::vector<double> PIrms(totBin);

for(int i=1; i<totBin+1; i++){
Prms[i] = pxMTP->GetBinError(i);
Krms[i] = pxMTK->GetBinError(i);
PIrms[i] = pxMTPI->GetBinError(i);
}


TH2D* allMTP = new TH2D("allMTP", "Before RMS cut: 1/B Meas minus Theory  vs. P (Proton)", pBin, 0, 5, betaBin, -HLmt, HLmt);
TH2D* someMTP = new TH2D("someMTP", "After RMS cut: 1/B Meas minus Theory  vs. P (Proton)", pBin, 0, 5, betaBin, -HLmt, HLmt);
std::vector< double > ProtonY;
std::vector< double > ProtonPT;


for(long int ievent=0; ievent<nevents; ievent++){
    tree.GetEntry(ievent);
    for(int itrk=0; itrk<tree.Reco_Track_p->size(); itrk++){
        bool passCut = false;
        if(
            (*tree.Reco_Track_beta_PV)[itrk] > 0.
            //&& fabs((*tree.Gen_Track_pdgId)[itrk]) < 2212 + err
            //&& fabs((*tree.Gen_Track_pdgId)[itrk]) > 2211 - err            
            ){passCut=true;}
        if(passCut){
            allMTP->Fill((*tree.Reco_Track_p)[itrk],  (  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mP,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5))  );
            for(int i=1; i<totBin+1; i++){
                auto lim1=pxMTP->GetBinLowEdge(i);
                auto lim2=pxMTP->GetBinLowEdge(i+1);
                if(
                    (*tree.Reco_Track_p)[itrk] < lim2 &&
                    (*tree.Reco_Track_p)[itrk] >= lim1
                    && fabs(  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mP,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5)) <= multrms*Prms[i]
                    && fabs(  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mK,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5)) > multrms*Krms[i]
                    && fabs(  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mPI,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5)) > multrms*PIrms[i]
                    ){
                    someMTP->Fill((*tree.Reco_Track_p)[itrk],  (  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mP,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5))  );
                    //preparing TGraph for Y vs pT
                    ProtonY.push_back(abs(yval(mP,(*tree.Reco_Track_eta)[itrk],(*tree.Reco_Track_pt)[itrk])));
                    ProtonPT.push_back((*tree.Reco_Track_pt)[itrk]); 
                     }
                }
            }
        }
    }

TCanvas* c7 = new TCanvas("c7", "", 1200, 800);
c7->Divide(2);
c7->cd(1);
allMTP->DrawCopy("COLZ");
c7->cd(2);
someMTP->DrawCopy("COLZ");



TH2D* allMTK = new TH2D("allMTK", "Before RMS cut: 1/B Meas minus Theory  vs. P (Kaon)", pBin, 0, 5, betaBin, -HLmt, HLmt);
TH2D* someMTK = new TH2D("someMTK", "After RMS cut: 1/B Meas minus Theory  vs. P (Kaon)", pBin, 0, 5, betaBin, -HLmt, HLmt);


for(long int ievent=0; ievent<nevents; ievent++){
    tree.GetEntry(ievent);
    for(int itrk=0; itrk<tree.Reco_Track_p->size(); itrk++){
        bool passCut = false;
        if(
            (*tree.Reco_Track_beta_PV)[itrk] > 0.
            //&& fabs((*tree.Gen_Track_pdgId)[itrk]) < 321 + err
            //&& fabs((*tree.Gen_Track_pdgId)[itrk]) > 321 - err            
            ){passCut=true;}
        if(passCut){
            allMTK->Fill((*tree.Reco_Track_p)[itrk],  (  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mK,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5))  );
            for(int i=1; i<totBin+1; i++){
                auto lim1=pxMTK->GetBinLowEdge(i);
                auto lim2=pxMTK->GetBinLowEdge(i+1);
                if(
                    (*tree.Reco_Track_p)[itrk] < lim2 &&
                    (*tree.Reco_Track_p)[itrk] >= lim1
                    && fabs(  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mP,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5)) > multrms*Prms[i]
                    && fabs(  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mK,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5)) <= multrms*Krms[i]
                    && fabs(  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mPI,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5)) > multrms*PIrms[i]
                    ){
                     someMTK->Fill((*tree.Reco_Track_p)[itrk],  (  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mK,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5))  );
                     }
                }
            }
        }
    }

TCanvas* c8 = new TCanvas("c8", "", 1200, 800);
c8->Divide(2);
c8->cd(1);
allMTK->DrawCopy("COLZ");
c8->cd(2);
someMTK->DrawCopy("COLZ");




TH2D* allMTPI = new TH2D("allMTPI", "Before RMS cut: 1/B Meas minus Theory  vs. P (Pion)", pBin, 0, 5, betaBin, -HLmt, HLmt);
TH2D* someMTPI = new TH2D("someMTPI", "After RMS cut:  1/B Meas minus Theory  vs. P (Pion)", pBin, 0, 5, betaBin, -HLmt, HLmt);


for(long int ievent=0; ievent<nevents; ievent++){
    tree.GetEntry(ievent);
    for(int itrk=0; itrk<tree.Reco_Track_p->size(); itrk++){
        bool passCut = false;
        if(
            (*tree.Reco_Track_beta_PV)[itrk] > 0.
            //&& fabs((*tree.Gen_Track_pdgId)[itrk]) < 211 + err 
            //&& fabs((*tree.Gen_Track_pdgId)[itrk]) > 211 - err            
            ){passCut=true;}
        if(passCut){
            allMTPI->Fill((*tree.Reco_Track_p)[itrk],  (  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mPI,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5))  );
            for(int i=1; i<totBin+1; i++){
                auto lim1=pxMTPI->GetBinLowEdge(i);
                auto lim2=pxMTPI->GetBinLowEdge(i+1);
                if(
                    (*tree.Reco_Track_p)[itrk] < lim2 &&
                    (*tree.Reco_Track_p)[itrk] >= lim1
                    && fabs(  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mP,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5)) > multrms*Prms[i]
                    && fabs(  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mK,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5)) > multrms*Krms[i]
                    && fabs(  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mPI,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5)) <= multrms*PIrms[i]
                    ){
                     someMTPI->Fill((*tree.Reco_Track_p)[itrk],  (  (1./(*tree.Reco_Track_beta_PV)[itrk])-pow((1+(pow(mPI,2))/pow((*tree.Reco_Track_p)[itrk],2)),0.5))  );
                     }
                }
            }
        }
    }

TCanvas* c9 = new TCanvas("c9", "", 1200, 800);
c9->Divide(2);
c9->cd(1);
allMTPI->DrawCopy("COLZ");
c9->cd(2);
someMTPI->DrawCopy("COLZ");
/*
TGraph* ProtonPTY = new TGraph(ProtonPT.size(),&ProtonY[0],&ProtonPT[0]);
TCanvas* c10 = new TCanvas("c10", "", 1200, 800);
ProtonPTY->Draw("A*");
*/



}
