#include <iostream>
#include<fstream>
#include<string>
#include "TObject.h"
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include <ctime>


using namespace std;

double PhiInRange(const double& phi) {
  double phiout = phi;

  if( phiout > 2*M_PI || phiout < -2*M_PI) {
    phiout = fmod( phiout, 2*M_PI);
  }
  if (phiout <= -M_PI) phiout += 2*M_PI;
  else if (phiout >  M_PI) phiout -= 2*M_PI;

  return phiout;
}

double delta2R(double eta1, double phi1, double eta2, double phi2) {
  return sqrt(pow(eta1 - eta2,2) +pow(PhiInRange(phi1 - phi2),2));
}


float leppt, lepy, b1pt, b2pt, b1y, b2y, bb_inv_mass, pho1pt, pho2pt, pho1y, pho2y, dipho_invmass, invmassbbgg, invmassbbgg_jesup, invmassbbgg_jesdn, invmassbbgg_resoup, invmassbbgg_resodn, pho1eta, pho1phi, pho1e, pho2eta, pho2phi, pho2e;
float b1_DeepFlv, b2_DeepFlv, pho1MVA, pho2MVA, gg_dR, MetPhi, delphi_ggEt;
double weight_nom, weight_PU_up, weight_leptonsf_up, weight_bSF_up, weight_trig_up, weight_PU_dn, weight_leptonsf_dn, weight_bSF_dn, weight_trig_dn;
int bb, ab;
bool isEle, isMu;
float sigmabb1, sigmah1;
float chi2;
double timetaken;


void sig_find_delphi_v1()
{

	clock_t start;
	clock_t end;

	start = clock();

	int n;
	cout << "Enter bin no: " << endl;
	cin >> n;


	double lumi_18 = 59730;


	TFile *f1 = new TFile("TTGJetsV2.root", "read");
        TFile *f2 = new TFile("DYJetsToLLM50V2.root", "read");
        TFile *f3 = new TFile("TTSLV2.root", "read");
        TFile *f4 = new TFile("TTTo2L2NuV2.root", "read");
	TFile *f5 = new TFile("WH_mA20_v7.root", "read");
        TFile *f6 = new TFile("WH_mA55_v7.root", "read");


	TTree *Tout1 = (TTree*)f1->Get("Tout");
        TTree *Tout2 = (TTree*)f2->Get("Tout");
        TTree *Tout3 = (TTree*)f3->Get("Tout");
        TTree *Tout4 = (TTree*)f4->Get("Tout");
	TTree *Tout5 = (TTree*)f5->Get("Tout");
	TTree *Tout6 = (TTree*)f6->Get("Tout");


	char name[1000]; 
	TFile *fout;
        sprintf(name,"sig_finder_v1/combined_2018_wh_v1_bin%i.root",n);
        fout = new TFile(name,"RECREATE");
  

    TH1F* tt_gjets = new TH1F("tt_gjets", "tt_gjets", n, 0.0, M_PI);
    tt_gjets->Sumw2();

    TH1F* dy = new TH1F("dy", "dy", n, 0.0, M_PI);
    dy->Sumw2();

    TH1F* tt_sl = new TH1F("tt_sl", "tt_sl", n, 0.0, M_PI);
    tt_sl->Sumw2();

    TH1F* tt_2l2nu = new TH1F("tt_2l2nu", "tt_2l2nu", n, 0.0, M_PI);
    tt_2l2nu->Sumw2();

    TH1F* data_obs = new TH1F("data_obs", "data_obs", n, 0.0, M_PI);
    data_obs->Sumw2();

    TH1F* sig_20 = new TH1F("sig_20", "sig_20", n, 0.0, M_PI);
    sig_20->Sumw2();

    TH1F* sig_55 = new TH1F("sig_55", "sig_55", n, 0.0, M_PI);
    sig_55->Sumw2();

    

    Tout1->SetBranchAddress("bb_inv_mass",&bb_inv_mass);
    Tout1->SetBranchAddress("dipho_invmass",&dipho_invmass);
    Tout1->SetBranchAddress("invmassbbgg",&invmassbbgg);
    Tout1->SetBranchAddress("weight_nom",&weight_nom);
    Tout1->SetBranchAddress("ab",&ab);
    Tout1->SetBranchAddress("bb",&bb);
    Tout1->SetBranchAddress("isEle",&isEle);
    Tout1->SetBranchAddress("isMu",&isMu);
    Tout1->SetBranchAddress("pho1pt",&pho1pt);
    Tout1->SetBranchAddress("pho2pt",&pho2pt);
    Tout1->SetBranchAddress("pho1y",&pho1y);
    Tout1->SetBranchAddress("pho2y",&pho2y);
    Tout1->SetBranchAddress("pho1eta",&pho1eta);
    Tout1->SetBranchAddress("pho1phi",&pho1phi);
    Tout1->SetBranchAddress("pho1e",&pho1e);
    Tout1->SetBranchAddress("pho2eta",&pho2eta);
    Tout1->SetBranchAddress("pho2phi",&pho2phi);
    Tout1->SetBranchAddress("pho2e",&pho2e);
    Tout1->SetBranchAddress("MetPhi",&MetPhi);

    Tout2->SetBranchAddress("bb_inv_mass",&bb_inv_mass);
    Tout2->SetBranchAddress("dipho_invmass",&dipho_invmass);
    Tout2->SetBranchAddress("invmassbbgg",&invmassbbgg);
    Tout2->SetBranchAddress("weight_nom",&weight_nom);
    Tout2->SetBranchAddress("ab",&ab);
    Tout2->SetBranchAddress("bb",&bb);
    Tout2->SetBranchAddress("isEle",&isEle);
    Tout2->SetBranchAddress("isMu",&isMu);
    Tout2->SetBranchAddress("pho1pt",&pho1pt);
    Tout2->SetBranchAddress("pho2pt",&pho2pt);
    Tout2->SetBranchAddress("pho1y",&pho1y);
    Tout2->SetBranchAddress("pho2y",&pho2y);
    Tout2->SetBranchAddress("pho1eta",&pho1eta);
    Tout2->SetBranchAddress("pho1phi",&pho1phi);
    Tout2->SetBranchAddress("pho1e",&pho1e);
    Tout2->SetBranchAddress("pho2eta",&pho2eta);
    Tout2->SetBranchAddress("pho2phi",&pho2phi);
    Tout2->SetBranchAddress("pho2e",&pho2e);
    Tout2->SetBranchAddress("MetPhi",&MetPhi);

    Tout3->SetBranchAddress("bb_inv_mass",&bb_inv_mass);
    Tout3->SetBranchAddress("dipho_invmass",&dipho_invmass);
    Tout3->SetBranchAddress("invmassbbgg",&invmassbbgg);
    Tout3->SetBranchAddress("weight_nom",&weight_nom);
    Tout3->SetBranchAddress("ab",&ab);
    Tout3->SetBranchAddress("bb",&bb);
    Tout3->SetBranchAddress("isEle",&isEle);
    Tout3->SetBranchAddress("isMu",&isMu);
    Tout3->SetBranchAddress("pho1pt",&pho1pt);
    Tout3->SetBranchAddress("pho2pt",&pho2pt);
    Tout3->SetBranchAddress("pho1y",&pho1y);
    Tout3->SetBranchAddress("pho2y",&pho2y);
    Tout3->SetBranchAddress("pho1eta",&pho1eta);
    Tout3->SetBranchAddress("pho1phi",&pho1phi);
    Tout3->SetBranchAddress("pho1e",&pho1e);
    Tout3->SetBranchAddress("pho2eta",&pho2eta);
    Tout3->SetBranchAddress("pho2phi",&pho2phi);
    Tout3->SetBranchAddress("pho2e",&pho2e);
    Tout3->SetBranchAddress("MetPhi",&MetPhi);

    Tout4->SetBranchAddress("bb_inv_mass",&bb_inv_mass);
    Tout4->SetBranchAddress("dipho_invmass",&dipho_invmass);
    Tout4->SetBranchAddress("invmassbbgg",&invmassbbgg);
    Tout4->SetBranchAddress("weight_nom",&weight_nom);
    Tout4->SetBranchAddress("ab",&ab);
    Tout4->SetBranchAddress("bb",&bb);
    Tout4->SetBranchAddress("isEle",&isEle);
    Tout4->SetBranchAddress("isMu",&isMu);
    Tout4->SetBranchAddress("pho1pt",&pho1pt);
    Tout4->SetBranchAddress("pho2pt",&pho2pt);
    Tout4->SetBranchAddress("pho1y",&pho1y);
    Tout4->SetBranchAddress("pho2y",&pho2y);
    Tout4->SetBranchAddress("pho1eta",&pho1eta);
    Tout4->SetBranchAddress("pho1phi",&pho1phi);
    Tout4->SetBranchAddress("pho1e",&pho1e);
    Tout4->SetBranchAddress("pho2eta",&pho2eta);
    Tout4->SetBranchAddress("pho2phi",&pho2phi);
    Tout4->SetBranchAddress("pho2e",&pho2e);
    Tout4->SetBranchAddress("MetPhi",&MetPhi);

    Tout5->SetBranchAddress("bb_inv_mass",&bb_inv_mass);
    Tout5->SetBranchAddress("dipho_invmass",&dipho_invmass);
    Tout5->SetBranchAddress("invmassbbgg",&invmassbbgg);
    Tout5->SetBranchAddress("weight_nom",&weight_nom);
    Tout5->SetBranchAddress("ab",&ab);
    Tout5->SetBranchAddress("bb",&bb);
    Tout5->SetBranchAddress("isEle",&isEle);
    Tout5->SetBranchAddress("isMu",&isMu);
    Tout5->SetBranchAddress("pho1pt",&pho1pt);
    Tout5->SetBranchAddress("pho2pt",&pho2pt);
    Tout5->SetBranchAddress("pho1y",&pho1y);
    Tout5->SetBranchAddress("pho2y",&pho2y);
    Tout5->SetBranchAddress("pho1eta",&pho1eta);
    Tout5->SetBranchAddress("pho1phi",&pho1phi);
    Tout5->SetBranchAddress("pho1e",&pho1e);
    Tout5->SetBranchAddress("pho2eta",&pho2eta);
    Tout5->SetBranchAddress("pho2phi",&pho2phi);
    Tout5->SetBranchAddress("pho2e",&pho2e);
    Tout5->SetBranchAddress("MetPhi",&MetPhi);

    Tout6->SetBranchAddress("bb_inv_mass",&bb_inv_mass);
    Tout6->SetBranchAddress("dipho_invmass",&dipho_invmass);
    Tout6->SetBranchAddress("invmassbbgg",&invmassbbgg);
    Tout6->SetBranchAddress("weight_nom",&weight_nom);
    Tout6->SetBranchAddress("ab",&ab);
    Tout6->SetBranchAddress("bb",&bb);
    Tout6->SetBranchAddress("isEle",&isEle);
    Tout6->SetBranchAddress("isMu",&isMu);
    Tout6->SetBranchAddress("pho1pt",&pho1pt);
    Tout6->SetBranchAddress("pho2pt",&pho2pt);
    Tout6->SetBranchAddress("pho1y",&pho1y);
    Tout6->SetBranchAddress("pho2y",&pho2y);
    Tout6->SetBranchAddress("pho1eta",&pho1eta);
    Tout6->SetBranchAddress("pho1phi",&pho1phi);
    Tout6->SetBranchAddress("pho1e",&pho1e);
    Tout6->SetBranchAddress("pho2eta",&pho2eta);
    Tout6->SetBranchAddress("pho2phi",&pho2phi);
    Tout6->SetBranchAddress("pho2e",&pho2e);
    Tout6->SetBranchAddress("MetPhi",&MetPhi);



    double fact1 = (lumi_18*4.078)/(2.66557e+07);
    double fact2 = (lumi_18*5343)/(9.14154e+07);
    double fact3 = (lumi_18*365.34*1.3)/(1.48113e+11);
    double fact4 = (lumi_18*88.29*1.3)/(1.03533e+10);
    double fact5 = (lumi_18*0.01)/(993834);
    double fact6 = (lumi_18*0.01)/(903853);



        int nevt1;
        nevt1=Tout1->GetEntries();

        for (int i = 0; i < nevt1; i++)
        {
                Tout1->GetEntry(i);

                if ((isEle || isMu) && bb == 1) 
		{

                	if (weight_nom <= 0) continue;

                        sigmabb1 = 0.20*dipho_invmass + 1.69;
                        sigmah1 = -0.05*dipho_invmass + 16.34;
                   //     if (sigmabb1 == 0) continue;
                   //     if (sigmah1 == 0) continue;

                        chi2 = (((invmassbbgg - 125)/sigmah1)*((invmassbbgg - 125)/sigmah1) + ((bb_inv_mass - dipho_invmass)/sigmabb1)*((bb_inv_mass - dipho_invmass)/sigmabb1));          // for bkg
                        chi2 = (((invmassbbgg - 125)/16.04)*((invmassbbgg - 125)/16.04) + ((bb_inv_mass - dipho_invmass)/5.69)*((bb_inv_mass - dipho_invmass)/5.69));                // for 20 GeV
                        chi2 = (((invmassbbgg - 125)/13.98)*((invmassbbgg - 125)/13.98) + ((bb_inv_mass - dipho_invmass)/12.61)*((bb_inv_mass - dipho_invmass)/12.61));                // for 55 GeV

                        TLorentzVector b1, b2, g1, g2;

                        g1.SetPtEtaPhiE(pho1pt, pho1eta, pho1phi, pho1e);
                        g2.SetPtEtaPhiE(pho2pt, pho2eta, pho2phi, pho2e);

                        delphi_ggEt = fabs(PhiInRange((g1+g2).Phi() - MetPhi));

			tt_gjets->Fill(delphi_ggEt, weight_nom*fact1);
                        data_obs->Fill(delphi_ggEt, weight_nom*fact1);

		}	
        }

	int nevt2;
        nevt2=Tout2->GetEntries();

        for (int i = 0; i < nevt2; i++)
        {
                Tout2->GetEntry(i);

                if ((isEle || isMu) && bb == 1) 
		{
			if (weight_nom <= 0) continue;

                        sigmabb1 = 0.20*dipho_invmass + 1.69;
                        sigmah1 = -0.05*dipho_invmass + 16.34;
                   //     if (sigmabb1 == 0) continue;
                   //     if (sigmah1 == 0) continue;

                        chi2 = (((invmassbbgg - 125)/sigmah1)*((invmassbbgg - 125)/sigmah1) + ((bb_inv_mass - dipho_invmass)/sigmabb1)*((bb_inv_mass - dipho_invmass)/sigmabb1));          // for bkg
                        chi2 = (((invmassbbgg - 125)/16.04)*((invmassbbgg - 125)/16.04) + ((bb_inv_mass - dipho_invmass)/5.69)*((bb_inv_mass - dipho_invmass)/5.69));                // for 20 GeV
                        chi2 = (((invmassbbgg - 125)/13.98)*((invmassbbgg - 125)/13.98) + ((bb_inv_mass - dipho_invmass)/12.61)*((bb_inv_mass - dipho_invmass)/12.61));                // for 55 GeV

                        TLorentzVector b1, b2, g1, g2;

                        g1.SetPtEtaPhiE(pho1pt, pho1eta, pho1phi, pho1e);
                        g2.SetPtEtaPhiE(pho2pt, pho2eta, pho2phi, pho2e);

                        delphi_ggEt = fabs(PhiInRange((g1+g2).Phi() - MetPhi));

  	                dy->Fill(delphi_ggEt, weight_nom*fact2); 
                }
        }

	int nevt3;
        nevt3=Tout3->GetEntries();

        for (int i = 0; i < nevt3; i++)
        {
                Tout3->GetEntry(i);

                if ((isEle || isMu) && bb == 1) 
		{
			if (weight_nom <= 0) continue;

                        sigmabb1 = 0.20*dipho_invmass + 1.69;
                        sigmah1 = -0.05*dipho_invmass + 16.34;
                   //     if (sigmabb1 == 0) continue;
                   //     if (sigmah1 == 0) continue;

                        chi2 = (((invmassbbgg - 125)/sigmah1)*((invmassbbgg - 125)/sigmah1) + ((bb_inv_mass - dipho_invmass)/sigmabb1)*((bb_inv_mass - dipho_invmass)/sigmabb1));          // for bkg
                        chi2 = (((invmassbbgg - 125)/16.04)*((invmassbbgg - 125)/16.04) + ((bb_inv_mass - dipho_invmass)/5.69)*((bb_inv_mass - dipho_invmass)/5.69));                // for 20 GeV
                        chi2 = (((invmassbbgg - 125)/13.98)*((invmassbbgg - 125)/13.98) + ((bb_inv_mass - dipho_invmass)/12.61)*((bb_inv_mass - dipho_invmass)/12.61));                // for 55 GeV

                        TLorentzVector b1, b2, g1, g2;

                        g1.SetPtEtaPhiE(pho1pt, pho1eta, pho1phi, pho1e);
                        g2.SetPtEtaPhiE(pho2pt, pho2eta, pho2phi, pho2e);

                        delphi_ggEt = fabs(PhiInRange((g1+g2).Phi() - MetPhi));

                        tt_sl->Fill(delphi_ggEt, weight_nom*fact3); 
                }
        }

	int nevt4;
        nevt4=Tout4->GetEntries();

        for (int i = 0; i < nevt4; i++)
        {
                Tout4->GetEntry(i);

                if ((isEle || isMu) && bb == 1) 
		{
			if (weight_nom <= 0) continue;

                        sigmabb1 = 0.20*dipho_invmass + 1.69;
                        sigmah1 = -0.05*dipho_invmass + 16.34;
                   //     if (sigmabb1 == 0) continue;
                   //     if (sigmah1 == 0) continue;

                        chi2 = (((invmassbbgg - 125)/sigmah1)*((invmassbbgg - 125)/sigmah1) + ((bb_inv_mass - dipho_invmass)/sigmabb1)*((bb_inv_mass - dipho_invmass)/sigmabb1));          // for bkg
                        chi2 = (((invmassbbgg - 125)/16.04)*((invmassbbgg - 125)/16.04) + ((bb_inv_mass - dipho_invmass)/5.69)*((bb_inv_mass - dipho_invmass)/5.69));                // for 20 GeV
                        chi2 = (((invmassbbgg - 125)/13.98)*((invmassbbgg - 125)/13.98) + ((bb_inv_mass - dipho_invmass)/12.61)*((bb_inv_mass - dipho_invmass)/12.61));                // for 55 GeV

                        TLorentzVector b1, b2, g1, g2;

                        g1.SetPtEtaPhiE(pho1pt, pho1eta, pho1phi, pho1e);
                        g2.SetPtEtaPhiE(pho2pt, pho2eta, pho2phi, pho2e);

                        delphi_ggEt = fabs(PhiInRange((g1+g2).Phi() - MetPhi));

                        tt_2l2nu->Fill(delphi_ggEt, weight_nom*fact4); 
                }
        }

	int nevt5;
        nevt5=Tout5->GetEntries();

        for (int i = 0; i < nevt5; i++)
        {
                Tout5->GetEntry(i);

                if ((isEle || isMu) && bb == 1) 
		{
			if (weight_nom <= 0) continue;

                        sigmabb1 = 0.20*dipho_invmass + 1.69;
                        sigmah1 = -0.05*dipho_invmass + 16.34;
                   //     if (sigmabb1 == 0) continue;
                   //     if (sigmah1 == 0) continue;

                        chi2 = (((invmassbbgg - 125)/sigmah1)*((invmassbbgg - 125)/sigmah1) + ((bb_inv_mass - dipho_invmass)/sigmabb1)*((bb_inv_mass - dipho_invmass)/sigmabb1));          // for bkg
                        chi2 = (((invmassbbgg - 125)/16.04)*((invmassbbgg - 125)/16.04) + ((bb_inv_mass - dipho_invmass)/5.69)*((bb_inv_mass - dipho_invmass)/5.69));                // for 20 GeV
                        chi2 = (((invmassbbgg - 125)/13.98)*((invmassbbgg - 125)/13.98) + ((bb_inv_mass - dipho_invmass)/12.61)*((bb_inv_mass - dipho_invmass)/12.61));                // for 55 GeV

                        TLorentzVector b1, b2, g1, g2;

                        g1.SetPtEtaPhiE(pho1pt, pho1eta, pho1phi, pho1e);
                        g2.SetPtEtaPhiE(pho2pt, pho2eta, pho2phi, pho2e);

                        delphi_ggEt = fabs(PhiInRange((g1+g2).Phi() - MetPhi));

                        sig_20->Fill(delphi_ggEt, weight_nom*fact5); 
                }
        }

	int nevt6;
        nevt6=Tout6->GetEntries();

        for (int i = 0; i < nevt6; i++)
        {
                Tout6->GetEntry(i);

                if ((isEle || isMu) && bb == 1) 
		{
			if (weight_nom <= 0) continue;

                        sigmabb1 = 0.20*dipho_invmass + 1.69;
                        sigmah1 = -0.05*dipho_invmass + 16.34;
                   //     if (sigmabb1 == 0) continue;
                   //     if (sigmah1 == 0) continue;

                        chi2 = (((invmassbbgg - 125)/sigmah1)*((invmassbbgg - 125)/sigmah1) + ((bb_inv_mass - dipho_invmass)/sigmabb1)*((bb_inv_mass - dipho_invmass)/sigmabb1));          // for bkg
                        chi2 = (((invmassbbgg - 125)/16.04)*((invmassbbgg - 125)/16.04) + ((bb_inv_mass - dipho_invmass)/5.69)*((bb_inv_mass - dipho_invmass)/5.69));                // for 20 GeV
                        chi2 = (((invmassbbgg - 125)/13.98)*((invmassbbgg - 125)/13.98) + ((bb_inv_mass - dipho_invmass)/12.61)*((bb_inv_mass - dipho_invmass)/12.61));                // for 55 GeV

                        TLorentzVector b1, b2, g1, g2;

                        g1.SetPtEtaPhiE(pho1pt, pho1eta, pho1phi, pho1e);
                        g2.SetPtEtaPhiE(pho2pt, pho2eta, pho2phi, pho2e);

                        delphi_ggEt = fabs(PhiInRange((g1+g2).Phi() - MetPhi));

                        sig_55->Fill(delphi_ggEt, weight_nom*fact6); 
                }
        }


	data_obs->Add(dy);
        data_obs->Add(tt_sl);
        data_obs->Add(tt_2l2nu);
	

	    
	f1->cd();
        delete Tout1;
        delete f1;

        f2->cd();
        delete Tout2;
        delete f2;

        f3->cd();
        delete Tout3;
        delete f3;

        f4->cd();
        delete Tout4;
        delete f4;

	f5->cd();
        delete Tout5;
        delete f5;

        f6->cd();
        delete Tout6;
        delete f6;


        fout->cd();
        fout->Write();
        fout->Close();

	end = clock();

	timetaken = (end - start) / (double)CLOCKS_PER_SEC;

	cout << "Time taken : " << timetaken << "s" << endl;

	}
