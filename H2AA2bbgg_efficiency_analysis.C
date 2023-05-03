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
#include <algorithm>
#include <vector>


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

double event_weight;

static const int njetmx = 100;
static const int npartmx = 200;

int npfjetAK4;
 float pfjetAK4pt[njetmx], pfjetAK4y[njetmx], pfjetAK4eta[njetmx], pfjetAK4phi[njetmx], pfjetAK4mass[njetmx], pfjetAK4JEC[njetmx], pfjetAK4reso[njetmx];
 bool pfjetAK4tightID[njetmx];
 float pfjetAK4btag_CMVA[njetmx], pfjetAK4btag_CSV[njetmx], pfjetAK4btag_DeepCSV[njetmx], pfjetAK4btag_DeepCSV2[njetmx], pfjetAK4btag_DeepFlav[njetmx], pfjetAK4btag_DeepQCD[njetmx];

float miset , misphi , sumEt;

int ngenjetAK4;
  float genjetAK4pt[njetmx], genjetAK4y[njetmx], genjetAK4phi[njetmx], genjetAK4btag[njetmx], genjetAK4mass[njetmx], genjetAK4sdmass[njetmx];
  
  int nmuons;
  float muonp[njetmx], muone[njetmx], muonpt[njetmx], muony[njetmx], muoneta[njetmx], muonphi[njetmx], muondrbm[njetmx], muondz[njetmx], muonpter[njetmx], muonchi[njetmx], muonecal[njetmx], muonhcal[njetmx], muonemiso[njetmx], muonhadiso[njetmx], muontkpt03[njetmx], muontkpt05[njetmx];
  float muonposmatch[njetmx], muontrkink[njetmx], muonsegcom[njetmx], muonpfiso[njetmx], muontrkvtx[njetmx], muonhit[njetmx], muonpixhit[njetmx], muonmst[njetmx], muontrklay[njetmx], muonvalfrac[njetmx];
  int muonndf[njetmx];
  bool muonisPF[njetmx], muonisGL[njetmx], muonisTRK[njetmx];
  bool muonisGoodGL[njetmx], muonisMed[njetmx], muonisLoose[njetmx];
  
  int nelecs;
  bool elmvaid[njetmx], elmvaid_noIso[njetmx];
  float elpt[njetmx], eleta[njetmx], ely[njetmx], elphi[njetmx], ele[njetmx], elp[njetmx], eldxy[njetmx],  eldxy_sv[njetmx], eldz[njetmx], elhovere[njetmx], elqovrper[njetmx], elchi[njetmx], elemiso03[njetmx], elhadiso03[njetmx], elemiso04[njetmx], elhadiso04[njetmx], elhadisodepth03[njetmx], eltkpt03[njetmx], eltkpt04[njetmx], eleoverp[njetmx], elietaieta[njetmx], eletain[njetmx], elphiin[njetmx], elfbrem[njetmx], elchhadiso03[njetmx], elchhadiso04[njetmx], elnohits[njetmx], elmisshits[njetmx] ;
  float elchhadiso[njetmx], elneuhadiso[njetmx], elphoiso[njetmx], elpuchhadiso[njetmx], elpfiso[njetmx], elconvdist[njetmx], elconvdoct[njetmx];
  int elndf[njetmx];
  float elsupcl_eta[njetmx], elsupcl_phi[njetmx], elsupcl_rawE[njetmx];
  
  int nphotons;
  bool phomvaid[njetmx];
  float phopt[njetmx], phoe[njetmx], phoIDMVA[njetmx], phoy[njetmx], phoeta[njetmx], phophi[njetmx], phoe1by9[njetmx], phoe9by25[njetmx], phohadbyem[njetmx], photrkiso[njetmx], phoemiso[njetmx];
  float phohadiso[njetmx], phochhadiso[njetmx], phoneuhadiso[njetmx], phoPUiso[njetmx], phophoiso[njetmx], phoietaieta[njetmx];
  
  int ngenparticles;
  int genpartstatus[npartmx], genpartpdg[npartmx], genpartmompdg[npartmx];
  float genpartpt[npartmx], genparteta[npartmx], genpartphi[npartmx], genpartm[npartmx];
  bool genpartfromhard[npartmx];
  
  double weightev;
  float ept, ey, bq1pt, bq1y, bq1eta, bq1phi, bq1e, bq2pt, bq2y, bq2eta, bq2phi, bq2e, bb_inv_mass, bb_DelR, mupt, muy;
  int bjet_no, pho_no, nonbjet_no;
  float leadphopt, subleadphopt, leadphoy, leadphoeta, leadphophi, leadphoe, subleadphoy, subleadphoeta, subleadphophi, subleadphoe, leadsubleadpho_invmass, leadsubleadpho_DelR, invmassbbgg;
  float M_T, Met;

  float bb_pt, bb_eta, bb_phi, gg_pt, gg_eta, gg_phi, bbgg_pt, bbgg_eta, bbgg_phi, bb_gg_dphi, bb_DelPhi, bb_DelEta, gg_DelPhi, gg_DelEta;
  double minimum, min_bgdelr;




    
int main(int argc, char *argv[])
{

        char fOut[50];
        string inputFile=argv[3];
        string path="/home/abala/cms/CMSSW_10_5_0/src/condor_job/";
	


        if(inputFile=="HAHMHToAA_To2B2G_MA_55GeV_v2.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/eff_old_trigger/55GeV/HAHMHToAA_To2B2G_MA_55GeV_v2_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="HAHMHToAA_To2B2G_MA_40GeV_v2.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/eff_old_trigger/40GeV/HAHMHToAA_To2B2G_MA_40GeV_v2_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="HAHMHToAA_To2B2G_MA_20GeV_v2.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/eff_old_trigger/20GeV/HAHMHToAA_To2B2G_MA_20GeV_v2_%s_%s.root",argv[1],argv[2]);
        }

        else{
                cout<<"Input file does not exist"<<endl;
                exit(0);
        }



        TFile *fout = new TFile(fOut,"RECREATE");
        TTree *Tout = new TTree("Tout", "Info");
        
        
        Tout->Branch("bjet_no", &bjet_no, "bjet_no/I");
        Tout->Branch("nonbjet_no", &nonbjet_no, "nonbjet_no/I");
        Tout->Branch("bq1pt", &bq1pt, "bq1pt/F");
        Tout->Branch("bq1y", &bq1y, "bq1y/F");
        Tout->Branch("bq1eta", &bq1eta, "bq1eta/F");
        Tout->Branch("bq1phi", &bq1phi, "bq1phi/F");
        Tout->Branch("bq1e", &bq1e, "bq1e/F");
        Tout->Branch("bq2pt", &bq2pt, "bq2pt/F");
        Tout->Branch("bq2y", &bq2y, "bq2y/F");
        Tout->Branch("bq2eta", &bq2eta, "bq2eta/F");
        Tout->Branch("bq2phi", &bq2phi, "bq2phi/F");
        Tout->Branch("bq2e", &bq2e, "bq2e/F");
        Tout->Branch("pho_no", &pho_no, "pho_no/I");
        Tout->Branch("leadphopt", &leadphopt, "leadphopt/F");
        Tout->Branch("subleadphopt", &subleadphopt, "subleadphopt/F");
	Tout->Branch("leadphoy", &leadphoy, "leadphoy/F");
        Tout->Branch("subleadphoy", &subleadphoy, "subleadphoy/F");
        Tout->Branch("leadphoeta", &leadphoeta, "leadphoeta/F");
        Tout->Branch("leadphophi", &leadphophi, "leadphophi/F");
        Tout->Branch("leadphoe", &leadphoe, "leadphoe/F");
        Tout->Branch("subleadphoeta", &subleadphoeta, "subleadphoeta/F");
        Tout->Branch("subleadphophi", &subleadphophi, "subleadphophi/F");
        Tout->Branch("subleadphoe", &subleadphoe, "subleadphoe/F");
        

        TH1F *cut_flow = new TH1F("cut_flow","cut_flow",10,0,10);
        cut_flow->Sumw2();

        
  
        int TotEvt=0,count=0;
   long int processedEvt=0;
   string fileName;
   ifstream infile;
   infile.open(path+argv[3]);
   while(!infile.eof()){
    count = count+1;
   getline(infile,fileName);

   int L_lim = stof(argv[1]);
   int H_lim = stof(argv[2]);
   if(count<=L_lim)continue;
   if(count>H_lim)continue;
   TFile *f = TFile::Open(fileName.data());
   if(f==0) continue;
   

   
   TTree *T1 = (TTree*)f->Get("T1");
   
   T1->SetBranchAddress("event_weight",&event_weight);
   
   T1->SetBranchAddress("npfjetAK4",&npfjetAK4);
    T1->SetBranchAddress("pfjetAK4tightID",pfjetAK4tightID);
    T1->SetBranchAddress("pfjetAK4pt",pfjetAK4pt);
    T1->SetBranchAddress("pfjetAK4y",pfjetAK4y);
    T1->SetBranchAddress("pfjetAK4eta", pfjetAK4eta);
    T1->SetBranchAddress("pfjetAK4phi",pfjetAK4phi);
    T1->SetBranchAddress("pfjetAK4mass",pfjetAK4mass);
    T1->SetBranchAddress("pfjetAK4JEC",pfjetAK4JEC);
    T1->SetBranchAddress("pfjetAK4reso",pfjetAK4reso);
    T1->SetBranchAddress("pfjetAK4btag_CMVA",pfjetAK4btag_CMVA);
    T1->SetBranchAddress("pfjetAK4btag_CSV",pfjetAK4btag_CSV);
    T1->SetBranchAddress("pfjetAK4btag_DeepCSV",pfjetAK4btag_DeepCSV);
    T1->SetBranchAddress("pfjetAK4btag_DeepFlav",pfjetAK4btag_DeepFlav);
    T1->SetBranchAddress("pfjetAK4btag_DeepQCD",pfjetAK4btag_DeepQCD);
    
     T1->SetBranchAddress("nelecs",&nelecs);
//  T1->SetBranchAddress("elsupcl_eta",&elsupcl_eta);
//  T1->SetBranchAddress("elsupcl_phi",&elsupcl_phi);
//  T1->SetBranchAddress("elsupcl_rawE",&elsupcl_rawE);
  T1->SetBranchAddress("elpt",elpt);
  T1->SetBranchAddress("eleta",eleta);
  T1->SetBranchAddress("ely",ely);
  T1->SetBranchAddress("elphi",elphi);
  T1->SetBranchAddress("elp",elp);
  T1->SetBranchAddress("ele",ele);
  T1->SetBranchAddress("elmvaid",elmvaid);
  T1->SetBranchAddress("elmvaid_noIso",elmvaid_noIso);
  T1->SetBranchAddress("eldxy",eldxy);
  T1->SetBranchAddress("eldxy_sv",eldxy_sv);
  T1->SetBranchAddress("eldz",eldz);
  T1->SetBranchAddress("elhovere",elhovere);
  T1->SetBranchAddress("elchi",elchi);
  T1->SetBranchAddress("elndf",elndf);
  T1->SetBranchAddress("eltkpt03",eltkpt03);
//  T1->SetBranchAddress("elemiso03",elemiso03);
//  T1->SetBranchAddress("elhadiso03",elhadiso03);
  T1->SetBranchAddress("eltkpt04",eltkpt04);
  T1->SetBranchAddress("elemiso04",elemiso04);
  T1->SetBranchAddress("elhadiso04",elhadiso04);
  T1->SetBranchAddress("eletain",eletain);
  T1->SetBranchAddress("elphiin",elphiin);
//T1->SetBranchAddress("elsceta",elsceta);
  T1->SetBranchAddress("elfbrem",elfbrem);
  T1->SetBranchAddress("elhadisodepth03",elhadisodepth03);
  T1->SetBranchAddress("eleoverp",eleoverp);
  T1->SetBranchAddress("elietaieta",elietaieta);
  T1->SetBranchAddress("elmisshits",elmisshits);
  T1->SetBranchAddress("elchhadiso",elchhadiso);
  T1->SetBranchAddress("elneuhadiso",elneuhadiso);
  T1->SetBranchAddress("elphoiso",elphoiso);
  T1->SetBranchAddress("elpuchhadiso",elpuchhadiso);
  T1->SetBranchAddress("elpfiso",elpfiso);
  T1->SetBranchAddress("elconvdist",elconvdist);
  T1->SetBranchAddress("elconvdoct",elconvdoct);
  
  T1->SetBranchAddress("nmuons",&nmuons);
  T1->SetBranchAddress("muonpt",muonpt);
  T1->SetBranchAddress("muony",muony);
  T1->SetBranchAddress("muoneta",muoneta);
  T1->SetBranchAddress("muonphi",muonphi);
  T1->SetBranchAddress("muone",muone);
  T1->SetBranchAddress("muonpfiso",muonpfiso);
  T1->SetBranchAddress("muonisPF",muonisPF);
  T1->SetBranchAddress("muonisMed",muonisMed);
  
  T1->SetBranchAddress("PFMET",&miset) ;
    T1->SetBranchAddress("PFMETPhi",&misphi) ;
    T1->SetBranchAddress("sumEt",&sumEt);
    
    T1->SetBranchAddress("nphotons",&nphotons);
    T1->SetBranchAddress("phopt",phopt);
    T1->SetBranchAddress("phoeta",phoeta);
    T1->SetBranchAddress("phoy",phoy);
    T1->SetBranchAddress("phophi",phophi);
    T1->SetBranchAddress("phoe",phoe);
    T1->SetBranchAddress("phoIDMVA",phoIDMVA);
    T1->SetBranchAddress("phochhadiso",phochhadiso);
    T1->SetBranchAddress("phoneuhadiso",phoneuhadiso);
    T1->SetBranchAddress("phophoiso",phophoiso);
    T1->SetBranchAddress("phohadbyem",phohadbyem);
    T1->SetBranchAddress("phoietaieta",phoietaieta);
    T1->SetBranchAddress("phoe9by25",phoe9by25);
    
    T1->SetBranchAddress("ngenparticles",&ngenparticles);
    T1->SetBranchAddress("genpartstatus",genpartstatus);
    T1->SetBranchAddress("genpartpdg",genpartpdg);
    T1->SetBranchAddress("genpartmompdg",genpartmompdg);
    T1->SetBranchAddress("genpartfromhard",genpartfromhard);
    T1->SetBranchAddress("genpartpt",genpartpt);
    T1->SetBranchAddress("genparteta",genparteta);
    T1->SetBranchAddress("genpartphi",genpartphi);
    T1->SetBranchAddress("genpartm",genpartm);
    
  
    
    
    int nevents;
    nevents=T1->GetEntries();
    
    for(int iev=0; iev<nevents; iev++)
    {
            
            T1->GetEntry(iev);
            
            cut_flow->Fill(0);



	//--------------------------------------------------------------------------------------------------------------//


                        int indx = 0;
                        int indx_totjet = npfjetAK4;

                        for(int i=0; i<npfjetAK4; i++)
                        {
                                if(!pfjetAK4tightID[i]) continue;
                                pfjetAK4pt[i] *= pfjetAK4JEC[i];
                                pfjetAK4mass[i] *= pfjetAK4JEC[i];
                                pfjetAK4pt[i] *= (1.+pfjetAK4reso[i]);
                                if(fabs(pfjetAK4y[i]) > 2.5) continue;
                                if(pfjetAK4pt[i]<20) continue;
                                if (pfjetAK4btag_DeepFlav[i] < 0.3033) continue;

                                pfjetAK4pt[indx] = pfjetAK4pt[i];
                                pfjetAK4y[indx] = pfjetAK4y[i];
                                pfjetAK4eta[indx] = pfjetAK4eta[i];
                                pfjetAK4phi[indx] = pfjetAK4phi[i];
                                pfjetAK4mass[indx] = pfjetAK4mass[i];
                                pfjetAK4btag_DeepFlav[indx] = pfjetAK4btag_DeepFlav[i];

                                if (++indx >= njetmx) break;
                        }

                        npfjetAK4 = indx;
                        int indx_b = npfjetAK4;
                        int indx_nonb = (indx_totjet - indx_b);

                        int indx1 = 0;
                        for (int i=0; i<indx_nonb; i++)
                        {
                                if(!pfjetAK4tightID[i]) continue;
                                pfjetAK4pt[i] *= pfjetAK4JEC[i];
                                pfjetAK4mass[i] *= pfjetAK4JEC[i];
                                pfjetAK4pt[i] *= (1.+pfjetAK4reso[i]);
                                if(fabs(pfjetAK4y[i]) > 2.5) continue;
                                if(pfjetAK4pt[i]<20) continue;

                                if (++indx1 >= njetmx) break;
                         }

                        if (npfjetAK4 < 2) continue;
                        cut_flow->Fill(1);	


                        for (int i=0; i<npfjetAK4; i++)
                        {
                               for (int j=i+1; j<npfjetAK4; j++)
                               {
                                        if(pfjetAK4btag_DeepFlav[i]<pfjetAK4btag_DeepFlav[j])
                                        {
                                                 float a = pfjetAK4btag_DeepFlav[i];
                                                 pfjetAK4btag_DeepFlav[i] = pfjetAK4btag_DeepFlav[j];
                                                 pfjetAK4btag_DeepFlav[j] = a;

                                                 a = pfjetAK4pt[i];
                                                 pfjetAK4pt[i] = pfjetAK4pt[j];
                                                 pfjetAK4pt[j] = a;

                                                 a = pfjetAK4y[i];
                                                 pfjetAK4y[i] = pfjetAK4y[j];
                                                 pfjetAK4y[j] = a;

                                                 a = pfjetAK4eta[i];
                                                 pfjetAK4eta[i] = pfjetAK4eta[j];
                                                 pfjetAK4eta[j] = a;

                                                 a = pfjetAK4phi[i];
                                                 pfjetAK4phi[i] = pfjetAK4phi[j];
                                                 pfjetAK4phi[j] = a;

                                                 a = pfjetAK4mass[i];
                                                 pfjetAK4mass[i] = pfjetAK4mass[j];
                                                 pfjetAK4mass[j] = a;
                                         }
                                }
                         }

                        TLorentzVector leadbvec, subleadbvec;

                        if (pfjetAK4pt[0] >= pfjetAK4pt[1])
                        {
                        leadbvec.SetPtEtaPhiM(pfjetAK4pt[0], pfjetAK4eta[0], pfjetAK4phi[0], pfjetAK4mass[0]);
                        subleadbvec.SetPtEtaPhiM(pfjetAK4pt[1], pfjetAK4eta[1], pfjetAK4phi[1], pfjetAK4mass[1]);
                        }

                        else
                        {
                        subleadbvec.SetPtEtaPhiM(pfjetAK4pt[0], pfjetAK4eta[0], pfjetAK4phi[0], pfjetAK4mass[0]);
                        leadbvec.SetPtEtaPhiM(pfjetAK4pt[1], pfjetAK4eta[1], pfjetAK4phi[1], pfjetAK4mass[1]);
                        }


//---------------------------------------------------------------------------------------------------------------------------//
				

		indx = 0;
		for(int i=0; i<nphotons; i++)
		{

                //        if (phopt[i] < 20) continue;
                        if (fabs(phoeta[i]) > 2.5) continue;
			if (fabs(phoeta[i]) > 1.44 && fabs(phoeta[i]) < 1.57) continue;
                        if (phohadbyem[i] >= 0.08) continue;
                        if (phoIDMVA[i] <= -0.9) continue;

                        TLorentzVector phovec;
                        phovec.SetPtEtaPhiE(phopt[i], phoeta[i], phophi[i], phoe[i]);

                        if ((phoe9by25[i] > 0.8 || phochhadiso[i] < 20 || phochhadiso[i]/phoe[i] < 0.3))
                        {
                                phopt[indx] = phopt[i];
                                phoy[indx] = phoy[i];
                                phoeta[indx] = phoeta[i];
                                phophi[indx] = phophi[i];
                                phoe[indx] = phoe[i];

                                if (++indx >= njetmx) break;
                        }
		}
		
		nphotons = indx;
		if (nphotons < 2) continue;
		cut_flow->Fill(2);
		
		for (int i=0; i<nphotons; i++)
		{
			for (int j=i+1; j<nphotons; j++)
			{
				if(phopt[i]<phopt[j])
				{
					float a = phopt[i];
                                        phopt[i] = phopt[j];
                                        phopt[j] = a;
                    
                                        a = phoy[i];
                                        phoy[i] = phoy[j];
                                        phoy[j] = a;   
                    
                                        a = phoeta[i];
                                        phoeta[i] = phoeta[j];
                                        phoeta[j] = a;
                    
                                        a = phophi[i];
                                        phophi[i] = phophi[j];
                                        phophi[j] = a;
                    
                                        a = phoe[i];
                                        phoe[i] = phoe[j];
                                        phoe[j] = a;
				}
			}
		}	
		
	       	TLorentzVector leadphovec, subleadphovec;
		leadphovec.SetPtEtaPhiE(phopt[0], phoeta[0], phophi[0], phoe[0]);
		subleadphovec.SetPtEtaPhiE(phopt[1], phoeta[1], phophi[1], phoe[1]);

                if (leadphovec.Pt() < 30) continue;
                if (subleadphovec.Pt() < 18) continue;

                cut_flow->Fill(3);

         
		

	//------------------------------------------------------------------------------------------------------------//	


                bjet_no = npfjetAK4;
                nonbjet_no = indx1;
                bq1pt = leadbvec.Pt();
                bq1y = leadbvec.Rapidity();
                bq1eta = leadbvec.Eta();
                bq1phi = leadbvec.Phi();
                bq1e = leadbvec.Energy();
                bq2pt = subleadbvec.Pt();
                bq2y = subleadbvec.Rapidity();
                bq2eta = subleadbvec.Eta();
                bq2phi = subleadbvec.Phi();
                bq2e = subleadbvec.Energy();
	
                pho_no = nphotons;
                leadphopt = leadphovec.Pt();
		leadphoy = leadphovec.Rapidity();
                leadphoeta = leadphovec.Eta();
                leadphophi = leadphovec.Phi();
                leadphoe = leadphovec.Energy();
                subleadphopt = subleadphovec.Pt();
		subleadphoy = subleadphovec.Rapidity();
                subleadphoeta = subleadphovec.Eta();
                subleadphophi = subleadphovec.Phi();
                subleadphoe = subleadphovec.Energy();


		Tout->Fill();
	}
	
	    f->cd();
            delete T1;
            delete f;

}
         infile.close();
        fout->cd();
        fout->Write();
        fout->Close();

	return 0;

}				
