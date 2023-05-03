#include<iostream>
#include<vector>
#include<fstream>
#include<string>
//#include "lester_mt2_bisect.h"
//#include "lester_mt2_bisect.cpp"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TObject.h"

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

int main(int argc, char *argv[]){
	//Declaration of the leaf types for reco objets
	cout<<"Program started"<<endl;

	char fOut[50];
	string inputFile=argv[3];
	string path="/home/abala/cms/CMSSW_10_5_0/src/condor_job/";


	if(inputFile=="TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_Autumn18.txt"){
		sprintf(fOut,"TTJets_sl_%s_%s.root",argv[1],argv[2]);
	}

	else{
		cout<<"Input file does not exist"<<endl;
		exit(0);
	}

    TFile *fout = new TFile(fOut,"RECREATE");
    
    const int njetmax = 30;
    
    
    TH1F *el_pt = new TH1F("el_pt","el_pt",100,0,500);
    el_pt->Sumw2();

    TH1F *el_eta = new TH1F("el_eta","el_eta",100,-3,3);                           // Reco level info.
    el_eta->Sumw2();

    TH1F *jet_pt = new TH1F("jet_pt","jet_pt",100,0,500);
    jet_pt->Sumw2();

    TH1F *jet_eta = new TH1F("jet_eta","jet_eta",100,-3,3);
    jet_eta->Sumw2();

    TH1F *met = new TH1F("met","met",100,0,500);
    met->Sumw2();

    TH1F *el_gen_pt = new TH1F("el_gen_pt","el_gen_pt",100,0,500);
    el_gen_pt->Sumw2();

    TH1F *el_gen_eta = new TH1F("el_gen_eta","el_gen_eta",100,-3,3);
    el_gen_eta->Sumw2();

    TH1F *jet_gen_pt = new TH1F("jet_gen_pt","jet_gen_pt",100,0,500);                   // GEN level info.
    jet_gen_pt->Sumw2();

    TH1F *jet_gen_eta = new TH1F("jet_gen_eta","jet_gen_eta",100,-3,3);
    jet_gen_eta->Sumw2();

    TH1F *met_gen = new TH1F("met_gen","met_gen",100,0,500);
    met_gen->Sumw2();

    TH1F *t_gen_pt = new TH1F("t_gen_pt","t_gen_pt",100,0,500);
    t_gen_pt->Sumw2();

    TH1F *t_gen_eta = new TH1F("t_gen_eta","t_gen_eta",100,-3,3);
    t_gen_eta->Sumw2();

    TH1F *no_gen_jet = new TH1F("no_gen_jet","no_gen_jet",20,0,20);
    no_gen_jet->Sumw2();

    TH1F *no_rec_jet = new TH1F("no_rec_jet","no_rec_jet",20,0,20);
    no_rec_jet->Sumw2();

    TH1F *no_gen_el = new TH1F("no_gen_el","no_gen_el",10,0,10);
    no_gen_el->Sumw2();

    TH1F *no_rec_el = new TH1F("no_rec_el","no_rec_el",10,0,10);
    no_rec_el->Sumw2();
    
    TH1F *trmass_rec = new TH1F("trmass_rec", "trmass_rec", 100, 0., 300.);
    trmass_rec->Sumw2();
    
    TH1F *trmass_gen = new TH1F("trmass_gen", "trmass_gen", 100, 0., 300.);
    trmass_gen->Sumw2();
    


    
    int TotEvt=0,count=0;
   long int processedEvt=0;
   double weight = 1;
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
   
   cout<<fileName<<endl;

   cout<<"processedEvt="<<processedEvt<<endl;


   Int_t           el_n_reco;
   vector<double>  *el_E_reco=0;
   vector<double>  *el_px_reco=0;
   vector<double>  *el_py_reco=0;
   vector<double>  *el_pz_reco=0;
   vector<double>  *el_missingInnerHits_reco=0;
   vector<double>  *el_passConversionVeto_reco=0;
   vector<double>  *el_mvaEleIDFall17isoV2wp80_reco=0;
   vector<double>  *el_byTightPFbasedCombinedRelativeDeltaBetaCorrdR03_reco=0;
   vector<double>  *el_dxy_reco=0;
   vector<double>  *el_dz_reco=0;

   Int_t           jet_n_reco;
   vector<double>  *jet_E_reco=0;                             // Reco level info.
   vector<double>  *jet_px_reco=0;
   vector<double>  *jet_py_reco=0;
   vector<double>  *jet_pz_reco=0;

   double     MET_E_reco;
   double     MET_phi_reco;

   Int_t           el_n_gen;
   vector<double>  *el_E_gen=0;
   vector<double>  *el_px_gen=0;
   vector<double>  *el_py_gen=0;
   vector<double>  *el_pz_gen=0;

   Int_t           jet_n_gen;                              // GEN level info.
   vector<double>  *jet_E_gen=0;
   vector<double>  *jet_px_gen=0;
   vector<double>  *jet_py_gen=0;
   vector<double>  *jet_pz_gen=0;

   double     MET_E_gen;
   double     MET_phi_gen;
   
   Int_t           t_n_gen;
   vector<double>  *t_E_gen=0;
   vector<double>  *t_px_gen=0;
   vector<double>  *t_py_gen=0;                            // GEN level info.
   vector<double>  *t_pz_gen=0;



   TTree *reco_tree;
   TTree *gen_tree;

   reco_tree = (TTree*)f->Get("miniaodsim2custom/reco");
   gen_tree = (TTree*)f->Get("miniaodsim2custom/gen");

   reco_tree->SetBranchAddress("el_n_reco", &el_n_reco);
   reco_tree->SetBranchAddress("el_E_reco", &el_E_reco);
   reco_tree->SetBranchAddress("el_px_reco", &el_px_reco);
   reco_tree->SetBranchAddress("el_py_reco", &el_py_reco);
   reco_tree->SetBranchAddress("el_pz_reco", &el_pz_reco);
   reco_tree->SetBranchAddress("el_missingInnerHits_reco", &el_missingInnerHits_reco);
   reco_tree->SetBranchAddress("el_passConversionVeto_reco", &el_passConversionVeto_reco);
   reco_tree->SetBranchAddress("el_mvaEleIDFall17isoV2wp80_reco", &el_mvaEleIDFall17isoV2wp80_reco);
   reco_tree->SetBranchAddress("el_byTightPFbasedCombinedRelativeDeltaBetaCorrdR03_reco", &el_byTightPFbasedCombinedRelativeDeltaBetaCorrdR03_reco);
   reco_tree->SetBranchAddress("el_dxy_reco", &el_dxy_reco);
   reco_tree->SetBranchAddress("el_dz_reco", &el_dz_reco);

   reco_tree->SetBranchAddress("jet_n_reco", &jet_n_reco);                             // RECO
   reco_tree->SetBranchAddress("jet_E_reco", &jet_E_reco);
   reco_tree->SetBranchAddress("jet_px_reco", &jet_px_reco);
   reco_tree->SetBranchAddress("jet_py_reco", &jet_py_reco);
   reco_tree->SetBranchAddress("jet_pz_reco", &jet_pz_reco);

   reco_tree->SetBranchAddress("MET_E_reco", &MET_E_reco);
   reco_tree->SetBranchAddress("MET_phi_reco", &MET_phi_reco);

   gen_tree->SetBranchAddress("el_n_gen", &el_n_gen);
   gen_tree->SetBranchAddress("el_E_gen", &el_E_gen);
   gen_tree->SetBranchAddress("el_px_gen", &el_px_gen);
   gen_tree->SetBranchAddress("el_py_gen", &el_py_gen);
   gen_tree->SetBranchAddress("el_pz_gen", &el_pz_gen);

   gen_tree->SetBranchAddress("jet_n_gen", &jet_n_gen);                          
   gen_tree->SetBranchAddress("jet_E_gen", &jet_E_gen);
   gen_tree->SetBranchAddress("jet_px_gen", &jet_px_gen);
   gen_tree->SetBranchAddress("jet_py_gen", &jet_py_gen);                         // GEN 
   gen_tree->SetBranchAddress("jet_pz_gen", &jet_pz_gen);

   gen_tree->SetBranchAddress("MET_E_gen", &MET_E_gen);
   gen_tree->SetBranchAddress("MET_phi_gen", &MET_phi_gen);

   gen_tree->SetBranchAddress("t_n_gen", &t_n_gen);
   gen_tree->SetBranchAddress("t_E_gen", &t_E_gen);
   gen_tree->SetBranchAddress("t_px_gen", &t_px_gen);
   gen_tree->SetBranchAddress("t_py_gen", &t_py_gen);
   gen_tree->SetBranchAddress("t_pz_gen", &t_pz_gen);


   int ent_rec = reco_tree->GetEntries();          // RECO  
   for(int i=0; i<ent_rec; ++i)              // event loop starts  
   {
	   reco_tree->GetEntry(i);

	   if((TotEvt+i)%10000==0)
	   {
		   cout<<i<<" Events Analysed out of "<<ent_rec<<endl;
	   }

	   int nrecel=0; int nrecjet=0; int count_1=0;
	   TLorentzVector ele[njetmax];
       
	   for(int iel=0; iel < el_px_reco->size(); ++iel)
	   {
		   if (el_mvaEleIDFall17isoV2wp80_reco->at(iel)==0) continue;
		   if (el_passConversionVeto_reco->at(iel)==0) continue;
		   if (el_missingInnerHits_reco->at(iel)>1) continue;
		   if(fabs(el_dxy_reco->at(iel)) > 0.045 || fabs(el_dz_reco->at(iel))>0.2) continue;
		   if (el_byTightPFbasedCombinedRelativeDeltaBetaCorrdR03_reco->at(iel)==0) continue;
                  nrecel++;
                  ele[nrecel].SetPxPyPzE(el_px_reco->at(iel),el_py_reco->at(iel),el_pz_reco->at(iel),fabs(el_E_reco->at(iel)));
	   }
           if(ele[nrecel].Pt()<20 || fabs(ele[nrecel].Eta())>2.5) continue;
	   no_rec_el->Fill(nrecel);
           if (nrecel==1)
           {
                 el_pt->Fill(ele[nrecel].Pt());
                 el_eta->Fill(ele[nrecel].Eta());
           }


	   for(int ijet=0; ijet<jet_px_reco->size(); ++ijet)
	   {
		   TLorentzVector je(jet_px_reco->at(ijet),jet_py_reco->at(ijet),jet_pz_reco->at(ijet),fabs(jet_E_reco->at(ijet)));
		   if(je.Pt()<20 || fabs(je.Eta())>2.5) continue;
		   jet_pt->Fill(je.Pt());
                   jet_eta->Fill(je.Eta());
		   nrecjet++;
	   }
	   no_rec_jet->Fill(nrecjet);

	   if(MET_E_reco>15)
	   {
		   met->Fill(MET_E_reco);
	   }
	   
	   if (nrecel==1)
	   {
		   double dphi = PhiInRange(ele[nrecel].Phi()-MET_phi_reco);
		   double Mt_rec = sqrt(2*(ele[nrecel].Pt())*(MET_E_reco)*(1-cos(dphi))); 
                   trmass_rec->Fill(Mt_rec);
	   }
	   
   }           // end of event loop


   int ent_gen = gen_tree->GetEntries();         // GEN    
   for(int j=0; j<ent_gen; ++j)              // event loop starts     
   {
           gen_tree->GetEntry(j);

           if((TotEvt+j)%10000==0)
           {
                   cout<<j<<" Events Analysed out of "<<ent_gen<<endl;
           }

	   int ngenel=0; int ngenjet=0; int count_2=0;
	   TLorentzVector eleg[njetmax];
	   
           for(int ielg=0; ielg < el_px_gen->size(); ++ielg)
           {
		   ngenel++;
                   eleg[ngenel].SetPxPyPzE(el_px_gen->at(ielg),el_py_gen->at(ielg),el_pz_gen->at(ielg),fabs(el_E_gen->at(ielg)));
           }
           if(eleg[ngenel].Pt()<20 || fabs(eleg[ngenel].Eta())>2.5) continue;
           no_gen_el->Fill(ngenel);
           if(ngenel==1)
           {
                 el_gen_pt->Fill(eleg[ngenel].Pt());
                 el_gen_eta->Fill(eleg[ngenel].Eta());
           }

           for(int ijetg=0; ijetg<jet_px_gen->size(); ++ijetg)
           {
                   TLorentzVector jeg(jet_px_gen->at(ijetg),jet_py_gen->at(ijetg),jet_pz_gen->at(ijetg),fabs(jet_E_gen->at(ijetg)));
                   if(jeg.Pt()<20 || fabs(jeg.Eta())>2.5) continue;
                   jet_gen_pt->Fill(jeg.Pt());
                   jet_gen_eta->Fill(jeg.Eta());
		   ngenjet++;
           }
	   no_gen_jet->Fill(ngenjet);

           if(MET_E_gen>15)
           {
                   met_gen->Fill(MET_E_gen);
           }
           
           if (ngenel==1)
	   {
		    double dphig = PhiInRange(eleg[ngenel].Phi()-MET_phi_gen);
		    double Mt_gen = sqrt(2*(eleg[ngenel].Pt())*(MET_E_gen)*(1-cos(dphig))); 
                    trmass_gen->Fill(Mt_gen);
	   }

	   for(int itg=0; itg < t_px_gen->size(); ++itg)
           {
                   TLorentzVector tg(t_px_gen->at(itg),t_py_gen->at(itg),t_pz_gen->at(itg),fabs(t_E_gen->at(itg)));
                   if(tg.Pt()<20 || fabs(tg.Eta())>2.5) continue;
                   t_gen_pt->Fill(tg.Pt());
                   t_gen_eta->Fill(tg.Eta());
           }
   }          // end of event loop 



   f->Close();
}

infile.close();
fout->cd();

el_pt->Write();
el_eta->Write();
jet_pt->Write();
jet_eta->Write();
met->Write();

el_gen_pt->Write();
el_gen_eta->Write();
jet_gen_pt->Write();
jet_gen_eta->Write();
met_gen->Write();

t_gen_pt->Write();
t_gen_eta->Write();

no_rec_el->Write();
no_rec_jet->Write();
no_gen_el->Write();
no_gen_jet->Write();

trmass_gen->Write();
trmass_rec->Write();


return 0;
}
