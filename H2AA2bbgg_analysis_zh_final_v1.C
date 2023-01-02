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
#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cpp"


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


float* Muon_SF(TFile *file_mu_sf, string id, float pt, float eta){

        char name[100];

        sprintf(name,"NUM_%sID_DEN_TrackerMuons_abseta_pt",id.c_str());
        TH2F *h_SF = (TH2F*)file_mu_sf->Get(name);
        sprintf(name,"NUM_%sID_DEN_TrackerMuons_abseta_pt_stat",id.c_str());
        TH2F *h_SF_stat = (TH2F*)file_mu_sf->Get(name);
        sprintf(name,"NUM_%sID_DEN_TrackerMuons_abseta_pt_syst",id.c_str());
        TH2F *h_SF_sys = (TH2F*)file_mu_sf->Get(name);

        int eta_bin_id = h_SF->GetXaxis()->FindBin(fabs(eta));
        int pt_bin_id = h_SF->GetYaxis()->FindBin(pt);

        float sf, sf_stat, sf_sys, sf_err;

        if(eta_bin_id>0 && eta_bin_id<=(h_SF->GetNbinsX()) && pt_bin_id>0 && pt_bin_id<=(h_SF->GetNbinsY())){
                sf = h_SF->GetBinContent(eta_bin_id,pt_bin_id);
                sf_err = h_SF->GetBinError(eta_bin_id,pt_bin_id);
                sf_stat = h_SF_stat->GetBinContent(eta_bin_id,pt_bin_id);
                sf_sys = h_SF_sys->GetBinContent(eta_bin_id,pt_bin_id);
        }else{
                sf = 1;
                sf_err = 0;
                sf_stat = sf_sys = 1;
                }

        static float sfvalues[5];
        sfvalues[0] = sf;
        sfvalues[1] = sf+sf_err;
        sfvalues[2] = sf-sf_err;
        sfvalues[3] = sf_stat;
        sfvalues[4] = sf_sys;

        return sfvalues;
}


float* Electron_SF(TFile *file_el_sf, float pt, float eta){

        char name[100];

        TH2F *h_SF = (TH2F*)file_el_sf->Get("EGamma_SF2D");
        TH2F *h_SF_statData = (TH2F*)file_el_sf->Get("statData");
        TH2F *h_SF_statMC = (TH2F*)file_el_sf->Get("statMC");
        TH2F *h_SF_altBkgModel = (TH2F*)file_el_sf->Get("altBkgModel");
        TH2F *h_SF_altSignalModel = (TH2F*)file_el_sf->Get("altSignalModel");
        TH2F *h_SF_altMCEff = (TH2F*)file_el_sf->Get("altMCEff");
        TH2F *h_SF_altTagSelection = (TH2F*)file_el_sf->Get("altTagSelection");

        int eta_bin_id = h_SF->GetXaxis()->FindBin(eta);
        int pt_bin_id = h_SF->GetYaxis()->FindBin(pt);

        static float sfvalues[9] = {-100,-100,-100,-100,-100,-100,-100,-100,-100};

        if(eta_bin_id>0 && eta_bin_id<=(h_SF->GetNbinsX()) && pt_bin_id>0 && pt_bin_id<=(h_SF->GetNbinsY())){
                sfvalues[0] = h_SF->GetBinContent(eta_bin_id,pt_bin_id);
                sfvalues[1] = sfvalues[0]+h_SF->GetBinError(eta_bin_id,pt_bin_id);
                sfvalues[2] = sfvalues[0]-h_SF->GetBinError(eta_bin_id,pt_bin_id);
                sfvalues[3] = h_SF_statData->GetBinContent(eta_bin_id,pt_bin_id);
                sfvalues[4] = h_SF_statMC->GetBinContent(eta_bin_id,pt_bin_id);
                sfvalues[5] = h_SF_altBkgModel->GetBinContent(eta_bin_id,pt_bin_id);
                sfvalues[6] = h_SF_altSignalModel->GetBinContent(eta_bin_id,pt_bin_id);
                sfvalues[7] = h_SF_altMCEff->GetBinContent(eta_bin_id,pt_bin_id);
                sfvalues[8] = h_SF_altTagSelection->GetBinContent(eta_bin_id,pt_bin_id);
        }
        else{
                sfvalues[0] = sfvalues[1] = sfvalues[2] = 1;
                sfvalues[3] = sfvalues[4] = sfvalues[5] = sfvalues[6] = sfvalues[7] = sfvalues[8] = 0;
                }

        return sfvalues;
}


float* Get_PU_Weights(TFile *file_pu_ratio, int npu){

        TH1F *h_data = (TH1F*)file_pu_ratio->Get("pileup_weight");
        TH1F *h_data_plus = (TH1F*)file_pu_ratio->Get("pileup_plus_weight");
        TH1F *h_data_minus = (TH1F*)file_pu_ratio->Get("pileup_minus_weight");

        int bin_id = h_data->FindBin(npu);

        static float puweight[3] = {0,0,0};
        if(bin_id>=0 && bin_id<100){
                puweight[0] = h_data->GetBinContent(bin_id);
                puweight[1] = h_data_plus->GetBinContent(bin_id);
                puweight[2] = h_data_minus->GetBinContent(bin_id);
        }
        return puweight;
}



  static const int njetmx = 100;
  static const int njetmxAK8 =100;
  static const int npartmx = 100;
  static const int nconsmax = 1000;
  static const int njetconsmax = 3;
  static const int ngenjetAK8mx =100;

  BTagCalibration calib_deepcsv, calib_deepflav;
  BTagCalibrationReader reader_deepcsv, reader_deepflav;

  double Generator_weight;
  double weights[njetmx];

  double prefiringweight, prefiringweightup, prefiringweightdown;

  bool isMC;
  bool isFastSIM;
  bool isDATA;

  int nPFJetAK4, PFJetAK4_hadronflav[njetmx], PFJetAK4_partonflav[njetmx];
  float PFJetAK4_pt[njetmx], PFJetAK4_eta[njetmx], PFJetAK4_y[njetmx], PFJetAK4_phi[njetmx], PFJetAK4_mass[njetmx], PFJetAK4_PUID[njetmx];
  float PFJetAK4_btag_DeepCSV[njetmx], PFJetAK4_btag_DeepFlav[njetmx];
  bool PFJetAK4_jetID[njetmx], PFJetAK4_jetID_tightlepveto[njetmx];
  float PFJetAK4_btag_DeepCSV_SF[njetmx], PFJetAK4_btag_DeepCSV_SF_up[njetmx], PFJetAK4_btag_DeepCSV_SF_dn[njetmx];
  float PFJetAK4_btag_DeepFlav_SF[njetmx], PFJetAK4_btag_DeepFlav_SF_up[njetmx], PFJetAK4_btag_DeepFlav_SF_dn[njetmx];
  float PFJetAK4_reso[njetmx], PFJetAK4_resoup[njetmx], PFJetAK4_resodn[njetmx];
  float PFJetAK4_JEC[njetmx];
  float PFJetAK4_jesup_Total[njetmx], PFJetAK4_jesdn_Total[njetmx];

  double LHE_weight;

  int nMuon;
  float Muon_minchiso[njetmx], Muon_minnhiso[njetmx], Muon_minphiso[njetmx], Muon_minisoall[njetmx];
  float Muon_charge[njetmx], Muon_p[njetmx], Muon_pt[njetmx], Muon_eta[njetmx], Muon_phi[njetmx], Muon_e[njetmx], Muon_dz[njetmx], Muon_ip3d[njetmx], Muon_ptErr[njetmx], Muon_chi[njetmx], Muon_ecal[njetmx], Muon_hcal[njetmx];
  float Muon_posmatch[njetmx], Muon_trkink[njetmx], Muon_segcom[njetmx], Muon_pfiso[njetmx], Muon_dxy[njetmx], Muon_dxyErr[njetmx], Muon_hit[njetmx], Muon_pixhit[njetmx], Muon_mst[njetmx], Muon_trklay[njetmx], Muon_valfrac[njetmx],Muon_dxy_sv[njetmx];
  int Muon_ndf[njetmx];
  bool Muon_isPF[njetmx], Muon_isGL[njetmx], Muon_isTRK[njetmx];
  bool Muon_isGoodGL[njetmx], Muon_isTight[njetmx], Muon_isHighPt[njetmx], Muon_isHighPttrk[njetmx], Muon_isMed[njetmx], Muon_isMedPr[njetmx], Muon_isLoose[njetmx], Muon_TightID[njetmx];
  float Muon_corrected_pt[njetmx], Muon_correctedUp_pt[njetmx], Muon_correctedDown_pt[njetmx];

  int nElectron;
  bool Electron_mvaid_Fallv2WP90[njetmx], Electron_mvaid_Fallv2WP90_noIso[njetmx];
  bool Electron_mvaid_Fallv2WP80[njetmx], Electron_mvaid_Fallv2WP80_noIso[njetmx];
  float Electron_charge[njetmx], Electron_pt[njetmx], Electron_eta[njetmx], Electron_phi[njetmx], Electron_e[njetmx], Electron_e_ECAL[njetmx], Electron_p[njetmx];
  float Electron_dxy[njetmx],  Electron_dxyErr[njetmx], Electron_dxy_sv[njetmx], Electron_dz[njetmx], Electron_dzErr[njetmx], Electron_ip3d[njetmx];
  float Electron_hovere[njetmx], Electron_qovrper[njetmx], Electron_chi[njetmx]; //Electron_emiso03[njetmx], Electron_hadiso03[njetmx], Electron_emiso04[njetmx], Electron_hadiso04[njetmx];
  float Electron_eoverp[njetmx], Electron_ietaieta[njetmx], Electron_etain[njetmx], Electron_phiin[njetmx], Electron_fbrem[njetmx];
  float Electron_nohits[njetmx], Electron_misshits[njetmx];
  float Electron_pfiso_drcor[njetmx];
  float Electron_pfiso_eacor[njetmx];
  float Electron_pfiso04_eacor[njetmx];
  int Electron_ndf[njetmx];
  float Electron_eccalTrkEnergyPostCorr[njetmx];
  float Electron_energyScaleValue[njetmx];
  float Electron_energyScaleUp[njetmx];
  float Electron_energyScaleDown[njetmx];
  float Electron_energySigmaValue[njetmx];
  float Electron_energySigmaUp[njetmx];
  float Electron_energySigmaDown[njetmx];
  float Electron_supcl_eta[njetmx];
  float Electron_supcl_phi[njetmx];
  float Electron_supcl_e[njetmx];
  float Electron_supcl_rawE[njetmx];
  float Electron_sigmaieta[njetmx];
  float Electron_sigmaiphi[njetmx];
  float Electron_r9full[njetmx];
  float Electron_supcl_etaw[njetmx];
  float Electron_supcl_phiw[njetmx];
  float Electron_hcaloverecal[njetmx];
  float Electron_cloctftrkn[njetmx];
  float Electron_cloctftrkchi2[njetmx];
  float Electron_e1x5bye5x5[njetmx];
  float Electron_normchi2[njetmx];
  float Electron_hitsmiss[njetmx];
  float Electron_trkmeasure[njetmx];
  float Electron_convtxprob[njetmx];
  float Electron_ecloverpout[njetmx];
  float Electron_ecaletrkmomentum[njetmx];
  float Electron_deltaetacltrkcalo[njetmx];
  float Electron_supcl_preshvsrawe[njetmx];
  bool Electron_convVeto[njetmx];
  float Electron_pfisolsumphet[njetmx];
  float Electron_pfisolsumchhadpt[njetmx];
  float Electron_pfsiolsumneuhadet[njetmx];
  float Electron_minchiso[njetmx];
  float Electron_minnhiso[njetmx];
  float Electron_minphiso[njetmx];
  float Electron_minisoall[njetmx];

  int nPhoton;
  bool Photon_mvaid_Fall17V2_WP90[njetmx];
  bool Photon_mvaid_Fall17V2_WP80[njetmx];
  float Photon_mvaid_Fall17V2_raw[njetmx];
  bool Photon_mvaid_Spring16V1_WP90[njetmx];
  bool Photon_mvaid_Spring16V1_WP80[njetmx];
  float Photon_e[njetmx];
  float Photon_pt[njetmx];
  float Photon_eta[njetmx];
  float Photon_phi[njetmx];
  float Photon_e1by9[njetmx];
  float Photon_e9by25[njetmx];
  float Photon_hadbyem[njetmx];
  float Photon_trkiso[njetmx];
  float Photon_emiso[njetmx];
  float Photon_hadiso[njetmx];
  float Photon_chhadiso[njetmx];
  float Photon_neuhadiso[njetmx];
  float Photon_PUiso[njetmx];
  float Photon_phoiso[njetmx];
  float Photon_ietaieta[njetmx];

  int nGenJetAK4;
  float GenJetAK4_pt[njetmx];
  float GenJetAK4_eta[njetmx];
  float GenJetAK4_phi[njetmx];
  float GenJetAK4_mass[njetmx];
  int GenJetAK4_hadronflav[njetmx];

  float Generator_qscale, Generator_x1, Generator_x2, Generator_xpdf1, Generator_xpdf2, Generator_scalePDF;
  int Generator_id1, Generator_id2;

  int npu_vert;
  int npu_vert_true;

  bool hlt_IsoMu24;
  bool hlt_Ele32_WPTight_Gsf;
  bool hlt_Mu17_Mu8;
  bool hlt_Ele23_Ele12;
  bool hlt_DoubleEle33;
  bool hlt_DoubleEle25_CaloIdL_MW;

  int trig_value;

  float lep1pt, lep1y, lep1eta, lep1phi, lep1e, lep2pt, lep2y, lep2eta, lep2phi, lep2e, dilep_invmass;

  int b1_hadFlv, b2_hadFlv, b1_partFlv, b2_partFlv;
  float bjet_no, b1pt, b1y, b1eta, b1phi, b1e, b2pt, b2y, b2eta, b2phi, b2e, bb_inv_mass, b1_DeepFlv, b2_DeepFlv, b1_PUid, b2_PUid;
  float b1pt_jesup, b1y_jesup, b1eta_jesup, b1phi_jesup, b1e_jesup, b2pt_jesup, b2y_jesup, b2eta_jesup, b2phi_jesup, b2e_jesup, bb_inv_mass_jesup;
  float b1pt_jesdn, b1y_jesdn, b1eta_jesdn, b1phi_jesdn, b1e_jesdn, b2pt_jesdn, b2y_jesdn, b2eta_jesdn, b2phi_jesdn, b2e_jesdn, bb_inv_mass_jesdn;
  float b1pt_resoup, b1y_resoup, b1eta_resoup, b1phi_resoup, b1e_resoup, b2pt_resoup, b2y_resoup, b2eta_resoup, b2phi_resoup, b2e_resoup, bb_inv_mass_resoup;
  float b1pt_resodn, b1y_resodn, b1eta_resodn, b1phi_resodn, b1e_resodn, b2pt_resodn, b2y_resodn, b2eta_resodn, b2phi_resodn, b2e_resodn, bb_inv_mass_resodn;

  float pho_no, pho1pt, pho1y, pho1eta, pho1phi, pho1e, pho2pt, pho2y, pho2eta, pho2phi, pho2e, dipho_invmass;

  float invmassbbgg, invmassbbgg_jesup, invmassbbgg_jesdn, invmassbbgg_resoup, invmassbbgg_resodn;
  float bb_gg_dphi, bb_gg_dphi_jesup, bb_gg_dphi_jesdn, bb_gg_dphi_resoup, bb_gg_dphi_resodn;

  double weight_nom, weight_PU_up, weight_leptonsf_up, weight_prefiring_up, weight_bSF_up, weight_trig_up, weight_PU_dn, weight_leptonsf_dn, weight_prefiring_dn, weight_bSF_dn, weight_trig_dn;

  float Jet_pt_nom[njetmx], Jet_mass_nom[njetmx], Jet_pt_jesup[njetmx], Jet_mass_jesup[njetmx], Jet_pt_jesdn[njetmx], Jet_mass_jesdn[njetmx], Jet_pt_resoup[njetmx], Jet_mass_resoup[njetmx], Jet_pt_resodn[njetmx], Jet_mass_resodn[njetmx];

  int aa, ab, bb;
  bool isEle, isMu;
  bool no_CR_1, no_CR_2;


  int year = 2018;
  string muon_id_name = "Tight";
  string electron_id_name = "wp90noiso";

  float DAK4_T = 0.71;
  float DAK4_M = 0.2783;
  float DAK4_L = 0.0490;



int main(int argc, char *argv[])
{

	isMC = true;
        isFastSIM = false;
        char fOut[50];
        string inputFile=argv[3];
        string path="/home/abala/cms/CMSSW_10_5_0/src/condor_job/";


        calib_deepflav = BTagCalibration("DeepJet", "BtagRecommendation106XUL18/DeepJet_106XUL18SF_WPonly_V1p1.csv");
        reader_deepflav = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"});
        reader_deepflav.load(calib_deepflav, BTagEntry::FLAV_B, "comb");
        reader_deepflav.load(calib_deepflav, BTagEntry::FLAV_C, "comb");
        reader_deepflav.load(calib_deepflav, BTagEntry::FLAV_UDSG, "incl");

	char name[1000];

        TFile *file_mu_sf;
        sprintf(name,"data/Efficiencies_muon_generalTracks_Z_Run%i_UL_ID.root",year);
        file_mu_sf = new TFile(name,"read");

        TFile *file_el_sf;
        sprintf(name,"data/egammaEffi.txt_Ele_%s_EGM2D_UL%i.root",electron_id_name.c_str(),year);
        file_el_sf = new TFile(name,"read");

        TFile *file_pu_ratio;
        sprintf(name,"data/pileup/RatioPileup-UL%i-100bins.root",year);
        file_pu_ratio = new TFile(name,"read");


	if(inputFile=="new_v1/ZH_mA_20_test.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/FastSIM/ZH_mA_20/ZH_mA_20_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="new_v1/ZH_mA_55_test.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/FastSIM/ZH_mA_55/ZH_mA_55_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="new_v1/TTGJets_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTGJets/2018/ZH/TTGJets_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="new_v1/TTGG_0Jets_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTGG/2018/ZH/TTGG_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="new_v1/DYJetsToLL_M50_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/DYJetsToLL/2018/ZH/DYJetsToLL_M50_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="new_v1/TTTo2L2Nu_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTTo2L2Nu/2018/ZH/TTTo2L2Nu_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="new_v1/TTToSemiLeptonic_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTSL/2018/ZH/TTSL_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="new_v1/WJetsToLNu_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/WJetsToLNu/2018/ZH/WJetsToLNu_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="new_v1/EGamma_2018A_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/egamma/egamma_2018A_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="new_v1/EGamma_2018B_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/egamma/egamma_2018B_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="new_v1/EGamma_2018C_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/egamma/egamma_2018C_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="new_v1/EGamma_2018D_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/egamma/egamma_2018D_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="new_v1/DoubleMuon_2018A_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/DoubMu/DoubMu_2018A_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="new_v1/DoubleMuon_2018B_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/DoubMu/DoubMu_2018B_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="new_v1/DoubleMuon_2018C_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/DoubMu/DoubMu_2018C_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="new_v1/DoubleMuon_2018D_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/DoubMu/DoubMu_2018D_%s_%s.root",argv[1],argv[2]);
        }

        else{
                cout<<"Input file does not exist"<<endl;
                exit(0);
        }


	TFile *fSF_dilepton;
        TH2F  *h_trigmumu; TH2F  *h_trigemu; TH2F  *h_trigee;

	if(!isDATA)
	{
		fSF_dilepton = TFile::Open("Trigger_SF/DileptonTriggerSF_UL_miniAODv2/TriggerSF_2018_ULv2.root");
		h_trigmumu = (TH2F*)fSF_dilepton->Get("h2D_SF_mumu_lepABpt_FullError");
        	h_trigee   = (TH2F*)fSF_dilepton->Get("h2D_SF_ee_lepABpt_FullError");
        	h_trigemu  = (TH2F*)fSF_dilepton->Get("h2D_SF_emu_lepABpt_FullError");
	}

        TFile *fout = new TFile(fOut,"RECREATE");
        TTree *Tout = new TTree("Tout", "Info");
        
        
	Tout->Branch("lep1pt", &lep1pt, "lep1pt/F");
        Tout->Branch("lep1y", &lep1y, "lep1y/F");
        Tout->Branch("lep1eta", &lep1eta, "lep1eta/F");
        Tout->Branch("lep1phi", &lep1phi, "lep1phi/F");
        Tout->Branch("lep1e", &lep1e, "lep1e/F");
	Tout->Branch("lep2pt", &lep2pt, "lep2pt/F");
        Tout->Branch("lep2y", &lep2y, "lep2y/F");
        Tout->Branch("lep2eta", &lep2eta, "lep2eta/F");
        Tout->Branch("lep2phi", &lep2phi, "lep2phi/F");
        Tout->Branch("lep2e", &lep2e, "lep2e/F");
	Tout->Branch("dilep_invmass", &dilep_invmass, "dilep_invmass/F");

	Tout->Branch("bjet_no", &bjet_no, "bjet_no/I");
        Tout->Branch("b1pt", &b1pt, "b1pt/F");
        Tout->Branch("b1y", &b1y, "b1y/F");
        Tout->Branch("b1eta", &b1eta, "b1eta/F");
        Tout->Branch("b1phi", &b1phi, "b1phi/F");
        Tout->Branch("b1e", &b1e, "b1e/F");
        Tout->Branch("b2pt", &b2pt, "b2pt/F");
        Tout->Branch("b2y", &b2y, "b2y/F");
        Tout->Branch("b2eta", &b2eta, "b2eta/F");
        Tout->Branch("b2phi", &b2phi, "b2phi/F");
        Tout->Branch("b2e", &b2e, "b2e/F");
        Tout->Branch("bb_inv_mass", &bb_inv_mass, "bb_inv_mass/F");
        Tout->Branch("b1_DeepFlv", &b1_DeepFlv, "b1_DeepFlv/F");
        Tout->Branch("b2_DeepFlv", &b2_DeepFlv, "b2_DeepFlv/F");
        Tout->Branch("b1_hadFlv", &b1_hadFlv, "b1_hadFlv/I");
        Tout->Branch("b2_hadFlv", &b2_hadFlv, "b2_hadFlv/I");
        Tout->Branch("b1_partFlv", &b1_partFlv, "b1_partFlv/I");
        Tout->Branch("b2_partFlv", &b2_partFlv, "b2_partFlv/I");
        Tout->Branch("b1_PUid", &b1_PUid, "b1_PUid/F");
        Tout->Branch("b2_PUid", &b2_PUid, "b2_PUid/F");

        Tout->Branch("b1pt_jesup", &b1pt_jesup, "b1pt_jesup/F");
        Tout->Branch("b1y_jesup", &b1y_jesup, "b1y_jesup/F");
        Tout->Branch("b1eta_jesup", &b1eta_jesup, "b1eta_jesup/F");
        Tout->Branch("b1phi_jesup", &b1phi_jesup, "b1phi_jesup/F");
        Tout->Branch("b1e_jesup", &b1e_jesup, "b1e_jesup/F");
        Tout->Branch("b2pt_jesup", &b2pt_jesup, "b2pt_jesup/F");
        Tout->Branch("b2y_jesup", &b2y_jesup, "b2y_jesup/F");
        Tout->Branch("b2eta_jesup", &b2eta_jesup, "b2eta_jesup/F");
        Tout->Branch("b2phi_jesup", &b2phi_jesup, "b2phi_jesup/F");
        Tout->Branch("b2e_jesup", &b2e_jesup, "b2e_jesup/F");
        Tout->Branch("bb_inv_mass_jesup", &bb_inv_mass_jesup, "bb_inv_mass_jesup/F");

	Tout->Branch("b1pt_jesdn", &b1pt_jesdn, "b1pt_jesdn/F");
        Tout->Branch("b1y_jesdn", &b1y_jesdn, "b1y_jesdn/F");
        Tout->Branch("b1eta_jesdn", &b1eta_jesdn, "b1eta_jesdn/F");
        Tout->Branch("b1phi_jesdn", &b1phi_jesdn, "b1phi_jesdn/F");
        Tout->Branch("b1e_jesdn", &b1e_jesdn, "b1e_jesdn/F");
        Tout->Branch("b2pt_jesdn", &b2pt_jesdn, "b2pt_jesdn/F");
        Tout->Branch("b2y_jesdn", &b2y_jesdn, "b2y_jesdn/F");
        Tout->Branch("b2eta_jesdn", &b2eta_jesdn, "b2eta_jesdn/F");
        Tout->Branch("b2phi_jesdn", &b2phi_jesdn, "b2phi_jesdn/F");
        Tout->Branch("b2e_jesdn", &b2e_jesdn, "b2e_jesdn/F");
        Tout->Branch("bb_inv_mass_jesdn", &bb_inv_mass_jesdn, "bb_inv_mass_jesdn/F");

        Tout->Branch("b1pt_resoup", &b1pt_resoup, "b1pt_resoup/F");
        Tout->Branch("b1y_resoup", &b1y_resoup, "b1y_resoup/F");
        Tout->Branch("b1eta_resoup", &b1eta_resoup, "b1eta_resoup/F");
        Tout->Branch("b1phi_resoup", &b1phi_resoup, "b1phi_resoup/F");
        Tout->Branch("b1e_resoup", &b1e_resoup, "b1e_resoup/F");
        Tout->Branch("b2pt_resoup", &b2pt_resoup, "b2pt_resoup/F");
        Tout->Branch("b2y_resoup", &b2y_resoup, "b2y_resoup/F");
        Tout->Branch("b2eta_resoup", &b2eta_resoup, "b2eta_resoup/F");
        Tout->Branch("b2phi_resoup", &b2phi_resoup, "b2phi_resoup/F");
        Tout->Branch("b2e_resoup", &b2e_resoup, "b2e_resoup/F");
        Tout->Branch("bb_inv_mass_resoup", &bb_inv_mass_resoup, "bb_inv_mass_resoup/F");

        Tout->Branch("b1pt_resodn", &b1pt_resodn, "b1pt_resodn/F");
        Tout->Branch("b1y_resodn", &b1y_resodn, "b1y_resodn/F");
        Tout->Branch("b1eta_resodn", &b1eta_resodn, "b1eta_resodn/F");
        Tout->Branch("b1phi_resodn", &b1phi_resodn, "b1phi_resodn/F");
        Tout->Branch("b1e_resodn", &b1e_resodn, "b1e_resodn/F");
        Tout->Branch("b2pt_resodn", &b2pt_resodn, "b2pt_resodn/F");
        Tout->Branch("b2y_resodn", &b2y_resodn, "b2y_resodn/F");
        Tout->Branch("b2eta_resodn", &b2eta_resodn, "b2eta_resodn/F");
        Tout->Branch("b2phi_resodn", &b2phi_resodn, "b2phi_resodn/F");
        Tout->Branch("b2e_resodn", &b2e_resodn, "b2e_resodn/F");
        Tout->Branch("bb_inv_mass_resodn", &bb_inv_mass_resodn, "bb_inv_mass_resodn/F");

	Tout->Branch("pho_no", &pho_no, "pho_no/I");
        Tout->Branch("pho1pt", &pho1pt, "pho1pt/F");
        Tout->Branch("pho2pt", &pho2pt, "pho2pt/F");
        Tout->Branch("pho1y", &pho1y, "pho1y/F");
        Tout->Branch("pho2y", &pho2y, "pho2y/F");
        Tout->Branch("pho1eta", &pho1eta, "pho1eta/F");
        Tout->Branch("pho1phi", &pho1phi, "pho1phi/F");
        Tout->Branch("pho1e", &pho1e, "pho1e/F");
        Tout->Branch("pho2eta", &pho2eta, "pho2eta/F");
        Tout->Branch("pho2phi", &pho2phi, "pho2phi/F");
        Tout->Branch("pho2e", &pho2e, "pho2e/F");
        Tout->Branch("dipho_invmass", &dipho_invmass, "dipho_invmass/F");

        Tout->Branch("invmassbbgg", &invmassbbgg, "invmassbbgg/F");
        Tout->Branch("invmassbbgg_jesup", &invmassbbgg_jesup, "invmassbbgg_jesup/F");
        Tout->Branch("invmassbbgg_jesdn", &invmassbbgg_jesdn, "invmassbbgg_jesdn/F");
        Tout->Branch("invmassbbgg_resoup", &invmassbbgg_resoup, "invmassbbgg_resoup/F");
        Tout->Branch("invmassbbgg_resodn", &invmassbbgg_resodn, "invmassbbgg_resodn/F");

        Tout->Branch("bb_gg_dphi", &bb_gg_dphi, "bb_gg_dphi/F");
        Tout->Branch("bb_gg_dphi_jesup", &bb_gg_dphi_jesup, "bb_gg_dphi_jesup/F");
        Tout->Branch("bb_gg_dphi_jesdn", &bb_gg_dphi_jesdn, "bb_gg_dphi_jesdn/F");
        Tout->Branch("bb_gg_dphi_resoup", &bb_gg_dphi_resoup, "bb_gg_dphi_resoup/F");
        Tout->Branch("bb_gg_dphi_resodn", &bb_gg_dphi_resodn, "bb_gg_dphi_resodn/F");

	Tout->Branch("weight_nom", &weight_nom, "weight_nom/D");
        Tout->Branch("weight_PU_up", &weight_PU_up, "weight_PU_up/D");
        Tout->Branch("weight_leptonsf_up", &weight_leptonsf_up, "weight_leptonsf_up/D");
	Tout->Branch("weight_bSF_up", &weight_bSF_up, "weight_bSF_up/D");
        Tout->Branch("weight_trig_up", &weight_trig_up, "weight_trig_up/D");
        Tout->Branch("weight_PU_dn", &weight_PU_dn, "weight_PU_dn/D");
        Tout->Branch("weight_leptonsf_dn", &weight_leptonsf_dn, "weight_leptonsf_dn/D");
	Tout->Branch("weight_bSF_dn", &weight_bSF_dn, "weight_bSF_dn/D");
        Tout->Branch("weight_trig_dn", &weight_trig_dn, "weight_trig_dn/D");

	Tout->Branch("ab", &ab, "ab/I");
        Tout->Branch("bb", &bb, "bb/I");

	Tout->Branch("isEle", &isEle, "isEle/O");
        Tout->Branch("isMu", &isMu, "isMu/O");

      
        TH1F *cut_flow_el = new TH1F("cut_flow_el","cut_flow_el",20,0,20);
        cut_flow_el->Sumw2();

        TH1F *cut_flow_mu = new TH1F("cut_flow_mu","cut_flow_mu",20,0,20);
        cut_flow_mu->Sumw2();
        

	Double_t AK4PtBINS[] =  {20, 30, 50, 70, 100, 140, 200, 300, 600, 1000, 7000};
        const int nAK4PtBINS = sizeof(AK4PtBINS)/sizeof(AK4PtBINS[0])-1;
        Double_t EtaBINS[] = {0,0.6,1.2,2.5};
        const int nEtaBINS = sizeof(EtaBINS)/sizeof(EtaBINS[0])-1;

        TH2F* h_Ak4_b_flv = new TH2F("h_Ak4_b_flv", "h_Ak4_b_flv", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_c_flv = new TH2F("h_Ak4_c_flv", "h_Ak4_c_flv", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_l_flv = new TH2F("h_Ak4_l_flv", "h_Ak4_l_flv", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);

        TH2F* h_Ak4_b_flv_pass_L = new TH2F("h_Ak4_b_flv_pass_L", "h_Ak4_b_flv_pass_L", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_c_flv_pass_L = new TH2F("h_Ak4_c_flv_pass_L", "h_Ak4_c_flv_pass_L", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_l_flv_pass_L = new TH2F("h_Ak4_l_flv_pass_L", "h_Ak4_l_flv_pass_L", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);

        TH2F* h_Ak4_b_flv_pass_M = new TH2F("h_Ak4_b_flv_pass_M", "h_Ak4_b_flv_pass_M", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_c_flv_pass_M = new TH2F("h_Ak4_c_flv_pass_M", "h_Ak4_c_flv_pass_M", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_l_flv_pass_M = new TH2F("h_Ak4_l_flv_pass_M", "h_Ak4_l_flv_pass_M", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);

        TH2F* h_Ak4_b_flv_pass_T = new TH2F("h_Ak4_b_flv_pass_T", "h_Ak4_b_flv_pass_T", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_c_flv_pass_T = new TH2F("h_Ak4_c_flv_pass_T", "h_Ak4_c_flv_pass_T", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_l_flv_pass_T = new TH2F("h_Ak4_l_flv_pass_T", "h_Ak4_l_flv_pass_T", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);

  
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


   TTree *T1 = (TTree*)f->Get("Events");

   T1->SetBranchAddress("trig_value",&trig_value);
   T1->SetBranchAddress("hlt_IsoMu24",&hlt_IsoMu24);
   T1->SetBranchAddress("hlt_Ele32_WPTight_Gsf",&hlt_Ele32_WPTight_Gsf);
   T1->SetBranchAddress("hlt_Mu17_Mu8",&hlt_Mu17_Mu8);
   T1->SetBranchAddress("hlt_Ele23_Ele12",&hlt_Ele23_Ele12);
   T1->SetBranchAddress("hlt_DoubleEle33",&hlt_DoubleEle33);
   T1->SetBranchAddress("hlt_DoubleEle25_CaloIdL_MW",&hlt_DoubleEle25_CaloIdL_MW);

   T1->SetBranchAddress("prefiringweight",&prefiringweight);
   T1->SetBranchAddress("prefiringweightup",&prefiringweightup);
   T1->SetBranchAddress("prefiringweightdown",&prefiringweightdown);

   T1->SetBranchAddress("nPFJetAK4",&nPFJetAK4);
   T1->SetBranchAddress("PFJetAK4_pt",PFJetAK4_pt);
   T1->SetBranchAddress("PFJetAK4_eta",PFJetAK4_eta);
   T1->SetBranchAddress("PFJetAK4_y",PFJetAK4_y);
   T1->SetBranchAddress("PFJetAK4_phi",PFJetAK4_phi);
   T1->SetBranchAddress("PFJetAK4_mass",PFJetAK4_mass);
   T1->SetBranchAddress("PFJetAK4_jetID",PFJetAK4_jetID);
   T1->SetBranchAddress("PFJetAK4_jetID_tightlepveto",PFJetAK4_jetID_tightlepveto);
   T1->SetBranchAddress("PFJetAK4_JEC",PFJetAK4_JEC);
   T1->SetBranchAddress("PFJetAK4_btag_DeepCSV",PFJetAK4_btag_DeepCSV);
   T1->SetBranchAddress("PFJetAK4_btag_DeepFlav",PFJetAK4_btag_DeepFlav);
   T1->SetBranchAddress("PFJetAK4_JER",PFJetAK4_reso);
   T1->SetBranchAddress("PFJetAK4_JERup",PFJetAK4_resoup);
   T1->SetBranchAddress("PFJetAK4_JERdn",PFJetAK4_resodn);
   T1->SetBranchAddress("PFJetAK4_jesup_Total",PFJetAK4_jesup_Total);
   T1->SetBranchAddress("PFJetAK4_jesdn_Total",PFJetAK4_jesdn_Total);
   T1->SetBranchAddress("PFJetAK4_btag_DeepCSV_SF",PFJetAK4_btag_DeepCSV_SF);
   T1->SetBranchAddress("PFJetAK4_btag_DeepCSV_SF_up",PFJetAK4_btag_DeepCSV_SF_up);
   T1->SetBranchAddress("PFJetAK4_btag_DeepCSV_SF_dn",PFJetAK4_btag_DeepCSV_SF_dn);
   T1->SetBranchAddress("PFJetAK4_btag_DeepFlav_SF",PFJetAK4_btag_DeepFlav_SF);
   T1->SetBranchAddress("PFJetAK4_btag_DeepFlav_SF_up",PFJetAK4_btag_DeepFlav_SF_up);
   T1->SetBranchAddress("PFJetAK4_btag_DeepFlav_SF_dn",PFJetAK4_btag_DeepFlav_SF_dn);
   T1->SetBranchAddress("PFJetAK4_hadronflav",PFJetAK4_hadronflav);
   T1->SetBranchAddress("PFJetAK4_partonflav",PFJetAK4_partonflav);
   T1->SetBranchAddress("PFJetAK4_PUID",PFJetAK4_PUID);

   T1->SetBranchAddress("nMuon",&nMuon);
  T1->SetBranchAddress("Muon_isPF",Muon_isPF);
  T1->SetBranchAddress("Muon_isGL",Muon_isGL);
  T1->SetBranchAddress("Muon_isTRK",Muon_isTRK);
  T1->SetBranchAddress("Muon_isLoose",Muon_isLoose);
  T1->SetBranchAddress("Muon_isGoodGL",Muon_isGoodGL);
  T1->SetBranchAddress("Muon_isMed",Muon_isMed);
  T1->SetBranchAddress("Muon_isMedPr",Muon_isMedPr);
  T1->SetBranchAddress("Muon_isTight",Muon_isTight);
  T1->SetBranchAddress("Muon_isHighPt",Muon_isHighPt);
  T1->SetBranchAddress("Muon_isHighPttrk",Muon_isHighPttrk);
  T1->SetBranchAddress("Muon_TightID",Muon_TightID);
  T1->SetBranchAddress("Muon_pt",Muon_pt);
  T1->SetBranchAddress("Muon_p",Muon_p);
  T1->SetBranchAddress("Muon_eta",Muon_eta);
  T1->SetBranchAddress("Muon_phi",Muon_phi);
//  T1->SetBranchAddress("Muon_e",Muon_e);
  T1->SetBranchAddress("Muon_minisoch", Muon_minchiso);
  T1->SetBranchAddress("Muon_minisonh", Muon_minnhiso);
  T1->SetBranchAddress("Muon_minisoph", Muon_minphiso);
  T1->SetBranchAddress("Muon_minisoall", Muon_minisoall);
  T1->SetBranchAddress("Muon_dxy",Muon_dxy);
  T1->SetBranchAddress("Muon_dz",Muon_dz);
  T1->SetBranchAddress("Muon_dxyErr",Muon_dxyErr);
  T1->SetBranchAddress("Muon_ip3d",Muon_ip3d);
  T1->SetBranchAddress("Muon_ptErr",Muon_ptErr);
  T1->SetBranchAddress("Muon_chi",Muon_chi);
  T1->SetBranchAddress("Muon_ndf",Muon_ndf);
  T1->SetBranchAddress("Muon_ecal",Muon_ecal);
  T1->SetBranchAddress("Muon_hcal",Muon_hcal);
  T1->SetBranchAddress("Muon_pfiso",Muon_pfiso);
  T1->SetBranchAddress("Muon_posmatch",Muon_posmatch);
  T1->SetBranchAddress("Muon_trkink",Muon_trkink);
  T1->SetBranchAddress("Muon_segcom",Muon_segcom);
  T1->SetBranchAddress("Muon_hit",Muon_hit);
  T1->SetBranchAddress("Muon_pixhit",Muon_pixhit);
  T1->SetBranchAddress("Muon_mst",Muon_mst);
  T1->SetBranchAddress("Muon_trklay",Muon_trklay);
  T1->SetBranchAddress("Muon_valfrac",Muon_valfrac);
  T1->SetBranchAddress("Muon_dxy_sv",Muon_dxy_sv);
  T1->SetBranchAddress("Muon_corrected_pt",Muon_corrected_pt);
  T1->SetBranchAddress("Muon_correctedUp_pt",Muon_correctedUp_pt);
  T1->SetBranchAddress("Muon_correctedDown_pt",Muon_correctedDown_pt);

  T1->SetBranchAddress("nElectron",&nElectron);
  T1->SetBranchAddress("Electron_pt",Electron_pt);
  T1->SetBranchAddress("Electron_eta",Electron_eta);
  T1->SetBranchAddress("Electron_phi",Electron_phi);
  T1->SetBranchAddress("Electron_p",Electron_p);
  T1->SetBranchAddress("Electron_e",Electron_e);
  T1->SetBranchAddress("Electron_e_ECAL",Electron_e_ECAL);
  T1->SetBranchAddress("Electron_mvaid_Fallv2WP90",Electron_mvaid_Fallv2WP90);
  T1->SetBranchAddress("Electron_mvaid_Fallv2WP90_noIso",Electron_mvaid_Fallv2WP90_noIso);
  T1->SetBranchAddress("Electron_mvaid_Fallv2WP80",Electron_mvaid_Fallv2WP80);
  T1->SetBranchAddress("Electron_mvaid_Fallv2WP80_noIso",Electron_mvaid_Fallv2WP80_noIso);
  T1->SetBranchAddress("Electron_dxy",Electron_dxy);
  T1->SetBranchAddress("Electron_dxyErr",Electron_dxyErr);
  T1->SetBranchAddress("Electron_dz",Electron_dz);
  T1->SetBranchAddress("Electron_dzErr",Electron_dzErr);
  T1->SetBranchAddress("Electron_ip3d",Electron_ip3d);
  T1->SetBranchAddress("Electron_dxy_sv",Electron_dxy_sv);
  T1->SetBranchAddress("Electron_hovere",Electron_hovere);
  T1->SetBranchAddress("Electron_chi",Electron_chi);
  T1->SetBranchAddress("Electron_ndf",Electron_ndf);
  T1->SetBranchAddress("Electron_eoverp",Electron_eoverp);
  T1->SetBranchAddress("Electron_ietaieta",Electron_ietaieta);
  T1->SetBranchAddress("Electron_misshits",Electron_misshits);
  T1->SetBranchAddress("Electron_pfiso_drcor",Electron_pfiso_drcor);
  T1->SetBranchAddress("Electron_pfiso_eacor",Electron_pfiso_eacor);
  T1->SetBranchAddress("Electron_pfiso04_eacor",Electron_pfiso04_eacor);
  T1->SetBranchAddress("Electron_r9full", Electron_r9full);
  T1->SetBranchAddress("Electron_hcaloverecal", Electron_hcaloverecal);
  T1->SetBranchAddress("Electron_hitsmiss", Electron_hitsmiss);
  T1->SetBranchAddress("Electron_ecloverpout", Electron_ecloverpout);
  T1->SetBranchAddress("Electron_convVeto", Electron_convVeto);

  T1->SetBranchAddress("nPhoton",&nPhoton);
  T1->SetBranchAddress("Photon_e",Photon_e);
  T1->SetBranchAddress("Photon_pt",Photon_pt);
  T1->SetBranchAddress("Photon_eta",Photon_eta);
  T1->SetBranchAddress("Photon_phi",Photon_phi);
  T1->SetBranchAddress("Photon_mvaid_Fall17V2_raw",Photon_mvaid_Fall17V2_raw);
  T1->SetBranchAddress("Photon_mvaid_Fall17V2_WP90",Photon_mvaid_Fall17V2_WP90);
  T1->SetBranchAddress("Photon_mvaid_Fall17V2_WP80",Photon_mvaid_Fall17V2_WP80);
  T1->SetBranchAddress("Photon_mvaid_Spring16V1_WP90",Photon_mvaid_Spring16V1_WP90);
  T1->SetBranchAddress("Photon_mvaid_Spring16V1_WP80",Photon_mvaid_Spring16V1_WP80);
  T1->SetBranchAddress("Photon_e1by9",Photon_e1by9);
  T1->SetBranchAddress("Photon_e9by25",Photon_e9by25);
  T1->SetBranchAddress("Photon_trkiso",Photon_trkiso);
  T1->SetBranchAddress("Photon_emiso",Photon_emiso);
  T1->SetBranchAddress("Photon_hadiso",Photon_hadiso);
  T1->SetBranchAddress("Photon_chhadiso",Photon_chhadiso);
  T1->SetBranchAddress("Photon_neuhadiso",Photon_neuhadiso);
  T1->SetBranchAddress("Photon_phoiso",Photon_phoiso);
  T1->SetBranchAddress("Photon_PUiso",Photon_PUiso);
  T1->SetBranchAddress("Photon_hadbyem",Photon_hadbyem);
  T1->SetBranchAddress("Photon_ietaieta",Photon_ietaieta);

  if(isMC){

  // generator-related info //

  T1->SetBranchAddress("Generator_weight", &Generator_weight);
  T1->SetBranchAddress("Generator_qscale",&Generator_qscale);
  T1->SetBranchAddress("Generator_x1",&Generator_x1);
  T1->SetBranchAddress("Generator_x2",&Generator_x2);
  T1->SetBranchAddress("Generator_xpdf1",&Generator_xpdf1);
  T1->SetBranchAddress("Generator_xpdf2",&Generator_xpdf2);
  T1->SetBranchAddress("Generator_id1",&Generator_id1);
  T1->SetBranchAddress("Generator_id2",&Generator_id2);
  T1->SetBranchAddress("Generator_scalePDF",&Generator_scalePDF);

  T1->SetBranchAddress("nGenJetAK4", &nGenJetAK4);
  T1->SetBranchAddress("GenJetAK4_pt",GenJetAK4_pt);
  T1->SetBranchAddress("GenJetAK4_eta",GenJetAK4_eta);
  T1->SetBranchAddress("GenJetAK4_phi",GenJetAK4_phi);
  T1->SetBranchAddress("GenJetAK4_mass",GenJetAK4_mass);
  T1->SetBranchAddress("GenJetAK4_hadronflav",GenJetAK4_hadronflav);

  T1->SetBranchAddress("npu_vert",&npu_vert);
  T1->SetBranchAddress("npu_vert_true",&npu_vert_true);

  T1->SetBranchAddress("LHE_weight",&LHE_weight);
  }
  
    
    
    int nevents;
    nevents=T1->GetEntries();
    
    for(int iev=0; iev<nevents; iev++)
    {
            
            T1->GetEntry(iev);
   
   	    cut_flow_el->Fill(0);
            cut_flow_mu->Fill(0);


        //-----------------------------------------------------------------------------event weight info--------------------------------------------------------------------------------//


          double event_weight;

              if(isMC){

                if (isFastSIM)
                {
                        event_weight = 1.0;
                }

                if(fabs(LHE_weight)>1.e-12) { event_weight = LHE_weight; }
                else { event_weight = Generator_weight; }
                      }
              else{
                event_weight = 1.;
                  }	    


	      //-----------------------------------------------------------------------AK4 histogram filling for btag SF-----------------------------------------------------------------//


              if(isMC){
                for(unsigned ijet=0; ijet<nPFJetAK4; ijet++){
                        double dR = 9999.9;
                        for(unsigned gjet=0; gjet<nGenJetAK4; gjet++)
                        {
                                double temp_dR = delta2R(PFJetAK4_eta[ijet],PFJetAK4_phi[ijet],GenJetAK4_eta[gjet],GenJetAK4_phi[gjet]) ;
                                if (temp_dR < dR )
                                {
                                        dR = temp_dR;
                                }
                        }
                        if(dR < 0.4)
                        {
                                if( abs(PFJetAK4_hadronflav[ijet]) == 5 )  {
                                        h_Ak4_b_flv->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight);
                                        if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_T  )  { h_Ak4_b_flv_pass_T->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
                                        if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_M )   { h_Ak4_b_flv_pass_M->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
                                        if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_L )   { h_Ak4_b_flv_pass_L->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
                                }
                                else if( abs(PFJetAK4_hadronflav[ijet]) == 4 )  {
                                        h_Ak4_c_flv->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight);
                                        if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_T  )  { h_Ak4_c_flv_pass_T->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
                                        if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_M )   { h_Ak4_c_flv_pass_M->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
                                        if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_L )   { h_Ak4_c_flv_pass_L->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
                                }
                                else if( abs(PFJetAK4_hadronflav[ijet]) == 0 )  {
                                        h_Ak4_l_flv->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight);
                                        if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_T  )  { h_Ak4_l_flv_pass_T->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
                                        if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_M )   { h_Ak4_l_flv_pass_M->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
                                        if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_L )   { h_Ak4_l_flv_pass_L->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
                                }
                        }
                }
                }


	  //--------------------------------------------------------------------------trigger info---------------------------------------------------------------------//


          if (!(hlt_Mu17_Mu8 || hlt_Ele23_Ele12)) continue;


          //---------------------------------------------------------------------------lepton info---------------------------------------------------------------------//



	  int nele = 0;
          int nmu = 0;

	  TLorentzVector lepton1vec, lepton2vec;


	  isEle = false;
            
	  if (hlt_Ele23_Ele12) {
          cut_flow_el->Fill(1);

          for(int i=0; i<nElectron; i++)
          {
                  if (Electron_hitsmiss[i]>1) continue;
                  if (fabs(Electron_dxy[i])>0.045) continue;
                  if (Electron_pfiso_drcor[i]>=0.15) continue;
                  if(fabs(Electron_eta[i]) > 2.5) continue;                                    //electron channel

		  Electron_pt[nele] = Electron_pt[i];
                  Electron_eta[nele] = Electron_eta[i];
                  Electron_phi[nele] = Electron_phi[i];
                  Electron_e[nele] = Electron_e[i];

                  if (++nele >= njetmx) break;
          }

	  nElectron = nele;
	  if (nElectron < 2) continue;

//	  if (nElectron >= 2) {
	  cut_flow_el->Fill(2);

                for (int i=0; i<nElectron; i++)
                {
                        for (int j=i+1; j<nElectron; j++)
                        {
                                if(fabs(Electron_pt[i])<fabs(Electron_pt[j]))
                                {
                                        float a = Electron_pt[i];
                                        Electron_pt[i] = Electron_pt[j];
                                        Electron_pt[j] = a;

                                        a = Electron_eta[i];
                                        Electron_eta[i] = Electron_eta[j];
                                        Electron_eta[j] = a;

                                        a = Electron_phi[i];
                                        Electron_phi[i] = Electron_phi[j];
                                        Electron_phi[j] = a;

                                        a = Electron_e[i];
                                        Electron_e[i] = Electron_e[j];
                                        Electron_e[j] = a;
                                }
                        }
                }

                if(fabs(Electron_pt[0])<24) continue;
                if(fabs(Electron_pt[1])<13) continue;
                
		cut_flow_el->Fill(3);

		isEle = true; }


	  isMu = false;

          if (hlt_Mu17_Mu8) {
          cut_flow_mu->Fill(1);

          for(int i=0; i<nMuon; i++)
          {
		  if (!Muon_isPF[i]) continue;
                  if (Muon_pfiso[i]>=0.15) continue;
                  if(fabs(Muon_eta[i]) > 2.5) continue;                                                                  //muon channel

                  Muon_pt[nmu] = Muon_pt[i];
                  Muon_eta[nmu] = Muon_eta[i];
                  Muon_phi[nmu] = Muon_phi[i];
                  Muon_p[nmu] = Muon_p[i];
	//	  Muon_e[nmu] = Muon_e[i];

                  if (++nmu >= njetmx) break;
          }

          nMuon = nmu;
          if (nMuon < 2) continue;

//	  if (nMuon >= 2) {
	  cut_flow_mu->Fill(2);

	  for (int i=0; i<nMuon; i++)
                {
                        for (int j=i+1; j<nMuon; j++)
                        {
                                if(fabs(Muon_pt[i])<fabs(Muon_pt[j]))
                                {
                                        float a = Muon_pt[i];
                                        Muon_pt[i] = Muon_pt[j];
                                        Muon_pt[j] = a;

                                        a = Muon_eta[i];
                                        Muon_eta[i] = Muon_eta[j];
                                        Muon_eta[j] = a;

                                        a = Muon_phi[i];
                                        Muon_phi[i] = Muon_phi[j];
                                        Muon_phi[j] = a;

                                        a = Muon_p[i];
                                        Muon_p[i] = Muon_p[j];
                                        Muon_p[j] = a;

		/*			a = Muon_e[i];
                                        Muon_e[i] = Muon_e[j];
                                        Muon_e[j] = a;         */
                                }
                        }
                }

                if(fabs(Muon_pt[0])<18) continue;
                if(fabs(Muon_pt[1])<9) continue;

                cut_flow_mu->Fill(3);

                isMu = true; }


                if (isEle) {
                lepton1vec.SetPtEtaPhiE(fabs(Electron_pt[0]), Electron_eta[0], Electron_phi[0], Electron_e[0]);
                lepton2vec.SetPtEtaPhiE(fabs(Electron_pt[1]), Electron_eta[1], Electron_phi[1], Electron_e[1]);
		}

		if (isMu) {
                lepton1vec.SetPtEtaPhiE(fabs(Muon_pt[0]), Muon_eta[0], Muon_phi[0], fabs(Muon_p[0]));
                lepton2vec.SetPtEtaPhiE(fabs(Muon_pt[1]), Muon_eta[1], Muon_phi[1], fabs(Muon_p[1]));
                }    

		


	//-------------------------------------------------------------------------------AK4 jet info-----------------------------------------------------------------------------------//


		
	int indx = 0;

          for(int i=0; i<nPFJetAK4; i++)
          {
                  Jet_pt_nom[i] = PFJetAK4_pt[i];
                  Jet_mass_nom[i] = PFJetAK4_mass[i];

                  Jet_pt_jesup[i] = PFJetAK4_pt[i];
                  Jet_mass_jesup[i] = PFJetAK4_mass[i];

                  Jet_pt_jesdn[i] = PFJetAK4_pt[i];
                  Jet_mass_jesdn[i] = PFJetAK4_mass[i];

                  Jet_pt_resoup[i] = PFJetAK4_pt[i];
                  Jet_mass_resoup[i] = PFJetAK4_mass[i];

                  Jet_pt_resodn[i] = PFJetAK4_pt[i];
                  Jet_mass_resodn[i] = PFJetAK4_mass[i];


                  Jet_pt_nom[i] *= PFJetAK4_JEC[i];
                  Jet_mass_nom[i] *= PFJetAK4_JEC[i];
                  Jet_pt_nom[i] *= (1.+PFJetAK4_reso[i]);
                  Jet_mass_nom[i] *= (1.+PFJetAK4_reso[i]);

                  Jet_pt_jesup[i] *= PFJetAK4_jesup_Total[i];
                  Jet_mass_jesup[i] *= PFJetAK4_jesup_Total[i];
                  Jet_pt_jesup[i] *= (1.+PFJetAK4_reso[i]);
                  Jet_mass_jesup[i] *= (1.+PFJetAK4_reso[i]);

                  Jet_pt_jesdn[i] *= PFJetAK4_jesdn_Total[i];
                  Jet_mass_jesdn[i] *= PFJetAK4_jesdn_Total[i];
                  Jet_pt_jesdn[i] *= (1.+PFJetAK4_reso[i]);
                  Jet_mass_jesdn[i] *= (1.+PFJetAK4_reso[i]);

                  Jet_pt_resoup[i] *= PFJetAK4_JEC[i];
                  Jet_mass_resoup[i] *= PFJetAK4_JEC[i];
                  Jet_pt_resoup[i] *= (1.+PFJetAK4_resoup[i]);
                  Jet_mass_resoup[i] *= (1.+PFJetAK4_resoup[i]);

                  Jet_pt_resodn[i] *= PFJetAK4_JEC[i];
                  Jet_mass_resodn[i] *= PFJetAK4_JEC[i];
                  Jet_pt_resodn[i] *= (1.+PFJetAK4_resodn[i]);
                  Jet_mass_resodn[i] *= (1.+PFJetAK4_resodn[i]);


		  if(!PFJetAK4_jetID_tightlepveto[i]) continue;
                  if(fabs(PFJetAK4_y[i]) > 2.5) continue;
       //           if(Jet_pt_nom[i]<20) continue;
          //        if (PFJetAK4_btag_DeepFlav[i] < 0.2783) continue;

                  TLorentzVector b_vec;
                  b_vec.SetPtEtaPhiM(Jet_pt_nom[i], PFJetAK4_eta[i], PFJetAK4_phi[i], Jet_mass_nom[i]);

                  if ((lepton1vec.DeltaR(b_vec) < 0.4) || (lepton2vec.DeltaR(b_vec) < 0.4)) continue;

                  if(isMC){

                        BTagEntry::JetFlavor btv_flav;
                        if(abs(PFJetAK4_hadronflav[i])==5){ btv_flav = BTagEntry::FLAV_B; }
                        else if (abs(PFJetAK4_hadronflav[i])==4){ btv_flav = BTagEntry::FLAV_C; }
                        else { btv_flav = BTagEntry::FLAV_UDSG; }

                        PFJetAK4_btag_DeepFlav_SF[i] = reader_deepflav.eval_auto_bounds("central",btv_flav,fabs(PFJetAK4_y[i]),Jet_pt_nom[i]);
                        PFJetAK4_btag_DeepFlav_SF_up[i] = reader_deepflav.eval_auto_bounds("up",btv_flav,fabs(PFJetAK4_y[i]),Jet_pt_nom[i]);
                        PFJetAK4_btag_DeepFlav_SF_dn[i] = reader_deepflav.eval_auto_bounds("down",btv_flav,fabs(PFJetAK4_y[i]),Jet_pt_nom[i]);

                  }

                  else{

                        PFJetAK4_btag_DeepFlav_SF[i] = PFJetAK4_btag_DeepFlav_SF_up[i] = PFJetAK4_btag_DeepFlav_SF_dn[i] = 0;

                        }


		  Jet_pt_nom[indx] = Jet_pt_nom[i];
                  Jet_mass_nom[indx] = Jet_mass_nom[i];
                  Jet_pt_jesup[indx] = Jet_pt_jesup[i];
                  Jet_mass_jesup[indx] = Jet_mass_jesup[i];
                  Jet_pt_jesdn[indx] = Jet_pt_jesdn[i];
                  Jet_mass_jesdn[indx] = Jet_mass_jesdn[i];
                  Jet_pt_resoup[indx] = Jet_pt_resoup[i];
                  Jet_mass_resoup[indx] = Jet_mass_resoup[i];
                  Jet_pt_resodn[indx] = Jet_pt_resodn[i];
                  Jet_mass_resodn[indx] = Jet_mass_resodn[i];

                  PFJetAK4_y[indx] = PFJetAK4_y[i];
                  PFJetAK4_eta[indx] = PFJetAK4_eta[i];
                  PFJetAK4_phi[indx] = PFJetAK4_phi[i];
                  PFJetAK4_btag_DeepFlav[indx] = PFJetAK4_btag_DeepFlav[i];
                  PFJetAK4_btag_DeepFlav_SF[indx] = PFJetAK4_btag_DeepFlav_SF[i];
                  PFJetAK4_btag_DeepFlav_SF_up[indx] = PFJetAK4_btag_DeepFlav_SF_up[i];
                  PFJetAK4_btag_DeepFlav_SF_dn[indx] = PFJetAK4_btag_DeepFlav_SF_dn[i];
		  PFJetAK4_hadronflav[indx] = PFJetAK4_hadronflav[i];
                  PFJetAK4_partonflav[indx] = PFJetAK4_partonflav[i];
                  PFJetAK4_PUID[indx] = PFJetAK4_PUID[i];

                  if (++indx >= njetmx) break;
          }


          nPFJetAK4 = indx;

          if (nPFJetAK4 < 2) continue;

          if (isEle)
          {
                  cut_flow_el->Fill(4);
          }
          if (isMu)
          {
                  cut_flow_mu->Fill(4);
          }


	  for (int i=0; i<nPFJetAK4; i++)
          {
                  for (int j=i+1; j<nPFJetAK4; j++)
                  {
                           if(PFJetAK4_btag_DeepFlav[i]<PFJetAK4_btag_DeepFlav[j])
                           {
                                    float a = PFJetAK4_btag_DeepFlav[i];
                                    PFJetAK4_btag_DeepFlav[i] = PFJetAK4_btag_DeepFlav[j];
                                    PFJetAK4_btag_DeepFlav[j] = a;

                                    a = Jet_pt_nom[i];
                                    Jet_pt_nom[i] = Jet_pt_nom[j];                   //nominal
                                    Jet_pt_nom[j] = a;

                                    a = PFJetAK4_y[i];
                                    PFJetAK4_y[i] = PFJetAK4_y[j];
                                    PFJetAK4_y[j] = a;

                                    a = PFJetAK4_eta[i];
                                    PFJetAK4_eta[i] = PFJetAK4_eta[j];
                                    PFJetAK4_eta[j] = a;

                                    a = PFJetAK4_phi[i];
                                    PFJetAK4_phi[i] = PFJetAK4_phi[j];
                                    PFJetAK4_phi[j] = a;

                                    a = Jet_mass_nom[i];
                                    Jet_mass_nom[i] = Jet_mass_nom[j];               //nominal
                                    Jet_mass_nom[j] = a;

                                    a = PFJetAK4_btag_DeepFlav_SF[i];
                                    PFJetAK4_btag_DeepFlav_SF[i] = PFJetAK4_btag_DeepFlav_SF[j];
                                    PFJetAK4_btag_DeepFlav_SF[j] = a;

                                    a = PFJetAK4_btag_DeepFlav_SF_up[i];
                                    PFJetAK4_btag_DeepFlav_SF_up[i] = PFJetAK4_btag_DeepFlav_SF_up[j];
                                    PFJetAK4_btag_DeepFlav_SF_up[j] = a;

                                    a = PFJetAK4_btag_DeepFlav_SF_dn[i];
                                    PFJetAK4_btag_DeepFlav_SF_dn[i] = PFJetAK4_btag_DeepFlav_SF_dn[j];
                                    PFJetAK4_btag_DeepFlav_SF_dn[j] = a;

				    a = PFJetAK4_hadronflav[i];
                                    PFJetAK4_hadronflav[i] = PFJetAK4_hadronflav[j];
                                    PFJetAK4_hadronflav[j] = a;

                                    a = PFJetAK4_partonflav[i];
                                    PFJetAK4_partonflav[i] = PFJetAK4_partonflav[j];
                                    PFJetAK4_partonflav[j] = a;

                                    a = PFJetAK4_PUID[i];
                                    PFJetAK4_PUID[i] = PFJetAK4_PUID[j];
                                    PFJetAK4_PUID[j] = a;

				    a = Jet_pt_jesup[i];
                                    Jet_pt_jesup[i] = Jet_pt_jesup[j];                    //jesup
                                    Jet_pt_jesup[j] = a;

                                    a = Jet_mass_jesup[i];
                                    Jet_mass_jesup[i] = Jet_mass_jesup[j];                //jesup
                                    Jet_mass_jesup[j] = a;

                                    a = Jet_pt_jesdn[i];
                                    Jet_pt_jesdn[i] = Jet_pt_jesdn[j];                 //jesdn
                                    Jet_pt_jesdn[j] = a;

                                    a = Jet_mass_jesdn[i];
                                    Jet_mass_jesdn[i] = Jet_mass_jesdn[j];               //jesdn
                                    Jet_mass_jesdn[j] = a;

                                    a = Jet_pt_resoup[i];
                                    Jet_pt_resoup[i] = Jet_pt_resoup[j];                    //resoup
                                    Jet_pt_resoup[j] = a;

                                    a = Jet_mass_resoup[i];
                                    Jet_mass_resoup[i] = Jet_mass_resoup[j];                 //resoup
                                    Jet_mass_resoup[j] = a;

                                    a = Jet_pt_resodn[i];
                                    Jet_pt_resodn[i] = Jet_pt_resodn[j];                  //resodn
                                    Jet_pt_resodn[j] = a;

                                    a = Jet_mass_resodn[i];
                                    Jet_mass_resodn[i] = Jet_mass_resodn[j];                      //resodn
                                    Jet_mass_resodn[j] = a;
                           }
                   }
           }


	  TLorentzVector leadbvec_nom, subleadbvec_nom, leadbvec_jesup, subleadbvec_jesup, leadbvec_jesdn, subleadbvec_jesdn, leadbvec_resoup, subleadbvec_resoup, leadbvec_resodn, subleadbvec_resodn;


          bb = 0;
          ab = 0;
          aa = 0;

          if (Jet_pt_nom[0] >= 20 && Jet_pt_nom[1] >= 20 && (PFJetAK4_btag_DeepFlav[0] >= 0.2783 && PFJetAK4_btag_DeepFlav[1] >= 0.2783)) {bb = 1;}
          if (Jet_pt_nom[0] >= 20 && Jet_pt_nom[1] >= 20 && (PFJetAK4_btag_DeepFlav[0] >= 0.2783 && PFJetAK4_btag_DeepFlav[1] > 0.0490 && PFJetAK4_btag_DeepFlav[1] < 0.2783)) {ab = 1;}
          if (Jet_pt_nom[0] >= 20 && Jet_pt_nom[1] >= 20 && (PFJetAK4_btag_DeepFlav[0] < 0.2783 && PFJetAK4_btag_DeepFlav[1] < 0.2783)) {aa = 1;}


          if (isEle)
          {
                  if (bb == 1) {cut_flow_el->Fill(5);}
                  if (ab == 1) {cut_flow_el->Fill(6);}
          }

          if (isMu)
          {
                  if (bb == 1) {cut_flow_mu->Fill(5);}
                  if (ab == 1) {cut_flow_mu->Fill(6);}
          }


	  if (Jet_pt_nom[0] >= 20 && Jet_pt_nom[1] >= 20) {
          if (Jet_pt_nom[0] >= Jet_pt_nom[1])
          {
                   leadbvec_nom.SetPtEtaPhiM(Jet_pt_nom[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_nom[0]);
                   subleadbvec_nom.SetPtEtaPhiM(Jet_pt_nom[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_nom[1]);
          }                                                                                                                    //nominal

          else
          {
                   subleadbvec_nom.SetPtEtaPhiM(Jet_pt_nom[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_nom[0]);
                   leadbvec_nom.SetPtEtaPhiM(Jet_pt_nom[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_nom[1]);
          } }

          if (Jet_pt_jesup[0] >=20 && Jet_pt_jesup[1] >= 20) {
          if (Jet_pt_jesup[0] >= Jet_pt_jesup[1])
          {
                   leadbvec_jesup.SetPtEtaPhiM(Jet_pt_jesup[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_jesup[0]);
                   subleadbvec_jesup.SetPtEtaPhiM(Jet_pt_jesup[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_jesup[1]);
          }                                                                                                                           //jesup

          else
          {
                   subleadbvec_jesup.SetPtEtaPhiM(Jet_pt_jesup[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_jesup[0]);
                   leadbvec_jesup.SetPtEtaPhiM(Jet_pt_jesup[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_jesup[1]);
          } }

          if (Jet_pt_jesdn[0] >= 20 && Jet_pt_jesdn[1] >= 20) {
          if (Jet_pt_jesdn[0] >= Jet_pt_jesdn[1])
          {
                   leadbvec_jesdn.SetPtEtaPhiM(Jet_pt_jesdn[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_jesdn[0]);
                   subleadbvec_jesdn.SetPtEtaPhiM(Jet_pt_jesdn[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_jesdn[1]);
          }                                                                                                                           //jesdn

          else
          {
                   subleadbvec_jesdn.SetPtEtaPhiM(Jet_pt_jesdn[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_jesdn[0]);
                   leadbvec_jesdn.SetPtEtaPhiM(Jet_pt_jesdn[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_jesdn[1]);
          } }

	  if (Jet_pt_resoup[0] >= 20 && Jet_pt_resoup[1] >= 20) {
          if (Jet_pt_resoup[0] >= Jet_pt_resoup[1])
          {
                   leadbvec_resoup.SetPtEtaPhiM(Jet_pt_resoup[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_resoup[0]);
                   subleadbvec_resoup.SetPtEtaPhiM(Jet_pt_resoup[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_resoup[1]);
          }                                                                                                                               //resoup

          else
          {
                   subleadbvec_resoup.SetPtEtaPhiM(Jet_pt_resoup[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_resoup[0]);
                   leadbvec_resoup.SetPtEtaPhiM(Jet_pt_resoup[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_resoup[1]);
          } }

          if (Jet_pt_resodn[0] >= 20 && Jet_pt_resodn[1] >= 20) {
          if (Jet_pt_resodn[0] >= Jet_pt_resodn[1])
          {
                   leadbvec_resodn.SetPtEtaPhiM(Jet_pt_resodn[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_resodn[0]);
                   subleadbvec_resodn.SetPtEtaPhiM(Jet_pt_resodn[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_resodn[1]);
          }                                                                                                                               //resodn

          else
          {
                   subleadbvec_resodn.SetPtEtaPhiM(Jet_pt_resodn[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_resodn[0]);
                   leadbvec_resodn.SetPtEtaPhiM(Jet_pt_resodn[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_resodn[1]);
          } }



	//-------------------------------------------------------------------------------photon info-----------------------------------------------------------------------------//
	
				

	  indx = 0;

          for(int i=0; i<nPhoton; i++)
          {
        //        if (!Photon_mvaid_Fall17V2_WP90[i]) continue;
                  if (Photon_pt[i] < 20) continue;
                  if (fabs(Photon_eta[i]) > 2.5) continue;
                  if (fabs(Photon_eta[i]) > 1.44 && fabs(Photon_eta[i]) < 1.57) continue;
                  if (Photon_hadbyem[i] >= 0.08) continue;
                  if (Photon_e9by25[i] < 0.8) continue;

                  TLorentzVector phovec;
                  phovec.SetPtEtaPhiE(Photon_pt[i], Photon_eta[i], Photon_phi[i], Photon_e[i]);

                  if ((lepton1vec.DeltaR(phovec) < 0.4) || (lepton2vec.DeltaR(phovec) < 0.4)) continue;
                  if ((leadbvec_nom.DeltaR(phovec) < 0.4) || (subleadbvec_nom.DeltaR(phovec) < 0.4) || (leadbvec_jesup.DeltaR(phovec) < 0.4) || (subleadbvec_jesup.DeltaR(phovec) < 0.4) || (leadbvec_jesdn.DeltaR(phovec) < 0.4) || (subleadbvec_jesdn.DeltaR(phovec) < 0.4) || (leadbvec_resoup.DeltaR(phovec) < 0.4) || (subleadbvec_resoup.DeltaR(phovec) < 0.4) || (leadbvec_resodn.DeltaR(phovec) < 0.4) || (subleadbvec_resodn.DeltaR(phovec) < 0.4)) continue;

                  Photon_pt[indx] = Photon_pt[i];
                  Photon_e[indx] = Photon_e[i];
                  Photon_eta[indx] = Photon_eta[i];
                  Photon_phi[indx] = Photon_phi[i];

                  if (++indx >= njetmx) break;
          }

          nPhoton = indx;
          if (nPhoton < 2) continue;

          if (isEle)
          {
                  if (bb == 1) {cut_flow_el->Fill(7);}
                  if (ab == 1) {cut_flow_el->Fill(8);}
          }

          if (isMu)
          {
                  if (bb == 1) {cut_flow_mu->Fill(7);}
                  if (ab == 1) {cut_flow_mu->Fill(8);}
          }


	  for (int i=0; i<nPhoton; i++)
          {
                  for (int j=i+1; j<nPhoton; j++)
                  {
                          if(Photon_pt[i]<Photon_pt[j])
                          {
                                  float a = Photon_e[i];
                                  Photon_e[i] = Photon_e[j];
                                  Photon_e[j] = a;

                                  a = Photon_pt[i];
                                  Photon_pt[i] = Photon_pt[j];
                                  Photon_pt[j] = a;

                                  a = Photon_eta[i];
                                  Photon_eta[i] = Photon_eta[j];
                                  Photon_eta[j] = a;

                                  a = Photon_phi[i];
                                  Photon_phi[i] = Photon_phi[j];
                                  Photon_phi[j] = a;
                          }
                  }
          }

          TLorentzVector leadphovec, subleadphovec;

          leadphovec.SetPtEtaPhiE(Photon_pt[0], Photon_eta[0], Photon_phi[0], Photon_e[0]);
          subleadphovec.SetPtEtaPhiE(Photon_pt[1], Photon_eta[1], Photon_phi[1], Photon_e[1]);

		

	//----------------------------------------------------------------------- all SFs --------------------------------------------------------------------------//	



	  double leptonsf_weight = 1.0;
          double leptonsf_weight_up = 1.0;
          double leptonsf_weight_dn = 1.0;
          double leptonsf_weight_stat = 1.0;
          double leptonsf_weight_syst = 1.0;

          double puWeight, puWeightup, puWeightdown;
          puWeight = puWeightup = puWeightdown = 1.0;

          if(isMC){

                if(npu_vert_true>=0 && npu_vert_true<100){
                        float *puweights = Get_PU_Weights(file_pu_ratio, npu_vert_true);
                        puWeight = puweights[0];
                        puWeightup = puweights[1];
                        puWeightdown = puweights[2];
                }

                if (isEle){
                for(unsigned lep=0; lep<nele; lep++){
                                float *sfvalues = Electron_SF(file_el_sf, Electron_pt[lep], Electron_eta[lep]);
                                leptonsf_weight *= sfvalues[0];
                                leptonsf_weight_up *= sfvalues[1];
                                leptonsf_weight_dn *= sfvalues[2];
                                leptonsf_weight_stat *= (sfvalues[0] + sqrt(sfvalues[3]*sfvalues[3] + sfvalues[4]*sfvalues[4]));  // like this for time being
                                leptonsf_weight_syst *= (sfvalues[0] + sqrt(sfvalues[5]*sfvalues[5] + sfvalues[6]*sfvalues[6] + sfvalues[7]*sfvalues[7] + sfvalues[8]*sfvalues[8]));  // like this for time being
                        } }

                if (isMu){
                for(unsigned lep=0; lep<nmu; lep++){
                                float *sfvalues;
                                sfvalues = Muon_SF(file_mu_sf, muon_id_name, Muon_pt[lep], Muon_eta[lep]);
                                leptonsf_weight *= *(sfvalues+0);
                                leptonsf_weight_up *= *(sfvalues+1);
                                leptonsf_weight_dn *= *(sfvalues+2);
                                leptonsf_weight_stat *= *(sfvalues+3);
                                leptonsf_weight_syst *= *(sfvalues+4);
                        } }

                }


	  double SF_Trig, SF_Trig_stat, SF_Trig_syst;
        SF_Trig = 1; SF_Trig_stat = 1; SF_Trig_syst = 1;

        if(!isDATA) {
                      if(isMu) {
                      if(lepton1vec.Pt() > 17 && lepton2vec.Pt() > 8)
                      {
			      int pt1bin = std::max(1, std::min(h_trigmumu->GetNbinsX(), h_trigmumu->GetXaxis()->FindBin(lepton1vec.Pt())));
				int pt2bin = std::max(1, std::min(h_trigmumu->GetNbinsY(), h_trigmumu->GetYaxis()->FindBin(lepton2vec.Pt())));
				SF_Trig    = TMath::Max(float(1.e-3),float(h_trigmumu->GetBinContent(pt1bin, pt2bin))) ;
				SF_Trig_stat    = TMath::Max(float(1.e-3),float(h_trigmumu->GetBinContent(pt1bin, pt2bin))) ;
				SF_Trig_syst    = TMath::Max(float(1.e-3),float(h_trigmumu->GetBinContent(pt1bin, pt2bin))) ;
        }
        }

                          if (isEle) {
                          if(lepton1vec.Pt() > 23 && lepton2vec.Pt() > 12)
                          {
				  int pt1bin = std::max(1, std::min(h_trigee->GetNbinsX(), h_trigee->GetXaxis()->FindBin(lepton1vec.Pt())));
                int pt2bin = std::max(1, std::min(h_trigee->GetNbinsY(), h_trigee->GetYaxis()->FindBin(lepton2vec.Pt())));
                SF_Trig    = TMath::Max(float(1.e-3),float(h_trigee->GetBinContent(pt1bin, pt2bin))) ;
				SF_Trig_stat    = TMath::Max(float(1.e-3),float(h_trigee->GetBinContent(pt1bin, pt2bin))) ;
				SF_Trig_syst    = TMath::Max(float(1.e-3),float(h_trigee->GetBinContent(pt1bin, pt2bin))) ;
        }
        }
}

        double SF_Trig_1_up = SF_Trig + sqrt((SF_Trig_stat - SF_Trig)*(SF_Trig_stat - SF_Trig) +  (SF_Trig_syst - SF_Trig)*(SF_Trig_syst - SF_Trig));
        double SF_Trig_1_dn = SF_Trig - sqrt((SF_Trig_stat - SF_Trig)*(SF_Trig_stat - SF_Trig) +  (SF_Trig_syst - SF_Trig)*(SF_Trig_syst - SF_Trig));

        double SF_Trig_2_up = SF_Trig + abs(SF_Trig-1);
        double SF_Trig_2_dn = SF_Trig - abs(SF_Trig-1);


	//-----------------------------------------------------------------------tree fill----------------------------------------------------------------------//
        

	  if (isEle || isMu) {
          lep1pt = lepton1vec.Pt();
          lep1y = lepton1vec.Rapidity();
          lep1eta = lepton1vec.Eta();
          lep1phi = lepton1vec.Phi();
          lep1e = lepton1vec.Energy();
	  lep2pt = lepton2vec.Pt();
          lep2y = lepton2vec.Rapidity();
          lep2eta = lepton2vec.Eta();
          lep2phi = lepton2vec.Phi();
          lep2e = lepton2vec.Energy();
	  dilep_invmass = (lepton1vec+lepton2vec).M();

	  b1pt = leadbvec_nom.Pt();
          b1y = leadbvec_nom.Rapidity();
          b1eta = leadbvec_nom.Eta();
          b1phi = leadbvec_nom.Phi();
          b1e = leadbvec_nom.Energy();
          b2pt = subleadbvec_nom.Pt();
          b2y = subleadbvec_nom.Rapidity();
          b2eta = subleadbvec_nom.Eta();
          b2phi = subleadbvec_nom.Phi();
          b2e = subleadbvec_nom.Energy();
          bb_inv_mass = (leadbvec_nom+subleadbvec_nom).M();
          b1_DeepFlv = PFJetAK4_btag_DeepFlav[0];
          b2_DeepFlv = PFJetAK4_btag_DeepFlav[1];
          b1_hadFlv = PFJetAK4_hadronflav[0];
          b2_hadFlv = PFJetAK4_hadronflav[1];
          b1_partFlv = PFJetAK4_partonflav[0];
          b2_partFlv = PFJetAK4_partonflav[1];
          b1_PUid = PFJetAK4_PUID[0];
          b2_PUid = PFJetAK4_PUID[1];

          b1pt_jesup = leadbvec_jesup.Pt();
          b1y_jesup = leadbvec_jesup.Rapidity();
          b1eta_jesup = leadbvec_jesup.Eta();
          b1phi_jesup = leadbvec_jesup.Phi();
          b1e_jesup = leadbvec_jesup.Energy();
          b2pt_jesup = subleadbvec_jesup.Pt();
          b2y_jesup = subleadbvec_jesup.Rapidity();
          b2eta_jesup = subleadbvec_jesup.Eta();
          b2phi_jesup = subleadbvec_jesup.Phi();
          b2e_jesup = subleadbvec_jesup.Energy();
          bb_inv_mass_jesup = (leadbvec_jesup+subleadbvec_jesup).M();

          b1pt_jesdn = leadbvec_jesdn.Pt();
          b1y_jesdn = leadbvec_jesdn.Rapidity();
          b1eta_jesdn = leadbvec_jesdn.Eta();
          b1phi_jesdn = leadbvec_jesdn.Phi();
          b1e_jesdn = leadbvec_jesdn.Energy();
          b2pt_jesdn = subleadbvec_jesdn.Pt();
          b2y_jesdn = subleadbvec_jesdn.Rapidity();
          b2eta_jesdn = subleadbvec_jesdn.Eta();
          b2phi_jesdn = subleadbvec_jesdn.Phi();
          b2e_jesdn = subleadbvec_jesdn.Energy();
          bb_inv_mass_jesdn = (leadbvec_jesdn+subleadbvec_jesdn).M();

	  b1pt_resoup = leadbvec_resoup.Pt();
          b1y_resoup = leadbvec_resoup.Rapidity();
          b1eta_resoup = leadbvec_resoup.Eta();
          b1phi_resoup = leadbvec_resoup.Phi();
          b1e_resoup = leadbvec_resoup.Energy();
          b2pt_resoup = subleadbvec_resoup.Pt();
          b2y_resoup = subleadbvec_resoup.Rapidity();
          b2eta_resoup = subleadbvec_resoup.Eta();
          b2phi_resoup = subleadbvec_resoup.Phi();
          b2e_resoup = subleadbvec_resoup.Energy();
          bb_inv_mass_resoup = (leadbvec_resoup+subleadbvec_resoup).M();

          b1pt_resodn = leadbvec_resodn.Pt();
          b1y_resodn = leadbvec_resodn.Rapidity();
          b1eta_resodn = leadbvec_resodn.Eta();
          b1phi_resodn = leadbvec_resodn.Phi();
          b1e_resodn = leadbvec_resodn.Energy();
          b2pt_resodn = subleadbvec_resodn.Pt();
          b2y_resodn = subleadbvec_resodn.Rapidity();
          b2eta_resodn = subleadbvec_resodn.Eta();
          b2phi_resodn = subleadbvec_resodn.Phi();
          b2e_resodn = subleadbvec_resodn.Energy();
          bb_inv_mass_resodn = (leadbvec_resodn+subleadbvec_resodn).M();

          invmassbbgg = (leadbvec_nom+subleadbvec_nom+leadphovec+subleadphovec).M();
          invmassbbgg_jesup = (leadbvec_jesup+subleadbvec_jesup+leadphovec+subleadphovec).M();
          invmassbbgg_jesdn = (leadbvec_jesdn+subleadbvec_jesdn+leadphovec+subleadphovec).M();
          invmassbbgg_resoup = (leadbvec_resoup+subleadbvec_resoup+leadphovec+subleadphovec).M();
          invmassbbgg_resodn = (leadbvec_resodn+subleadbvec_resodn+leadphovec+subleadphovec).M();

          bb_gg_dphi = (leadbvec_nom+subleadbvec_nom).DeltaPhi(leadphovec+subleadphovec);
          bb_gg_dphi_jesup = (leadbvec_jesup+subleadbvec_jesup).DeltaPhi(leadphovec+subleadphovec);
          bb_gg_dphi_jesdn = (leadbvec_jesdn+subleadbvec_jesdn).DeltaPhi(leadphovec+subleadphovec);
          bb_gg_dphi_resoup = (leadbvec_resoup+subleadbvec_resoup).DeltaPhi(leadphovec+subleadphovec);
          bb_gg_dphi_resodn = (leadbvec_resodn+subleadbvec_resodn).DeltaPhi(leadphovec+subleadphovec);

	  if(isMC){

 //         event_weight = 1.0;

          weight_nom = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * SF_Trig;

          weight_PU_up = event_weight * puWeightup * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1];
          weight_leptonsf_up = event_weight * puWeight * leptonsf_weight_up * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1];
          weight_bSF_up = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF_up[0] * PFJetAK4_btag_DeepFlav_SF_up[1];
          weight_trig_up = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * SF_Trig_1_up * SF_Trig_2_up;

          weight_PU_dn = event_weight * puWeightdown * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1];
          weight_leptonsf_dn = event_weight * puWeight * leptonsf_weight_dn * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1];
          weight_bSF_dn = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF_dn[0] * PFJetAK4_btag_DeepFlav_SF_dn[1];
          weight_trig_dn = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * SF_Trig_1_dn * SF_Trig_2_dn;
          }

          pho_no = nPhoton;
          pho1pt = leadphovec.Pt();
          pho1y = leadphovec.Rapidity();
          pho1eta = leadphovec.Eta();
          pho1phi = leadphovec.Phi();
          pho1e = leadphovec.Energy();
          pho2pt = subleadphovec.Pt();
          pho2y = subleadphovec.Rapidity();
          pho2eta = subleadphovec.Eta();
          pho2phi = subleadphovec.Phi();
          pho2e = subleadphovec.Energy();
          dipho_invmass = (leadphovec+subleadphovec).M(); }


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
