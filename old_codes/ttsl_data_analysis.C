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

    static const int njetmx = 20;
	static const int njetmxAK8 =10;
	static const int npartmx = 50;
	
	
 int npfjetAK8;
 float pfjetAK8pt[njetmxAK8], pfjetAK8y[njetmxAK8], pfjetAK8eta[njetmxAK8], pfjetAK8phi[njetmxAK8], pfjetAK8mass[njetmxAK8], pfjetAK8JEC[njetmxAK8], pfjetAK8reso[njetmxAK8], qscale, pfjetAK8DeepTag_TvsQCD[njetmxAK8], pfjetAK8sdmass[njetmxAK8], pfjetAK8tau1[njetmxAK8], pfjetAK8tau2[njetmxAK8], pfjetAK8tau3[njetmxAK8], pfjetAK8sub1btag[njetmxAK8], pfjetAK8sub2btag[njetmxAK8];
 float pfjetAK8sub1pt[njetmxAK8], pfjetAK8sub1eta[njetmxAK8], pfjetAK8sub1phi[njetmxAK8], pfjetAK8sub1mass[njetmxAK8];
 float pfjetAK8sub2pt[njetmxAK8], pfjetAK8sub2eta[njetmxAK8], pfjetAK8sub2phi[njetmxAK8], pfjetAK8sub2mass[njetmxAK8];
 bool pfjetAK8tightID[njetmxAK8], pfjetAK8hashadtop[njetmxAK8], pfjetAK8hashadtop_alldecay[njetmxAK8], pfjetAK8hasW_alldecay[njetmxAK8], pfjetAK8has_NoOtherTopContInTop[njetmxAK8], pfjetAK8has_NoOtherTopContInW[njetmxAK8];
 float pfjetAK8subjet1pt[njetmxAK8], pfjetAK8subjet1eta[njetmxAK8], pfjetAK8subjet1phi[njetmxAK8], pfjetAK8subjet1mass[njetmxAK8], pfjetAK8subjet2pt[njetmxAK8], pfjetAK8subjet2eta[njetmxAK8], pfjetAK8subjet2phi[njetmxAK8], pfjetAK8subjet2mass[njetmxAK8], pfjetAK8subjet3pt[njetmxAK8], pfjetAK8subjet3eta[njetmxAK8], pfjetAK8subjet3phi[njetmxAK8], pfjetAK8subjet3mass[njetmxAK8];

 int npfjetAK4;
 float pfjetAK4pt[njetmx], pfjetAK4y[njetmx], pfjetAK4eta[njetmx], pfjetAK4phi[njetmx], pfjetAK4mass[njetmx], pfjetAK4JEC[njetmx], pfjetAK4reso[njetmx];
 bool pfjetAK4tightID[njetmx];
 float pfjetAK4btag_CMVA[njetmx], pfjetAK4btag_CSV[njetmx], pfjetAK4btag_DeepCSV[njetmx], pfjetAK4btag_DeepCSV2[njetmx], pfjetAK4btag_DeepFlav[njetmx], pfjetAK4btag_DeepQCD[njetmx];

 float miset , misphi , sumEt, misetsig, chmiset, chmisphi, chmisetsig;

 int nelecs;
  bool elmvaid[njetmx], elmvaid_noIso[njetmx];
  float elpt[njetmx], eleta[njetmx], elphi[njetmx], ele[njetmx], elp[njetmx], eldxy[njetmx],  eldxy_sv[njetmx], eldz[njetmx], elhovere[njetmx], elqovrper[njetmx], elchi[njetmx], elemiso03[njetmx], elhadiso03[njetmx], elemiso04[njetmx], elhadiso04[njetmx], elhadisodepth03[njetmx], eltkpt03[njetmx], eltkpt04[njetmx], eleoverp[njetmx], elietaieta[njetmx], eletain[njetmx], elphiin[njetmx], elfbrem[njetmx], elchhadiso03[njetmx], elchhadiso04[njetmx], elnohits[njetmx], elmisshits[njetmx] ;
  float elchhadiso[njetmx], elneuhadiso[njetmx], elphoiso[njetmx], elpuchhadiso[njetmx], elpfiso[njetmx], elconvdist[njetmx], elconvdoct[njetmx];
  int elndf[njetmx];

  int nmuons;
  float muonp[njetmx], muone[njetmx], muonpt[njetmx], muoneta[njetmx], muonphi[njetmx], muondrbm[njetmx], muondz[njetmx], muonpter[njetmx], muonchi[njetmx], muonecal[njetmx], muonhcal[njetmx], muonemiso[njetmx], muonhadiso[njetmx], muontkpt03[njetmx], muontkpt05[njetmx];
  float muonposmatch[njetmx], muontrkink[njetmx], muonsegcom[njetmx], muonpfiso[njetmx], muontrkvtx[njetmx], muonhit[njetmx], muonpixhit[njetmx], muonmst[njetmx], muontrklay[njetmx], muonvalfrac[njetmx];
  int muonndf[njetmx];
  bool muonisPF[njetmx], muonisGL[njetmx], muonisTRK[njetmx];
  bool muonisGoodGL[njetmx], muonisMed[njetmx], muonisLoose[njetmx];

  int ngenparticles;
  int genpartstatus[npartmx], genpartpdg[npartmx], genpartmompdg[npartmx], genpartmomid[npartmx], genpartdaugno[npartmx];
  float genpartpt[npartmx], genparteta[npartmx], genpartphi[npartmx], genpartm[npartmx], genpartq[npartmx];
  bool genpartfromhard[npartmx], genpartfromhardbFSR[npartmx], genpartisPromptFinalState[npartmx], genpartisLastCopyBeforeFSR[npartmx];

  double weightev;

  int ihlt03;

  int nele, ak4j, bj, ak8j;
  float nelept, neleeta, ak4jpt, ak4jeta, ak8jpt, ak8jeta;
  float bjpt, bjeta;
  float mEt, TrMass;
  float Delphi;



int main(int argc, char *argv[])
{
	cout<<"Program started"<<endl;

        char fOut[50];
        string inputFile=argv[3];
        string path="/home/abala/cms/CMSSW_10_5_0/src/condor_job/";


        if(inputFile=="EGamma_2018A_v2_0000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018A_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018A_v2_0001.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018A_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018A_v2_0002.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018A_0002_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018A_v2_0003.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018A_0003_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018A_v2_0004.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018A_0004_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018A_0000_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/egamma_2018A_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018A_0001_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/egamma_2018A_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="egamma_2018B_0000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018B_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018B_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/egamma_2018B_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="egamma_2018C_0000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018C_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="egamma_2018C_0001.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018C_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018C_0000_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/egamma_2018C_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018C_0001_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/egamma_2018C_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="egamma_2018D_0000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018D_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="egamma_2018D_0001.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018D_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="egamma_2018D_0002.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018D_0002_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="egamma_2018D_0003.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018D_0003_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="egamma_2018D_0004.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018D_0004_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="egamma_2018D_0005.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018D_0005_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="egamma_2018D_0006.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018D_0006_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="egamma_2018D_0007.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018D_0007_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="egamma_2018D_0008.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018D_0008_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="egamma_2018D_0009.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/15_05_21/egamma_2018D_0009_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018D_0000_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/egamma_2018D_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018D_0001_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/egamma_2018D_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018D_0002_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/egamma_2018D_0002_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018D_0003_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/egamma_2018D_0003_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018D_0004_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/egamma_2018D_0004_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018D_0005_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/egamma_2018D_0005_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018D_0006_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/egamma_2018D_0006_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018D_0007_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/egamma_2018D_0007_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018D_0008_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/egamma_2018D_0008_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="EGamma_2018D_0009_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/egamma_2018D_0009_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="TTSL_0000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ttsl/06_05_21/TTSL_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="TTSL_0001.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ttsl/06_05_21/TTSL_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="DYJetsToLL_0001.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/DYJetsToLL/04_05_21/DYJetsToLL_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="DYJetsToLL_0000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/DYJetsToLL/04_05_21/DYJetsToLL_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="DYJetsToLL_0000_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/DYJetsToLL/26_05_21/DYJetsToLL_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="DYJetsToLL_0001_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/DYJetsToLL/26_05_21/DYJetsToLL_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="ST_t_channel_antitop_0000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_t_channel_antitop_4f_InclusiveDecays/04_05_21/ST_t_channel_antitop_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="ST_t_channel_antitop_0001.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_t_channel_antitop_4f_InclusiveDecays/04_05_21/ST_t_channel_antitop_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="ST_t-channel_antitop_0000_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_t_channel_antitop_4f_InclusiveDecays/01_06_21/ST_t_channel_antitop_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="ST_t-channel_antitop_0001_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_t_channel_antitop_4f_InclusiveDecays/01_06_21/ST_t_channel_antitop_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="ST_t_channel_top_0000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_t_channel_top_4f_InclusiveDecays/04_05_21/ST_t_channel_top_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="ST_t_channel_top_0001.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_t_channel_top_4f_InclusiveDecays/04_05_21/ST_t_channel_top_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="ST_t_channel_top_0002.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_t_channel_top_4f_InclusiveDecays/04_05_21/ST_t_channel_top_0002_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="ST_t-channel_top_0000_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_t_channel_top_4f_InclusiveDecays/01_06_21/ST_t_channel_top_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="ST_t-channel_top_0001_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_t_channel_top_4f_InclusiveDecays/01_06_21/ST_t_channel_top_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="ST_t-channel_top_0002_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_t_channel_top_4f_InclusiveDecays/01_06_21/ST_t_channel_top_0002_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="WJetsToLNu.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/WJetsToLNu/04_05_21/WJetsToLNu_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="WJetsToLNu_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/WJetsToLNu/26_05_21/WJetsToLNu_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="ST_tW_antitop_5f_inclusiveDecays.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_tW_antitop_5f_inclusiveDecays/04_05_21/ST_tW_antitop_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="ST_tW_antitop_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_tW_antitop_5f_inclusiveDecays/01_06_21/ST_tW_antitop_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="ST_tW_top_5f_inclusiveDecays.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_tW_top_5f_inclusiveDecays/04_05_21/ST_tW_top_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="ST_tW_top_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_tW_top_5f_inclusiveDecays/01_06_21/ST_tW_top_0000_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="TTToHadronic_0000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/TTToHadronic/29_03_21/TTToHadronic_0000_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="TTToHadronic_0001.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/TTToHadronic/29_03_21/TTToHadronic_0001_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="TTToHadronic_0002.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/TTToHadronic/29_03_21/TTToHadronic_0002_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="TTToHadronic_0003.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/TTToHadronic/29_03_21/TTToHadronic_0003_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="TTTo2L2Nu.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/TTTo2L2Nu/10_04_21/TTTo2L2Nu_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="SinEle_2017B_280521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/SinEle_2017B_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="SinEle_2017C_0000_280521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/SinEle_2017C_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="SinEle_2017C_0001_280521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/SinEle_2017C_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="SinEle_2017D_280521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/SinEle_2017D_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="SinEle_2017E_0000_280521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/SinEle_2017E_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="SinEle_2017E_0001_280521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/SinEle_2017E_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="SinEle_2017F_0000_280521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/SinEle_2017F_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="SinEle_2017F_0001_280521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/egamma/01_06_21/SinEle_2017F_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="DY_2017_280521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/DYJetsToLL/28_05_21/DY_2017_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="t_antitop_2017_280521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_t_channel_antitop_4f_InclusiveDecays/01_06_21/t_antitop_2017_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="t_top_2017_280521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_t_channel_top_4f_InclusiveDecays/01_06_21/t_top_2017_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="tw_antitop_2017_280521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_tW_antitop_5f_inclusiveDecays/01_06_21/tw_antitop_2017_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="tw_top_2017_280521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ST_tW_top_5f_inclusiveDecays/01_06_21/tw_top_2017_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="TTSL_2017_280521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ttsl/01_06_21/TTSL_2017_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="WJet_2017_280521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/WJetsToLNu/28_05_21/WJet_2017_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="TTToSemiLeptonic_0000_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ttsl/01_06_21/TTbar_SemiLepton_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="TTToSemiLeptonic_0001_250521.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ttsl/01_06_21/TTbar_SemiLepton_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="TTbar_SemiLepton_z005beta005_0000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ttsl/z005beta005/TTbar_SemiLepton_z005beta005_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="TTbar_SemiLepton_z005beta005_0001.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ttsl/z005beta005/TTbar_SemiLepton_z005beta005_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="TTbar_SemiLepton_z01beta005_0000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ttsl/z01beta005/TTbar_SemiLepton_z01beta005_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="TTbar_SemiLepton_z01beta005_0001.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ttsl/z01beta005/TTbar_SemiLepton_z01beta005_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="TTbar_SemiLepton_z02beta000_0000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ttsl/z02beta000/TTbar_SemiLepton_z02beta000_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="TTbar_SemiLepton_z02beta000_0001.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ttsl/z02beta000/TTbar_SemiLepton_z02beta000_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="TTbar_SemiLepton_z01beta1_0000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ttsl/z01beta1/TTbar_SemiLepton_z01beta1_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="TTbar_SemiLepton_z01beta1_0001.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/ttsl/z01beta1/TTbar_SemiLepton_z01beta1_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="QCD_HT300to500.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/QCD_HT_300to500/15_06_21/QCD_HT_300to500_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="QCD_HT500to700_0000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/QCD_HT_500to700/15_06_21/QCD_HT_500to700_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="QCD_HT500to700_0001.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/QCD_HT_500to700/15_06_21/QCD_HT_500to700_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="QCD_HT700to1000_0000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/QCD_HT_700to1000/15_06_21/QCD_HT_700to1000_0000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="QCD_HT700to1000_0001.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/QCD_HT_700to1000/15_06_21/QCD_HT_700to1000_0001_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="QCD_HT1000to1500.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/QCD_HT_1000to1500/15_06_21/QCD_HT_1000to1500_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="QCD_HT1500to2000.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/QCD_HT_1500to2000/15_06_21/QCD_HT_1500to2000_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="QCD_HT2000toInf.log"){
                sprintf(fOut,"/home/abala/t3store3/top_tag/QCD_HT_2000toInf/15_06_21/QCD_HT_2000toInf_%s_%s.root",argv[1],argv[2]);
        }

        else{
                cout<<"Input file does not exist"<<endl;
                exit(0);
        }




        TFile *fout = new TFile(fOut,"RECREATE");
        TTree *Tout = new TTree("Tout", "Info");

        Tout->Branch("weightev", &weightev, "weightev/D");

	Tout->Branch("nele", &nele, "nele/I");
	Tout->Branch("nelept", &nelept, "nelept/F");
	Tout->Branch("neleeta", &neleeta, "neleeta/F");

	Tout->Branch("mEt", &mEt, "mEt/F");
	Tout->Branch("TrMass", &TrMass, "TrMass/F");

	Tout->Branch("ak4j", &ak4j, "ak4j/I");
	Tout->Branch("ak4jpt", &ak4jpt, "ak4jpt/F");
	Tout->Branch("ak4jeta", &ak4jeta, "ak4jeta/F");

	Tout->Branch("ak8j", &ak8j, "ak8j/I");
        Tout->Branch("ak8jpt", &ak8jpt, "ak8jpt/F");
        Tout->Branch("ak8jeta", &ak8jeta, "ak8jeta/F");

	Tout->Branch("bj", &bj, "bj/I");
        Tout->Branch("bjpt", &bjpt, "bjpt/F");
        Tout->Branch("bjeta", &bjeta, "bjeta/F");

        Tout->Branch("Delphi", &Delphi, "Delphi/F");




	TH1F *el_no = new TH1F("el_no","el_no",20,0,20);
        el_no->Sumw2();

        TH1F *el_pt = new TH1F("el_pt","el_pt",100,0,500);
        el_pt->Sumw2();

	TH1F *el_eta = new TH1F("el_eta","el_eta",100,-5,5);
        el_eta->Sumw2();

        TH1F *pfMET = new TH1F("pfMET","pfMET",100,0,500);
        pfMET->Sumw2();

        TH1F *trmass_rec = new TH1F("trmass_rec","trmass_rec",100,0,300);
        trmass_rec->Sumw2();

        TH1F *trmass_nocut = new TH1F("trmass_nocut","trmass_nocut",100,0,300);
        trmass_nocut->Sumw2();

	TH1F *pfjetAK4_pt = new TH1F("pfjetAK4_pt","pfjetAK4_pt",100,0,1500);
        pfjetAK4_pt->Sumw2();

        TH1F *pfjetAK4_mass = new TH1F("pfjetAK4_mass","pfjetAK4_mass",100,0,500);
        pfjetAK4_mass->Sumw2();

        TH1F *pfjetAK4_eta = new TH1F("pfjetAK4_eta","pfjetAK4_eta",100,-5,5);
        pfjetAK4_eta->Sumw2();

        TH1F *pfjetAK4_no = new TH1F("pfjetAK4_no","pfjetAK4_no",20,0,20);
        pfjetAK4_no->Sumw2();

	TH1F *pfjetAK8_pt = new TH1F("pfjetAK8_pt","pfjetAK8_pt",50,0,1500);
        pfjetAK8_pt->Sumw2();

        TH1F *pfjetAK8_mass = new TH1F("pfjetAK8_mass","pfjetAK8_mass",50,0,500);
        pfjetAK8_mass->Sumw2();

        TH1F *pfjetAK8_eta = new TH1F("pfjetAK8_eta","pfjetAK8_eta",50,-5,5);
        pfjetAK8_eta->Sumw2();

        TH1F *pfjetAK8_rap = new TH1F("pfjetAK8_rap","pfjetAK8_rap",50,-5,5);
        pfjetAK8_rap->Sumw2();

        TH1F *pfjetAK8_no = new TH1F("pfjetAK8_no","pfjetAK8_no",20,0,20);
        pfjetAK8_no->Sumw2();

	TH1F *ak8invmass = new TH1F("ak8invmass","ak8invmass",50,0,1500);
        ak8invmass->Sumw2();

        TH1F *topjet_pt = new TH1F("topjet_pt","topjet_pt",50,0,1500);
        topjet_pt->Sumw2();

        TH1F *subwpt = new TH1F("subwpt","subwpt",50,0,1500);
        subwpt->Sumw2();

        TH1F *subbpt = new TH1F("subbpt","subbpt",50,0,1500);
        subbpt->Sumw2();

        TH1F *subq1pt = new TH1F("subq1pt","subq1pt",50,0,1500);
        subq1pt->Sumw2();

        TH1F *subq2pt = new TH1F("subq2pt","subq2pt",50,0,1500);
        subq2pt->Sumw2();

        TH1F *ptfrac1 = new TH1F("ptfrac1","ptfrac1",20,0,1);
        ptfrac1->Sumw2();

        TH1F *ptfrac2 = new TH1F("ptfrac2","ptfrac2",20,0,1);
        ptfrac2->Sumw2();

        TH1F *sub1ptfrac = new TH1F("sub1ptfrac","sub1ptfrac",20,0,1);
        sub1ptfrac->Sumw2();

        TH1F *sub1mass = new TH1F("sub1mass","sub1mass",30,0,300);
        sub1mass->Sumw2();

        TH1F *sub2ptfrac = new TH1F("sub2ptfrac","sub2ptfrac",20,0,1);
        sub2ptfrac->Sumw2();

        TH1F *sub2mass = new TH1F("sub2mass","sub2mass",30,0,300);
        sub2mass->Sumw2();

        TH1F *sub3ptfrac = new TH1F("sub3ptfrac","sub3ptfrac",20,0,1);
        sub3ptfrac->Sumw2();

        TH1F *sub3mass = new TH1F("sub3mass","sub3mass",30,0,300);
        sub3mass->Sumw2();

        TH1F *sub1efrac = new TH1F("sub1efrac","sub1efrac",20,0,1);
        sub1efrac->Sumw2();

        TH1F *sub2efrac = new TH1F("sub2efrac","sub2efrac",20,0,1);
        sub2efrac->Sumw2();

        TH1F *sub3efrac = new TH1F("sub3efrac","sub3efrac",20,0,1);
        sub3efrac->Sumw2();

        TH1F *sub1efrac_topmasscut = new TH1F("sub1efrac_topmasscut","sub1efrac_topmasscut",20,0,1);
        sub1efrac_topmasscut->Sumw2();

        TH1F *sub2efrac_topmasscut = new TH1F("sub2efrac_topmasscut","sub2efrac_topmasscut",20,0,1);
        sub2efrac_topmasscut->Sumw2();

        TH1F *sub3efrac_topmasscut = new TH1F("sub3efrac_topmasscut","sub3efrac_topmasscut",20,0,1);
        sub3efrac_topmasscut->Sumw2();

	TH1F *efrac_j1_j2 = new TH1F("efrac_j1_j2","efrac_j1_j2",20,0,1);
        efrac_j1_j2->Sumw2();

        TH1F *efrac_j1j2 = new TH1F("efrac_j1j2","efrac_j1j2",20,0,1);
        efrac_j1j2->Sumw2();

        TH1F *zkdist = new TH1F("zkdist","zkdist",20,0,1);
        zkdist->Sumw2();

        TH1F *topjet_eta = new TH1F("topjet_eta","topjet_eta",50,-5,5);
        topjet_eta->Sumw2();

        TH1F *topjet_rap = new TH1F("topjet_rap","topjet_rap",50,-5,5);
        topjet_rap->Sumw2();

        TH1F *topjet_no = new TH1F("topjet_no","topjet_no",20,0,20);
        topjet_no->Sumw2();

        TH1F *topjet_mass = new TH1F("topjet_mass","topjet_mass",30,0,300);
        topjet_mass->Sumw2();

        TH1F *ak8mass_mat_top = new TH1F("ak8mass_mat_top","ak8mass_mat_top",50,0,500);
        ak8mass_mat_top->Sumw2();

	TH1F *recon_top_mass = new TH1F("recon_top_mass","recon_top_mass",30,0,300);
        recon_top_mass->Sumw2();
	 
	TH1F *recon_w_mass = new TH1F("recon_w_mass","recon_w_mass",30,0,300);
        recon_w_mass->Sumw2();

	TH1F *bj1_invmass = new TH1F("bj1_invmass","bj1_invmass",30,0,300);
        bj1_invmass->Sumw2();

	TH1F *bj2_invmass = new TH1F("bj2_invmass","bj2_invmass",30,0,300);
        bj2_invmass->Sumw2();

        TH1F *w_mass = new TH1F("w_mass","w_mass",30,0,300);
        w_mass->Sumw2();

        TH1F *submass1 = new TH1F("submass1","submass1",30,0,300);
        submass1->Sumw2();

        TH1F *submass2 = new TH1F("submass2","submass2",30,0,300);
        submass2->Sumw2();

        TH1F *topjet_sdmass = new TH1F("topjet_sdmass","topjet_sdmass",100,0,300);
        topjet_sdmass->Sumw2();

        TH1F *topjet_sdtaumass = new TH1F("topjet_sdtaumass","topjet_sdtaumass",100,0,300);
        topjet_sdtaumass->Sumw2();

	TH1F *bjet_pt = new TH1F("bjet_pt","bjet_pt",100,0,1500);
        bjet_pt->Sumw2();

        TH1F *bjet_eta = new TH1F("bjet_eta","bjet_eta",100,-5,5);
        bjet_eta->Sumw2();

        TH1F *nbjet = new TH1F("nbjet","nbjet",20,0,20);
        nbjet->Sumw2();

	TH1F *cut_flow = new TH1F("cut_flow","cut_flow",20,0,20);
        cut_flow->Sumw2();

        TH1F *el_dz = new TH1F("el_dz","el_dz",50,-2,2);
        el_dz->Sumw2();

        TH1F *cut_flow_1 = new TH1F("cut_flow_1","cut_flow_1",10,0,10);
        cut_flow_1->Sumw2();

        TH1F *cut_flow_top = new TH1F("cut_flow_top","cut_flow_top",10,0,10);
        cut_flow_top->Sumw2();

        TH1F *cut_flow_leading_w_or_b = new TH1F("cut_flow_leading_w_or_b","cut_flow_leading_w_or_b",10,0,10);
        cut_flow_leading_w_or_b->Sumw2();

	TH1F *cutflow_j1b_j2b = new TH1F("cutflow_j1b_j2b","cutflow_j1b_j2b",10,0,10);
        cutflow_j1b_j2b->Sumw2();

        TH1F *top_pt_check = new TH1F("top_pt_check","top_pt_check",10,0,10);
        top_pt_check->Sumw2();

        TH1F *q1q2_dphi = new TH1F("q1q2_dphi","q1q2_dphi",50,0,2);
        q1q2_dphi->Sumw2();

        TH1F *q1q2_dphi_topmasscut = new TH1F("q1q2_dphi_topmasscut","q1q2_dphi_topmasscut",50,0,2);
        q1q2_dphi_topmasscut->Sumw2();

        TH1F *reco_top_mass_ratio = new TH1F("reco_top_mass_ratio","reco_top_mass_ratio",50,-2,2);
        reco_top_mass_ratio->Sumw2();

        TH1F *z12dist = new TH1F("z12dist","z12dist",50,0,1);
        z12dist->Sumw2();

        TH1F *delR12dist = new TH1F("delR12dist","delR12dist",50,0,1);
        delR12dist->Sumw2();

        TH1F *ak8tau32 = new TH1F("ak8tau32","ak8tau32",50,0,1);
        ak8tau32->Sumw2();

        TH1F *ak8tau32withcut = new TH1F("ak8tau32withcut","ak8tau32withcut",20,0,1);
        ak8tau32withcut->Sumw2();

 //       TH2F* z12vsdelR12 = new TH2F("z12vsdelR12", "z12vsdelR12", 50, 0, 1, 50, 0, 1);
 //       z12vsdelR12->Sumw2();



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

   cout<<fileName<<endl;

//   cout<<"processedEvt="<<processedEvt<<endl;



        TTree *T1 = (TTree*)f->Get("T1");



	T1->SetBranchAddress("event_weight",&event_weight);
        T1->SetBranchAddress("ihlt03",&ihlt03);

	T1->SetBranchAddress("npfjetAK8",&npfjetAK8);
    T1->SetBranchAddress("pfjetAK8tightID",pfjetAK8tightID);
    T1->SetBranchAddress("pfjetAK8pt",pfjetAK8pt);
    T1->SetBranchAddress("pfjetAK8y",pfjetAK8y);
    T1->SetBranchAddress("pfjetAK8eta", pfjetAK8eta);
    T1->SetBranchAddress("pfjetAK8phi",pfjetAK8phi);
    T1->SetBranchAddress("pfjetAK8mass",pfjetAK8mass);
    T1->SetBranchAddress("pfjetAK8sdmass",pfjetAK8sdmass);
    T1->SetBranchAddress("pfjetAK8JEC",pfjetAK8JEC);
    T1->SetBranchAddress("pfjetAK8reso",pfjetAK8reso);
    T1->SetBranchAddress("pfjetAK8DeepTag_TvsQCD",pfjetAK8DeepTag_TvsQCD);
    T1->SetBranchAddress("qscale",&qscale);
    T1->SetBranchAddress("pfjetAK8tau1",pfjetAK8tau1);
    T1->SetBranchAddress("pfjetAK8tau2",pfjetAK8tau2);
    T1->SetBranchAddress("pfjetAK8tau3",pfjetAK8tau3);
    T1->SetBranchAddress("pfjetAK8sub1btag",pfjetAK8sub1btag);
    T1->SetBranchAddress("pfjetAK8sub2btag",pfjetAK8sub2btag);
    T1->SetBranchAddress("pfjetAK8sub1pt",pfjetAK8sub1pt);
    T1->SetBranchAddress("pfjetAK8sub1eta",pfjetAK8sub1eta);
    T1->SetBranchAddress("pfjetAK8sub1phi",pfjetAK8sub1phi);
    T1->SetBranchAddress("pfjetAK8sub1mass",pfjetAK8sub1mass);
    T1->SetBranchAddress("pfjetAK8sub2pt",pfjetAK8sub2pt);
    T1->SetBranchAddress("pfjetAK8sub2eta",pfjetAK8sub2eta);
    T1->SetBranchAddress("pfjetAK8sub2phi",pfjetAK8sub2phi);
    T1->SetBranchAddress("pfjetAK8sub2mass",pfjetAK8sub2mass);
    T1->SetBranchAddress("pfjetAK8subjet1pt",pfjetAK8subjet1pt);
    T1->SetBranchAddress("pfjetAK8subjet1eta",pfjetAK8subjet1eta);
    T1->SetBranchAddress("pfjetAK8subjet1phi",pfjetAK8subjet1phi);
    T1->SetBranchAddress("pfjetAK8subjet1mass",pfjetAK8subjet1mass);
    T1->SetBranchAddress("pfjetAK8subjet2pt",pfjetAK8subjet2pt);
    T1->SetBranchAddress("pfjetAK8subjet2eta",pfjetAK8subjet2eta);
    T1->SetBranchAddress("pfjetAK8subjet2phi",pfjetAK8subjet2phi);
    T1->SetBranchAddress("pfjetAK8subjet2mass",pfjetAK8subjet2mass);
    T1->SetBranchAddress("pfjetAK8subjet3pt",pfjetAK8subjet3pt);
    T1->SetBranchAddress("pfjetAK8subjet3mass",pfjetAK8subjet3mass);

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
  T1->SetBranchAddress("muoneta",muoneta);
  T1->SetBranchAddress("muonphi",muonphi);
  T1->SetBranchAddress("muone",muone);
  T1->SetBranchAddress("muonpfiso",muonpfiso);
  T1->SetBranchAddress("muonisPF",muonisPF);
  T1->SetBranchAddress("muonisMed",muonisMed);

    T1->SetBranchAddress("PFMET",&miset) ;
    T1->SetBranchAddress("PFMETPhi",&misphi) ;
    T1->SetBranchAddress("MisEtSig",&misetsig);/*
    T1->SetBranchAddress("PFCHMET",&chmiset) ;
    T1->SetBranchAddress("PFCHMETPhi",&chmisphi);
    T1->SetBranchAddress("CHMisEtSig",&chmisetsig);*/
    T1->SetBranchAddress("sumEt",&sumEt);

    T1->SetBranchAddress("ngenparticles",&ngenparticles);
  T1->SetBranchAddress("genpartstatus",genpartstatus);
  T1->SetBranchAddress("genpartpdg",genpartpdg);
  T1->SetBranchAddress("genpartmompdg",genpartmompdg);
//  T1->Branch("genpartmomid",genpartmomid);
  T1->SetBranchAddress("genpartdaugno",genpartdaugno);
  T1->SetBranchAddress("genpartfromhard",genpartfromhard);
  T1->SetBranchAddress("genpartfromhardbFSR",genpartfromhardbFSR);
  T1->SetBranchAddress("genpartisPromptFinalState",genpartisPromptFinalState);
  T1->SetBranchAddress("genpartisLastCopyBeforeFSR",genpartisLastCopyBeforeFSR);
  T1->SetBranchAddress("genpartpt",genpartpt);
  T1->SetBranchAddress("genparteta",genparteta);
  T1->SetBranchAddress("genpartphi",genpartphi);
  T1->SetBranchAddress("genpartm",genpartm);
  T1->SetBranchAddress("genpartq",genpartq);




  int nevents;

    double weight;
    double sumw = 0.0;
    double wt;
    int totevt = 0;

    nevents=T1->GetEntries();


//    cout << "tot no. of event = "<< totevt << endl;


    for(int iev=0; iev<nevents; iev++)
    {
	    T1->GetEntry(iev);
//	    double fwt = weight*fact;
	    if((TotEvt+iev)%10000==0)
	    {
                   cout<<iev<<" Events Analysed out of "<<nevents<<endl;
            }

            if (event_weight < 0) continue;
            weightev = event_weight;
            weight = event_weight;

            cut_flow->Fill(0);

            if (ihlt03==0) continue;
            cut_flow->Fill(1);



            vector<TLorentzVector> elecCol;
	    int EleIndex=-10;
            for (int iel=0; iel<nelecs; iel++)
            {
                    TLorentzVector v1;
		    v1.SetPtEtaPhiE(fabs(elpt[iel]), eleta[iel], elphi[iel], ele[iel]);
                    cut_flow_1->Fill(0);
                    el_dz->Fill(eldz[iel]);
                    if (!elmvaid[iel]) continue;
                    cut_flow_1->Fill(1);
                    if (elmisshits[iel]>1) continue;
                    cut_flow_1->Fill(2);
                    if (fabs(eldxy[iel])>0.045) continue;
                    cut_flow_1->Fill(3);
                    if (elpfiso[iel]>=0.15) continue;
                    cut_flow_1->Fill(4);
                    if(fabs(v1.Eta()) > 2.5) continue;
                    cut_flow_1->Fill(5);
                    if(v1.Pt()<33) continue;
                    cut_flow_1->Fill(6);
                    elecCol.push_back(v1);
		    EleIndex=iel;
            }


	    if(elecCol.size()>1 || elecCol.size()<1) continue;
	    cut_flow->Fill(2);

            double leadElePt=0, leadEleEta=-100, leadElePhi=-100;

            for(int iel=0; iel<elecCol.size(); iel++)
            {
                        if(leadElePt<elecCol[iel].Pt())
                        {
                                leadElePt=elecCol[iel].Pt();
                                leadEleEta=elecCol[iel].Eta();
                                leadElePhi=elecCol[iel].Phi();
                        }
            }

	    double nLooseEle=0;
	    for (int iel=0; iel<nelecs; iel++)
            {
		    if(iel==EleIndex)continue;
                    TLorentzVector v1;
                    v1.SetPtEtaPhiE(fabs(elpt[iel]), eleta[iel], elphi[iel], ele[iel]);
                    if (!elmvaid[iel]) continue;
                    if (elmisshits[iel]>1) continue;
                    if (fabs(eldxy[iel])>0.045) continue;
                    if (elpfiso[iel]>=0.2) continue;
                    if(fabs(v1.Eta()) > 2.5) continue;
                    if(v1.Pt()<15) continue;
		    nLooseEle++;
            }
	
	    if(nLooseEle>0)continue;
	
	    cut_flow->Fill(3);

            double nLooseMu=0;
            for (int imu=0; imu<nmuons; imu++)
            {
                    TLorentzVector v1;
                    v1.SetPtEtaPhiE(muonpt[imu], muoneta[imu], muonphi[imu], muone[imu]);
                    if (!muonisPF[imu]) continue;
                    if (muonpfiso[imu]>=0.2) continue;
                    if(fabs(v1.Eta()) > 2.5) continue;
                    if(v1.Pt()<15) continue;
                    nLooseMu++;
            }
 
            if(nLooseMu>0)continue;

            cut_flow->Fill(4);

	    if (miset < 50.0) continue;
	    cut_flow->Fill(5);

	    float dphi = PhiInRange(leadElePhi-misphi);
            float mt_rec = sqrt(2*leadElePt*miset*(1-cos(dphi)));
            trmass_nocut->Fill(mt_rec);
	    if (mt_rec < 30.0 || mt_rec > 130.0) continue;
	    cut_flow->Fill(6);


	    vector<TLorentzVector> jetCol;
            vector<TLorentzVector> bjetCol;
            vector<TLorentzVector> alljetCol;


	    for(int ijet=0; ijet<npfjetAK4; ijet++)
	    {
                    if(!pfjetAK4tightID[ijet]) continue;
                    pfjetAK4pt[ijet] *= pfjetAK4JEC[ijet];
                    pfjetAK4mass[ijet] *= pfjetAK4JEC[ijet];
                    pfjetAK4pt[ijet] *= (1.+pfjetAK4reso[ijet]);

		    TLorentzVector jet;
		    jet.SetPtEtaPhiM(pfjetAK4pt[ijet],pfjetAK4eta[ijet],pfjetAK4phi[ijet],pfjetAK4mass[ijet]);
                    if(jet.Pt() < 20 || fabs(jet.Eta()) > 2.5)continue;
                    alljetCol.push_back(jet);

		    if (pfjetAK4btag_DeepCSV[ijet] > 0.4184)
		    {
			    bjetCol.push_back(jet);
		    }
		    else
		    {
			    jetCol.push_back(jet);
		    }

	    }


	    if(bjetCol.size()<2)continue;
	    cut_flow->Fill(7);

            if(alljetCol.size()<4)continue;
	    cut_flow->Fill(8);


	    double leadJetPt=0, leadJetEta=-100, leadJetPhi=-100, leadJetMass=-100;

            for(int ijet=0; ijet<jetCol.size(); ijet++)
	    {
                        if(leadJetPt<jetCol[ijet].Pt())
			{
                                leadJetPt=jetCol[ijet].Pt();
                                leadJetEta=jetCol[ijet].Eta();
                                leadJetPhi=jetCol[ijet].Phi();
                                leadJetMass=jetCol[ijet].M();
                        }
            }


	    double leadbJetPt=0, leadbJetEta=-100, leadbJetPhi=-100;

	    for(int ijet=0; ijet<bjetCol.size(); ijet++)
            {
                        if(leadbJetPt<bjetCol[ijet].Pt())
                        {
                                leadbJetPt=bjetCol[ijet].Pt();
                                leadbJetEta=bjetCol[ijet].Eta();
                                leadbJetPhi=bjetCol[ijet].Phi();
                        }
            }

	    double leadallJetPt=0, leadallJetEta=-100, leadallJetPhi=-100, leadallJetMass=-100;

	    for(int ijet=0; ijet<alljetCol.size(); ijet++)
            {
                        if(leadallJetPt<alljetCol[ijet].Pt())
                        {
                                leadallJetPt=alljetCol[ijet].Pt();
                                leadallJetEta=alljetCol[ijet].Eta();
                                leadallJetPhi=alljetCol[ijet].Phi();
                                leadallJetMass=alljetCol[ijet].M();
                        }
            }


	    float delphi = PhiInRange(leadJetPhi-leadbJetPhi);


            vector<TLorentzVector> GenTopCol;
            vector<TLorentzVector> GenAntitopCol;

            for (int ij=0; ij<ngenparticles; ij++)
            {
                    TLorentzVector genv;
                    genv.SetPtEtaPhiM(genpartpt[ij], genparteta[ij], genpartphi[ij], genpartm[ij]);
                    if(genv.Pt() < 300 || fabs(genv.Eta()) > 2.5)continue;
                    if (genpartpdg[ij] == 6)
                    {
                            if(abs(genpartmompdg[ij])!=6)
                            {
                                    GenTopCol.push_back(genv);
                            }
                    }

                    if (genpartpdg[ij] == -6)
                    {
                            if(abs(genpartmompdg[ij])!=6)
                            {
                                    GenAntitopCol.push_back(genv);
                            }
                    }
            }

            double leadGenTopPt=0, leadGenTopEta=-100, leadGenTopPhi=-100, leadGenTopMass=-100;


            for(int ij=0; ij<GenTopCol.size(); ij++)
            {
                        if(leadGenTopPt<GenTopCol[ij].Pt())
                        {
                                leadGenTopPt=GenTopCol[ij].Pt();
                                leadGenTopEta=GenTopCol[ij].Eta();
                                leadGenTopPhi=GenTopCol[ij].Phi();
                                leadGenTopMass=GenTopCol[ij].M();
                        }
            }

            double leadGenAntiTopPt=0, leadGenAntiTopEta=-100, leadGenAntiTopPhi=-100, leadGenAntiTopMass=-100;


            for(int ij=0; ij<GenAntitopCol.size(); ij++)
            {
                        if(leadGenAntiTopPt<GenAntitopCol[ij].Pt())
                        {
                                leadGenAntiTopPt=GenAntitopCol[ij].Pt();
                                leadGenAntiTopEta=GenAntitopCol[ij].Eta();
                                leadGenAntiTopPhi=GenAntitopCol[ij].Phi();
                                leadGenAntiTopMass=GenAntitopCol[ij].M();
                        }
            }


            if (GenTopCol.size()==1)
            {
                       cut_flow_top->Fill(3);
            }


	    vector<TLorentzVector> AK8jetCol;
            vector<TLorentzVector> topjetCol;
            vector<TLorentzVector> sdtopjetCol;
            vector<TLorentzVector> sdtautopjetCol;



	    for(int ijet=0; ijet<npfjetAK8; ijet++)
            {
                    if(!pfjetAK8tightID[ijet]) continue;
                    pfjetAK8pt[ijet] *= pfjetAK8JEC[ijet];
                    pfjetAK8mass[ijet] *= pfjetAK8JEC[ijet];
                    pfjetAK8pt[ijet] *= (1.+pfjetAK8reso[ijet]);

                    TLorentzVector AK8jet;
                    AK8jet.SetPtEtaPhiM(pfjetAK8pt[ijet],pfjetAK8eta[ijet],pfjetAK8phi[ijet],pfjetAK8mass[ijet]);

                    if(AK8jet.Pt() < 20 || fabs(AK8jet.Eta()) > 2.5)continue;

                    AK8jetCol.push_back(AK8jet);

                    if (((pfjetAK8tau3[ijet]/pfjetAK8tau2[ijet]) < 0.54) && (pfjetAK8sdmass[ijet] > 105 && pfjetAK8sdmass[ijet] < 220) && (pfjetAK8sub1btag[ijet] > 0.4184 || pfjetAK8sub2btag[ijet] > 0.4184))
                    {
                            topjetCol.push_back(AK8jet);
                    }


                    if (pfjetAK8sdmass[ijet] > 105 && pfjetAK8sdmass[ijet] < 220)
                    {
                            sdtopjetCol.push_back(AK8jet);
                    }
        
                    if (((pfjetAK8tau3[ijet]/pfjetAK8tau2[ijet]) < 0.54) && (pfjetAK8sdmass[ijet] > 105 && pfjetAK8sdmass[ijet] < 220))
                    {
                            sdtautopjetCol.push_back(AK8jet);
                    }

              }

	    double leadAK8JetPt=0, leadAK8JetEta=-100, leadAK8JetPhi=-100, leadAK8JetMass=-100, leadAK8JetRap=-100;

            for(int ijet=0; ijet<AK8jetCol.size(); ijet++)
            {
                        if(leadAK8JetPt<AK8jetCol[ijet].Pt())
                        {
                                leadAK8JetPt=AK8jetCol[ijet].Pt();
                                leadAK8JetEta=AK8jetCol[ijet].Eta();
                                leadAK8JetRap=AK8jetCol[ijet].Rapidity();
                                leadAK8JetPhi=AK8jetCol[ijet].Phi();
                                leadAK8JetMass=AK8jetCol[ijet].M();
                        }
            }


            int njet = 0;

            for(int ijet=0; ijet<npfjetAK8; ijet++)
            {

                        if(!pfjetAK8tightID[ijet]) continue;
                        pfjetAK8pt[ijet] *= pfjetAK8JEC[ijet];
                        pfjetAK8mass[ijet] *= pfjetAK8JEC[ijet];
                        pfjetAK8pt[ijet] *= (1.+pfjetAK8reso[ijet]);
                        if(fabs(pfjetAK8eta[ijet]) > 2.5) continue;
                        if(pfjetAK8pt[ijet]<20) continue;

                        pfjetAK8pt[njet]  = pfjetAK8pt[ijet];
                        pfjetAK8eta[njet] = pfjetAK8eta[ijet];
                        pfjetAK8phi[njet]  = pfjetAK8phi[ijet];
                        pfjetAK8mass[njet]  = pfjetAK8mass[ijet];
                        pfjetAK8tau3[njet] = pfjetAK8tau3[ijet];
                        pfjetAK8tau2[njet] = pfjetAK8tau2[ijet];
                        pfjetAK8sdmass[njet] = pfjetAK8sdmass[ijet];
                        pfjetAK8sub1btag[njet] = pfjetAK8sub1btag[ijet];
                        pfjetAK8sub2btag[njet] = pfjetAK8sub2btag[ijet];
                        pfjetAK8sub1pt[njet] = pfjetAK8sub1pt[ijet];
                        pfjetAK8sub2pt[njet] = pfjetAK8sub2pt[ijet];
                        pfjetAK8sub2mass[njet] = pfjetAK8sub2mass[ijet];
                        pfjetAK8sub1mass[njet] = pfjetAK8sub1mass[ijet];
                        pfjetAK8subjet1pt[njet] = pfjetAK8subjet1pt[ijet];
                        pfjetAK8subjet1mass[njet] = pfjetAK8subjet1mass[ijet];
                        pfjetAK8subjet2pt[njet] = pfjetAK8subjet2pt[ijet];
                        pfjetAK8subjet2mass[njet] = pfjetAK8subjet2mass[ijet];
                        pfjetAK8subjet3pt[njet] = pfjetAK8subjet3pt[ijet];
                        pfjetAK8subjet3mass[njet] = pfjetAK8subjet3mass[ijet];


                        ak8tau32->Fill(pfjetAK8tau3[njet]/pfjetAK8tau2[njet]);


                        if (((pfjetAK8tau3[njet]/pfjetAK8tau2[njet]) < 0.54) && (pfjetAK8sdmass[njet] > 105 && pfjetAK8sdmass[njet] < 220) && (pfjetAK8sub1btag[njet] > 0.4184 || pfjetAK8sub2btag[njet] > 0.4184))
                        {
                                cut_flow_leading_w_or_b->Fill(0);
                          
                                if (pfjetAK8sub1mass[njet]>=65 && pfjetAK8sub1mass[njet]<=95)
                                {
                                        cut_flow_leading_w_or_b->Fill(1);
                                }
                                else
                                {
                                        cut_flow_leading_w_or_b->Fill(2);
                                }

                                ak8tau32withcut->Fill(pfjetAK8tau3[njet]/pfjetAK8tau2[njet]);    

				TLorentzVector pfjetAK8vec;
                                pfjetAK8vec.SetPtEtaPhiM(pfjetAK8pt[njet], pfjetAK8eta[njet], pfjetAK8phi[njet], pfjetAK8mass[njet]);
                                

                                subwpt->Fill(pfjetAK8sub1pt[njet]);
                                subbpt->Fill(pfjetAK8sub2pt[njet]);
                                ptfrac1->Fill(pfjetAK8sub1pt[njet]/pfjetAK8pt[njet]);
                                submass1->Fill(pfjetAK8sub1mass[njet]);
                                ptfrac2->Fill(pfjetAK8sub2pt[njet]/pfjetAK8pt[njet]);
                                submass2->Fill(pfjetAK8sub2mass[njet]);
                                sub1ptfrac->Fill(pfjetAK8subjet1pt[njet]/pfjetAK8pt[njet]);
                                sub2ptfrac->Fill(pfjetAK8subjet2pt[njet]/pfjetAK8pt[njet]);
                                sub3ptfrac->Fill(pfjetAK8sub2pt[njet]/pfjetAK8pt[njet]);
                                subq1pt->Fill(pfjetAK8subjet1pt[njet]);
                                subq2pt->Fill(pfjetAK8subjet2pt[njet]);
                                sub1mass->Fill(pfjetAK8subjet1mass[njet]);
                                sub2mass->Fill(pfjetAK8subjet2mass[njet]);
                                sub3mass->Fill(pfjetAK8sub2mass[njet]);

                                TLorentzVector topsub1, topsub2, topsub3, topvec, wvec, Mbj1, Mbj2;

                                topsub1.SetPtEtaPhiM(pfjetAK8subjet1pt[njet], pfjetAK8subjet1eta[njet], pfjetAK8subjet1phi[njet], pfjetAK8subjet1mass[njet]);
                                topsub2.SetPtEtaPhiM(pfjetAK8subjet2pt[njet], pfjetAK8subjet2eta[njet], pfjetAK8subjet2phi[njet], pfjetAK8subjet2mass[njet]);
                                topsub3.SetPtEtaPhiM(pfjetAK8sub2pt[njet], pfjetAK8sub2eta[njet], pfjetAK8sub2phi[njet], pfjetAK8sub2mass[njet]);

                                wvec=topsub1+topsub2;
                                topvec=wvec+topsub3;
				Mbj1=topsub1+topsub3;
				Mbj2=topsub2+topsub3;

                                recon_top_mass->Fill(topvec.M());
                                recon_w_mass->Fill(wvec.M());
				bj1_invmass->Fill(Mbj1.M());
				bj2_invmass->Fill(Mbj2.M());

                                sub1efrac->Fill(topsub1.E()/pfjetAK8vec.E());
                                sub2efrac->Fill(topsub2.E()/pfjetAK8vec.E());
                                sub3efrac->Fill(topsub3.E()/pfjetAK8vec.E());

                                q1q2_dphi->Fill(abs(topsub1.DeltaPhi(topsub2)));

                                if (pfjetAK8mass[njet] >= 150 && pfjetAK8mass[njet] <= 200)
                                {
                                         sub1efrac_topmasscut->Fill(topsub1.E()/pfjetAK8vec.E());
                                         sub2efrac_topmasscut->Fill(topsub2.E()/pfjetAK8vec.E());
                                         sub3efrac_topmasscut->Fill(topsub3.E()/pfjetAK8vec.E());
                                         q1q2_dphi_topmasscut->Fill(abs(topsub1.DeltaPhi(topsub2)));
                                         reco_top_mass_ratio->Fill(1-topvec.M()/pfjetAK8vec.M());
                                }

                                double d12 = min(pow(pfjetAK8subjet1pt[njet], 2), pow(pfjetAK8subjet2pt[njet], 2))*pow(delta2R(pfjetAK8subjet1eta[njet], pfjetAK8subjet1phi[njet], pfjetAK8subjet2eta[njet], pfjetAK8subjet2phi[njet]), 2);

                                double d23 = min(pow(pfjetAK8subjet2pt[njet], 2), pow(pfjetAK8sub2pt[njet], 2))*pow(delta2R(pfjetAK8subjet2eta[njet], pfjetAK8subjet2phi[njet], pfjetAK8sub2eta[njet], pfjetAK8sub2phi[njet]), 2);

                                double d13 = min(pow(pfjetAK8subjet1pt[njet], 2), pow(pfjetAK8sub2pt[njet], 2))*pow(delta2R(pfjetAK8subjet1eta[njet], pfjetAK8subjet1phi[njet], pfjetAK8sub2eta[njet], pfjetAK8sub2phi[njet]), 2);

                                if ((d12 < d23) && (d12 < d13))
                                {
                                        double zk12 = max(topsub1.E(), topsub2.E())/pfjetAK8vec.E();
                                        zkdist->Fill(zk12);
                                }

                                if ((d23 < d12) && (d23 < d13))
                                {
                                        double zk23 = max(topsub2.E(), topsub3.E())/pfjetAK8vec.E();
                                        zkdist->Fill(zk23);
                                }

                                if ((d13 < d23) && (d13 < d12))
                                {
                                        double zk13 = max(topsub1.E(), topsub3.E())/pfjetAK8vec.E();
                                        zkdist->Fill(zk13);
                                }

                                double z12 = min(pfjetAK8subjet1pt[njet], pfjetAK8subjet2pt[njet])/(pfjetAK8subjet1pt[njet] + pfjetAK8subjet2pt[njet]);
                                double delR12 = delta2R(pfjetAK8subjet1eta[njet], pfjetAK8subjet1phi[njet], pfjetAK8subjet2eta[njet], pfjetAK8subjet2phi[njet])/0.8;
                                z12dist->Fill(z12);
                                delR12dist->Fill(delR12);
                      //          z12vsdelR12->Fill(delR12,z12);

				cutflow_j1b_j2b->Fill(0);

				if (Mbj1.M() < Mbj2.M())
				{
				        efrac_j1_j2->Fill(topsub1.E()/pfjetAK8vec.E());
					cutflow_j1b_j2b->Fill(1);
				}

				else
				{
					efrac_j1_j2->Fill(topsub2.E()/pfjetAK8vec.E());
					cutflow_j1b_j2b->Fill(2);
				}   


                                if (Mbj1.M() > Mbj2.M())
                                {
                                        efrac_j1j2->Fill(topsub1.E()/pfjetAK8vec.E());
                                }

                                else
                                {
                                        efrac_j1j2->Fill(topsub2.E()/pfjetAK8vec.E());
                                }    	

                        }

                        if (((pfjetAK8tau3[njet]/pfjetAK8tau2[njet]) < 0.54) && (pfjetAK8sdmass[njet] > 105 && pfjetAK8sdmass[njet] < 220) && (pfjetAK8sub1btag[njet] > 0.4184))
                        {
                                w_mass->Fill(pfjetAK8sub2mass[njet]);
                        }
        
                        if (((pfjetAK8tau3[njet]/pfjetAK8tau2[njet]) < 0.54) && (pfjetAK8sdmass[njet] > 105 && pfjetAK8sdmass[njet] < 220) && (pfjetAK8sub2btag[njet] > 0.4184))
                        {
                                w_mass->Fill(pfjetAK8sub1mass[njet]);
                        }


                        if(++njet>=njetmxAK8) break;
             }
	   
             npfjetAK8 = njet;
 
             for (int ijet=0; ijet<npfjetAK8; ijet++)
             {
                        for(int jjet=ijet+1; jjet<npfjetAK8; jjet++)
                        {
                                  if(pfjetAK8pt[ijet]<pfjetAK8pt[jjet])
                                  {
                                           float b = pfjetAK8pt[ijet];
                                           pfjetAK8pt[ijet] = pfjetAK8pt[jjet];
                                           pfjetAK8pt[jjet] = b;

                                           b = pfjetAK8eta[ijet];
                                           pfjetAK8eta[ijet] = pfjetAK8eta[jjet];
                                           pfjetAK8eta[jjet] = b;

                                           b = pfjetAK8phi[ijet];
                                           pfjetAK8phi[ijet] = pfjetAK8phi[jjet];
                                           pfjetAK8phi[jjet] = b;

                                           b=pfjetAK8mass[ijet];
                                           pfjetAK8mass[ijet]=pfjetAK8mass[jjet];
                                           pfjetAK8mass[jjet]=b;

                                           b=pfjetAK8tau3[ijet];
                                           pfjetAK8tau3[ijet]=pfjetAK8tau3[jjet];
                                           pfjetAK8tau3[jjet]=b;

                                           b=pfjetAK8sdmass[ijet];
                                           pfjetAK8sdmass[ijet]=pfjetAK8sdmass[jjet];
                                           pfjetAK8sdmass[jjet]=b;

                                           b=pfjetAK8tau2[ijet];
                                           pfjetAK8tau2[ijet]=pfjetAK8tau2[jjet];
                                           pfjetAK8tau2[jjet]=b;

                                           b=pfjetAK8sub1btag[ijet];
                                           pfjetAK8sub1btag[ijet]=pfjetAK8sub1btag[jjet];
                                           pfjetAK8sub1btag[jjet]=b;

                                           b=pfjetAK8sub2btag[ijet];
                                           pfjetAK8sub2btag[ijet]=pfjetAK8sub2btag[jjet];
                                           pfjetAK8sub2btag[jjet]=b;

                                           b=pfjetAK8sub1pt[ijet];
                                           pfjetAK8sub1pt[ijet]=pfjetAK8sub1pt[jjet];
                                           pfjetAK8sub1pt[jjet]=b;
 
                                           b=pfjetAK8sub2pt[ijet];
                                           pfjetAK8sub2pt[ijet]=pfjetAK8sub2pt[jjet];
                                           pfjetAK8sub2pt[jjet]=b;

                                           b=pfjetAK8sub2mass[ijet];
                                           pfjetAK8sub2mass[ijet]=pfjetAK8sub2mass[jjet];
                                           pfjetAK8sub2mass[jjet]=b;

                                           b=pfjetAK8sub1mass[ijet];
                                           pfjetAK8sub1mass[ijet]=pfjetAK8sub1mass[jjet];
                                           pfjetAK8sub1mass[jjet]=b;

                                           b=pfjetAK8subjet1pt[ijet];
                                           pfjetAK8subjet1pt[ijet]=pfjetAK8subjet1pt[jjet];
                                           pfjetAK8subjet1pt[jjet]=b;

                                           b=pfjetAK8subjet1mass[ijet];
                                           pfjetAK8subjet1mass[ijet]=pfjetAK8subjet1mass[jjet];
                                           pfjetAK8subjet1mass[jjet]=b;

                                           b=pfjetAK8subjet2pt[ijet];
                                           pfjetAK8subjet2pt[ijet]=pfjetAK8subjet2pt[jjet];
                                           pfjetAK8subjet2pt[jjet]=b;

                                           b=pfjetAK8subjet2mass[ijet];
                                           pfjetAK8subjet2mass[ijet]=pfjetAK8subjet2mass[jjet];
                                           pfjetAK8subjet2mass[jjet]=b;

                                           b=pfjetAK8subjet3pt[ijet];
                                           pfjetAK8subjet3pt[ijet]=pfjetAK8subjet3pt[jjet];
                                           pfjetAK8subjet3pt[jjet]=b;

                                           b=pfjetAK8subjet3mass[ijet];
                                           pfjetAK8subjet3mass[ijet]=pfjetAK8subjet3mass[jjet];
                                           pfjetAK8subjet3mass[jjet]=b;
                                   }
                         }
             }


	    TLorentzVector ak8jet1, ak8jet2, res1;
	    ak8jet1.SetPtEtaPhiM(pfjetAK8pt[0], pfjetAK8eta[0], pfjetAK8phi[0], pfjetAK8mass[0]);
            ak8jet2.SetPtEtaPhiM(pfjetAK8pt[1], pfjetAK8eta[1], pfjetAK8phi[1], pfjetAK8mass[1]);
            res1=ak8jet1+ak8jet2;


            double leadTopJetPt=0, leadTopJetEta=-100, leadTopJetPhi=-100, leadTopJetMass=-100, leadTopJetRap=-100;

        
            for(int ijet=0; ijet<topjetCol.size(); ijet++)
            {
                        if(leadTopJetPt<topjetCol[ijet].Pt())
                        {
                                leadTopJetPt=topjetCol[ijet].Pt();
                                leadTopJetEta=topjetCol[ijet].Eta();
                                leadTopJetRap=topjetCol[ijet].Rapidity();
				leadTopJetPhi=topjetCol[ijet].Phi();
                                leadTopJetMass=topjetCol[ijet].M();
                        }
                        
            }

           
            njet = 0;

            for(int ijet=0; ijet<npfjetAK8; ijet++)
            {

                        if(!pfjetAK8tightID[ijet]) continue;
                        pfjetAK8pt[ijet] *= pfjetAK8JEC[ijet];
                        pfjetAK8mass[ijet] *= pfjetAK8JEC[ijet];
                        pfjetAK8pt[ijet] *= (1.+pfjetAK8reso[ijet]);
                        if(fabs(pfjetAK8eta[ijet]) > 2.5) continue;
                        if(pfjetAK8pt[ijet]<20) continue;

                        pfjetAK8pt[njet]  = pfjetAK8pt[ijet];
                        pfjetAK8eta[njet] = pfjetAK8eta[ijet];
                        pfjetAK8phi[njet]  = pfjetAK8phi[ijet];
                        pfjetAK8mass[njet]  = pfjetAK8mass[ijet];
                        pfjetAK8tau3[njet] = pfjetAK8tau3[ijet];
                        pfjetAK8tau2[njet] = pfjetAK8tau2[ijet];
                        pfjetAK8sdmass[njet] = pfjetAK8sdmass[ijet];
                        pfjetAK8sub1btag[njet] = pfjetAK8sub1btag[ijet];
                        pfjetAK8sub2btag[njet] = pfjetAK8sub2btag[ijet];
                        pfjetAK8sub1pt[njet] = pfjetAK8sub1pt[ijet];
                        pfjetAK8sub2pt[njet] = pfjetAK8sub2pt[ijet];
                        pfjetAK8sub2mass[njet] = pfjetAK8sub2mass[ijet];
                        pfjetAK8sub1mass[njet] = pfjetAK8sub1mass[ijet];
                        pfjetAK8subjet1pt[njet] = pfjetAK8subjet1pt[ijet];
                        pfjetAK8subjet1mass[njet] = pfjetAK8subjet1mass[ijet];
                        pfjetAK8subjet2pt[njet] = pfjetAK8subjet2pt[ijet];
                        pfjetAK8subjet2mass[njet] = pfjetAK8subjet2mass[ijet];
                        pfjetAK8subjet3pt[njet] = pfjetAK8subjet3pt[ijet];
                        pfjetAK8subjet3mass[njet] = pfjetAK8subjet3mass[ijet];

                        if (delta2R(leadTopJetEta, leadTopJetPhi, pfjetAK8eta[njet], pfjetAK8phi[njet]) < 0.4)
                        {
                                  ak8mass_mat_top->Fill(pfjetAK8mass[njet]);
                        }

                        if(++njet>=njetmxAK8) break;
            }


            double leadsdTopJetPt=0, leadsdTopJetEta=-100, leadsdTopJetPhi=-100, leadsdTopJetMass=-100;


            for(int ijet=0; ijet<sdtopjetCol.size(); ijet++)
            {
                        if(leadsdTopJetPt<sdtopjetCol[ijet].Pt())
                        {
                                leadsdTopJetPt=sdtopjetCol[ijet].Pt();
                                leadsdTopJetEta=sdtopjetCol[ijet].Eta();
                                leadsdTopJetPhi=sdtopjetCol[ijet].Phi();
                                leadsdTopJetMass=sdtopjetCol[ijet].M();
                        }
            }

            double leadsdtauTopJetPt=0, leadsdtauTopJetEta=-100, leadsdtauTopJetPhi=-100, leadsdtauTopJetMass=-100;


            for(int ijet=0; ijet<sdtautopjetCol.size(); ijet++)
            {
                        if(leadsdtauTopJetPt<sdtautopjetCol[ijet].Pt())
                        {
                                leadsdtauTopJetPt=sdtautopjetCol[ijet].Pt();
                                leadsdtauTopJetEta=sdtautopjetCol[ijet].Eta();
                                leadsdtauTopJetPhi=sdtautopjetCol[ijet].Phi();
                                leadsdtauTopJetMass=sdtautopjetCol[ijet].M();
                        }
            }

	   
            if (topjetCol.size()==1)
            {
                    cut_flow_top->Fill(0); 

	    if ((delta2R(leadTopJetEta, leadTopJetPhi, leadGenTopEta, leadGenTopPhi) < 0.4) || (delta2R(leadTopJetEta, leadTopJetPhi, leadGenAntiTopEta, leadGenAntiTopPhi) < 0.4))
	    {
		    cut_flow_top->Fill(1);
	    }
	    else
	    {
		    cut_flow_top->Fill(2);
	    }
            }

    
            top_pt_check->Fill(0);
            if (leadTopJetPt >= 300 && leadTopJetPt <= 500)
            {
                    top_pt_check->Fill(1);
            }
            if (leadTopJetPt > 500 && leadTopJetPt <= 800)
            {
                    top_pt_check->Fill(2);
            }
            if (leadTopJetPt > 800)
            {
                    top_pt_check->Fill(3);
            }



	    bjet_pt->Fill(leadbJetPt);
            bjet_eta->Fill(leadbJetEta);
	    bjpt = leadbJetPt;
            bjeta = leadbJetEta;
	    nbjet->Fill(bjetCol.size());
            bj = bjetCol.size();
	    Delphi = delphi;
  

            pfjetAK4_no->Fill(alljetCol.size());
            ak4j = alljetCol.size();
            pfjetAK4_pt->Fill(leadallJetPt);
            pfjetAK4_eta->Fill(leadallJetEta);
            pfjetAK4_mass->Fill(leadallJetMass);
            ak4jpt = leadallJetPt;
            ak4jeta = leadallJetEta;
           

            pfjetAK8_no->Fill(AK8jetCol.size());
            ak8j = AK8jetCol.size();
            pfjetAK8_pt->Fill(leadAK8JetPt);
            pfjetAK8_eta->Fill(leadAK8JetEta);
            pfjetAK8_rap->Fill(leadAK8JetRap);
            pfjetAK8_mass->Fill(leadAK8JetMass);
            ak8jpt = leadAK8JetPt;
            ak8jeta = leadAK8JetEta;
	    ak8invmass->Fill(res1.M());


            topjet_no->Fill(topjetCol.size());
            topjet_pt->Fill(leadTopJetPt);
            topjet_eta->Fill(leadTopJetEta);
            topjet_rap->Fill(leadTopJetRap);
            topjet_mass->Fill(leadTopJetMass);
            topjet_sdmass->Fill(leadsdTopJetMass);
            topjet_sdtaumass->Fill(leadsdtauTopJetMass);

	    el_no->Fill(elecCol.size());
	    el_pt->Fill(leadElePt);
	    el_eta->Fill(leadEleEta);

	    nele = elecCol.size();
	    nelept = leadElePt;
	    neleeta = leadEleEta;

	    pfMET->Fill(miset);
	    trmass_rec->Fill(mt_rec);

	    mEt = miset;
	    TrMass = mt_rec;

	

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

//cout<<count1<<"		"<<count2<<"	"<<count3<<"	"<<count4<<"	"<<count5<<"	"<<count6<<endl;
	return 0;
}
