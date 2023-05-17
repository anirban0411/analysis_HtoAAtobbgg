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


void sig_find_delphi_1D_v2()
{

	clock_t start;
        clock_t end;

        start = clock();

	TFile* f = TFile::Open("combined_2018_wh_v1_bin50.root");

	TH1F *h1 = (TH1F*) f->Get("sig_20");
        TH1F *h2 = (TH1F*) f->Get("data_obs");
        TH1F *h3 = (TH1F*) f->Get("sig_55");

	int dim = h1->GetNbinsX();
	float significance20[dim], significance55[dim];

	for(int i = 1; i < dim; i++)
	{
		float sum20 = 0.0;
        	float sum55 = 0.0;
		float bkg = 0.0;

		for(int j = i+1; j <= dim; j++)
		{
			float s20 = h1->GetBinContent(j);
                	float b = h2->GetBinContent(j);
                	float s55 = h3->GetBinContent(j);

			sum20 = sum20 + s20;
			sum55 = sum55 + s55;
			bkg = bkg + b;
		}

		if ((bkg) <= 0) break;
		if ((bkg) <= 0) break;

		significance20[i] = sum20/sqrt(bkg);
		significance55[i] = sum55/sqrt(bkg);

		cout << sum20 << "	" << bkg << "	" << significance20[i]  << "	" << (M_PI/dim)*i << endl;
	}

	float max20 = -999.0;
        int index20 = -1;
        float max55 = -999.0;
        int index55 = -1;

        for (int i = 1; i < dim; i++)
        {
                if(significance20[i]>max20)
                {
                        max20 = significance20[i];
                        index20 = i;
                }

                if(significance55[i]>max55)
                {
                        max55 = significance55[i];
                        index55 = i;
                }
        }


        cout << "max sensitivity at 20 GeV: " << max20 << endl;
	cout << "corresponding delphi at 20 GeV: " << (M_PI/dim)*index20 << endl;
        cout << "corresponding bin no at 20 GeV: " << index20 << endl;
        cout << "max sensitivity at 55 GeV: " << max55 << endl;
	cout << "corresponding delphi at 55 GeV: " << (M_PI/dim)*index55 << endl;
        cout << "corresponding bin no at 55 GeV: " << index55 << endl;


	end = clock();

        timetaken = (end - start) / (double)CLOCKS_PER_SEC;

        cout << "Time taken : " << timetaken << "s" << endl;

}
