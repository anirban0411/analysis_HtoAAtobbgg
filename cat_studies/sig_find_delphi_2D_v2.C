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
#include <iomanip>



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


void sig_find_delphi_2D_v2()
{

	clock_t start;
        clock_t end;

        start = clock();

	TFile* f = TFile::Open("combined_2018_wh_v1_bin20.root");

	TH1F *h1 = (TH1F*) f->Get("sig_20");
        TH1F *h2 = (TH1F*) f->Get("data_obs");
        TH1F *h3 = (TH1F*) f->Get("sig_55");

	int dim = h1->GetNbinsX();
	float significance20[dim], significance55[dim], significance20_1[dim*dim], significance55_1[dim*dim], sig_sum_20[dim*dim], sig_sum_55[dim*dim];


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

		if ((sum20+bkg) <= 0) continue;
		if ((sum55+bkg) <= 0) continue;

		significance20[i] = sum20/sqrt(sum20+bkg);
		significance55[i] = sum55/sqrt(sum55+bkg);

		for (int k = 1; k < i; k++)
		{
			float sum20_1 = 0.0;
                	float sum55_1 = 0.0;
                	float bkg_1 = 0.0;

			for(int l = k+1; l <= i; l++)
			{
				float s20_1 = h1->GetBinContent(l);
                        	float b_1 = h2->GetBinContent(l);
                        	float s55_1 = h3->GetBinContent(l);

				sum20_1 = sum20_1 + s20_1;
                        	sum55_1 = sum55_1 + s55_1;
                        	bkg_1 = bkg_1 + b_1;
			}

			if ((sum20_1+bkg_1) <= 0) continue;
                	if ((sum55_1+bkg_1) <= 0) continue;

                	significance20_1[k] = sum20_1/sqrt(sum20_1+bkg_1);
                	significance55_1[k] = sum55_1/sqrt(sum55_1+bkg_1);
		}


	}


	int indx = 1;
	for (int i = 1; i < dim; i++)
        {
		for (int k = 1; k < i; k++)
		{
			sig_sum_20[indx] = sqrt(significance20_1[k]*significance20_1[k] + significance20[i]*significance20[i]);
			sig_sum_55[indx] = sqrt(significance55_1[k]*significance55_1[k] + significance55[i]*significance55[i]);
			indx++;
		}
        }

	float max20 = -999.0;
        float max55 = -999.0;

        for (int i = 1; i <= indx; i++)
        {
                if(sig_sum_20[i]>max20)
                {
                        max20 = sig_sum_20[i];
                }

                if(sig_sum_55[i]>max55)
                {
                        max55 = sig_sum_55[i];
                }
        }


	cout<<setprecision(5)<<fixed;
	for (int i = 1; i < dim; i++)
        {
                for (int k = 1; k < i; k++)
                {
                        float sig_sum20 = sqrt(significance20_1[k]*significance20_1[k] + significance20[i]*significance20[i]);
                        float sig_sum55 = sqrt(significance55_1[k]*significance55_1[k] + significance55[i]*significance55[i]);

			cout << "sensitivity at 20 GeV: " << sig_sum20 << "	" << "sensitivity at 55 GeV: " << sig_sum55 << "	" << "lower bin: " << k << "	" << "higher bin: " << i << endl;
                }
        }


	int a0 = 0; 
	int b0 = 0; 
	int a1 = 0; 
	int b1 = 0;

	for (int i = 1; i < dim; i++)
        {
                for (int k = 1; k < i; k++)
                {
                        float sig_sum20 = sqrt(significance20_1[k]*significance20_1[k] + significance20[i]*significance20[i]);
			if (sig_sum20 == max20)
			{
				a0 = k;
				b0 = i;
				break;
			}

			float sig_sum55 = sqrt(significance55_1[k]*significance55_1[k] + significance55[i]*significance55[i]);
                        if (sig_sum55 == max55)
                        {
				a1 = k;
				b1 = i;
				break;
                        }
                }
        }


	cout << "\n";
	cout << "corresponding bins no at 20 GeV: " << a0 << "   " << b0 << endl;
	cout << "corresponding delphi at 20 GeV: " << (M_PI/dim)*a0 << "        " << (M_PI/dim)*b0 << endl;
	cout << "corresponding bins no at 55 GeV: " << a1 << "   " << b1 << endl;
	cout << "corresponding delphi at 55 GeV: " << (M_PI/dim)*a1 << "	" << (M_PI/dim)*b1 << endl << "\n";
        cout << "max sensitivity at 20 GeV: " << max20 << endl;
        cout << "max sensitivity at 55 GeV: " << max55 << endl;


	end = clock();

        timetaken = (end - start) / (double)CLOCKS_PER_SEC;

        cout << "Time taken : " << timetaken << "s" << endl;

}
