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



void file_out_1()
{

	TFile* f = TFile::Open("combined_2018_wh_v1_cat15.root");

	TH1F *h1 = (TH1F*) f->Get("sig_20");
	TH1F *h2 = (TH1F*) f->Get("data_obs");
	TH1F *h3 = (TH1F*) f->Get("sig_55");

	float sum20 = 0.0;
	float sum55 = 0.0;

	for(int i = 1; i <= h1->GetNbinsX(); i++)
	{
		float s20 = h1->GetBinContent(i);
		float b = h2->GetBinContent(i);
		float s55 = h3->GetBinContent(i);

		if (s20+b <= 0) continue;
		if (s55+b <= 0) continue;

		sum20 = sum20 + (s20/sqrt(s20+b))*(s20/sqrt(s20+b));
		sum55 = sum55 + (s55/sqrt(s55+b))*(s55/sqrt(s55+b));
	}

	cout << "significance at 20 GeV = " << sqrt(sum20) << endl;
	cout << "significance at 55 GeV = " << sqrt(sum55) << endl;


}
