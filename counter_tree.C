#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

using namespace std;

void counter_tree(string filename)
{

int nevent_total;
double tot_weight;

char rootfiles[100];
char infile[1000];
char datafile[1000];

//cout <<"Give the input file name"<<endl;
//cin>> rootfiles;
sprintf(rootfiles,"%s",filename.c_str());

ifstream file_db;
file_db.open(rootfiles);

double count[2]={0};

while(!(file_db.eof())){

	file_db >> datafile;
	if(file_db.eof()) break;

	sprintf(infile, "%s", datafile);
	cout<<"infile "<<infile<<endl;

        TFile* fileIn = TFile::Open(infile);
	double sumweight = 0;
	double sum_pos =0, sum_neg = 0;

	double event_weight_LHE;

	TTree* T1;
	T1 = (TTree*)fileIn->Get("Events_All");
        //                     T1->SetBranchAddress("event_weight", &event_weight_LHE);
	T1->SetBranchAddress("Generator_weight", &event_weight_LHE);

	int nevt = T1->GetEntries();

	for(int iev=0; iev<nevt; iev++){
		T1->GetEntry(iev);
		count[0]++;
		count[1] += event_weight_LHE;
        //	cout<<iev<<" weight "<<event_weight_LHE<<endl;
		//if((iev+1)%100000 == 1) { cout<<iev<<" weight "<<event_weight_LHE<<" event sum "<<count[1]<<" weight sum "<<count[2]<<endl; }

		//cout<<"sum "<<sumweight<<endl;
	}

	fileIn->cd();
	delete T1;
	delete fileIn;
}

cout<<" number of entries "<<count[0]<<" sumofweights "<<count[1]<<endl;

}
/*
int main()
{
counter_tree();
}
*/
