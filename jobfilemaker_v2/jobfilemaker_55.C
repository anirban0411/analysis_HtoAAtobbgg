#include<stdio.h>
#include <bits/stdc++.h>
using namespace std;

int main()
{
fstream file3;
char name_buffer3[512];
file3.open("condor_file/condor_WH_mA_55.sh",ios::out);

for(int i=1; i <= 9; i =i+2)
{
        int j = i + 2;
        fstream file;
        char name_buffer[512];
        sprintf(name_buffer,"exe_file/exe_WH_mA_55_%d_%d.sh",i,j);
        file.open(name_buffer,ios::out);
   if(!file)
   {
       cout<<"Error in creating file!!!";
       return 0;
   }

   file << "#!/bin/bash\ncd /home/abala/cms/CMSSW_10_5_0/src/condor_job/ \nexport VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\nexport SCRAM_ARCH=slc6_amd64_gcc630\nsource /cvmfs/cms.cern.ch/cmsset_default.sh\n#cmsenv\nexport X509_USER_PROXY=/home/abala/x509up_u56615\neval\n./analysis  "<< i << " " << j << "   all_log_files/WHmA55.log\n" ;
   cout<<"exe created successfully." << endl;
   file.close();

         fstream file1;
         char name_buffer1[512];
         sprintf(name_buffer1,"sub_file/sub_WH_mA_55_%d_%d.sh",i,j);
         file1.open(name_buffer1,ios::out);
   if(!file1)
   {
       cout<<"Error in creating file!!!";
       return 0;
   }
         file1 << "universe = vanilla\nexecutable = exe_file/exe_WH_mA_55_" << i <<"_"<< j << ".sh\ngetenv = TRUE\nlog =/home/abala/t3store3/Higgs/others/WH_mA_55_condor_2018_" << i <<"To"<<j<<".log\noutput =/home/abala/t3store3/Higgs/others/WH_mA_55_condor_2018_" << i <<"To"<<j<<".out\nerror =/home/abala/t3store3/Higgs/others/WH_mA_55_condor_2018_" << i <<"To"<<j<<".error\nnotification = never\nshould_transfer_files = YES\nwhen_to_transfer_output = ON_EXIT\nqueue";
   cout<<"sub created successfully." << endl;
   file1.close();


   file3 << "condor_submit " << name_buffer1 << endl;
    
}
   return 0;
}

