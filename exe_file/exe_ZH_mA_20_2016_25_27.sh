#!/bin/bash
cd /home/abala/cms/CMSSW_10_5_0/src/condor_job/ 
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
export SCRAM_ARCH=slc6_amd64_gcc630
source /cvmfs/cms.cern.ch/cmsset_default.sh
#cmsenv
export X509_USER_PROXY=/home/abala/x509up_u56615
eval
./analysis  25 27   new_v3_2016/ZHmA20_2016PostAPV_new_v1.log
