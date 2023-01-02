#!/bin/bash

OPT=$1
FILE=$2
GPP_FLAGS_EXTRA=$3

NAME="${FILE##*/}"
EXT="${FILE##*.}"
NAME="${FILE%.*}"

echo "Full file name:" $FILE
echo "File name:" $NAME
echo "File extension:" $EXT
echo "Extra g++ flags:" $GPP_FLAGS_EXTRA


ROOT_FLAGS=`root-config --libs --cflags --ldflags`" -lRooFitCore -lMinuit -lTMVA -lTMVAGui"
#CLHEP_FLAGS=`clhep-config --include --ldflags --libs`
#PYTHIA8_FLAGS=`pythia8-config --cflags --ldflags`
#HEPMC_FLAGS="-I/opt/HepMC-2.06.09/include -L/opt/HepMC-2.06.09/lib -lHepMC -Wl,-rpath=/opt/HepMC-2.06.09/lib"
#LHAPDF_FLAGS=`lhapdf-config --cflags --libs`
#FASTJET_FLAGS=`fastjet-config --cxxflags --libs`

#DELPHES_PATH="/opt/Delphes-3.4.1"
#DELPHES_PATH="/opt/delphes"
#DELPHES_FLAGS="-I$DELPHES_PATH -L$DELPHES_PATH -lDelphes -Wl,-rpath=$DELPHES_PATH"

#HEPTOPTAGGER_PATH="/opt/HEPTopTagger2"
#HEPTOPTAGGER_FLAGS="-L$HEPTOPTAGGER_PATH/Nsubjettiness -lNsubjettiness -L$HEPTOPTAGGER_PATH/qjets/lib -lQjets -I/opt/HEPTopTagger2"
#HEPTOPTAGGER_FILES="$HEPTOPTAGGER_PATH/HEPTopTagger.cc"

#CUSTOM_HEADERS=`find /media/soham/D/Programs/TIFR/PhD/misc/ | grep ".cc"`" "`find /media/soham/D/Programs/TIFR/PhD/analysis/CMSSW_8_0_20/src/stopPair/HeaderFiles/ | grep ".cc"`

#CUSTOM_HEADER_PATH="/media/soham/D/Programs/TIFR/PhD/analysis/CMSSW_8_0_20/src/stopPair"
CUSTOM_HEADER_PATH="/home/abala/cms/CMSSW_10_5_0/src/condor_job"
#CUSTOM_HEADERS=`find $CUSTOM_HEADER_PATH/HeaderFiles/ | grep ".cc$"`

#g++ -std=c++11 $FILE.cc $CUSTOM_HEADERS -o $FILE $ROOT_FLAGS $CLHEP_FLAGS $PYTHIA8_FLAGS $FASTJET_FLAGS $HEPMC_FLAGS -ldl


GPP_FLAGS="-std=c++11 -O2 -w $GPP_FLAGS_EXTRA"

if [ "$OPT" == "n" ]; then
    g++ $GPP_FLAGS $NAME"."$EXT -o $NAME
elif [ "$OPT" == "h" ]; then
    g++ $GPP_FLAGS $HEPTOPTAGGER_FILES $NAME"."$EXT -o $NAME -ldl
elif [ "$OPT" == "ch" ]; then
    g++ $GPP_FLAGS $NAME"."$EXT $CUSTOM_HEADERS -o $NAME $ROOT_FLAGS -ldl -I$CUSTOM_HEADER_PATH
else
    echo "Wrong g++mod option."
fi



#if [ "$OPT" == "n" ]; then
#    g++ $GPP_FLAGS $NAME"."$EXT -o $NAME
#elif [ "$OPT" == "h" ]; then
#    g++ $GPP_FLAGS $HEPTOPTAGGER_FILES $NAME"."$EXT -o $NAME $ROOT_FLAGS $CLHEP_FLAGS $PYTHIA8_FLAGS $FASTJET_FLAGS $HEPMC_FLAGS $DELPHES_FLAGS $HEPTOPTAGGER_FLAGS -ldl
#elif [ "$OPT" == "ch" ]; then
#    g++ $GPP_FLAGS $NAME"."$EXT $CUSTOM_HEADERS -o $NAME $ROOT_FLAGS $CLHEP_FLAGS $PYTHIA8_FLAGS $FASTJET_FLAGS $HEPMC_FLAGS $DELPHES_FLAGS $HEPTOPTAGGER_FLAGS -ldl -I$CUSTOM_HEADER_PATH
#else
#    echo "Wrong g++mod option."
#fi
