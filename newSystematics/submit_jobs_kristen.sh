#!/bin/bash

modelName=$1

export BATCH_DEBUG=true
echo $BATCH_DEBUG
echo $modelName 

if [[ $modelName == "" ]]
then exit
fi

#DIR=$(pwd)
DIR=/afs/cern.ch/user/w/winstrom/RA2b_limit
cd $DIR
echo $DIR

rm -f submit_scripts/$modelName*kristen_sub.sh

for cut in ge1bTight ge1bLoose ge2bTight ge2bLoose ge3bLoose
do
    modelFile='/afs/cern.ch/user/j/joshmt/public/RA2bFall2011/signalSyst.'$modelName'.'$cut'.dat'
    echo $modelFile
    declare -i nLines=0
    declare -i iLine=0
    #cat $modelFile | while read m0 m12 rest
    cat reducedModel.T1bbbb.txt | while read m0 m12
    do 
	echo '#! /bin/bash -f' > 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_frequentist_sub.sh'
	echo 'export BATCH_DEBUG=true' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_frequentist_sub.sh'
	echo 'if [ ! -d "/tmp/ra2b" ]; then mkdir /tmp/ra2b; fi;' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_frequentist_sub.sh'
	#echo 'if [ -d "/tmp/ra2b" ]; then echo "found directory!"; fi;' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_frequentist_sub.sh'
	echo 'cd '$DIR >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_frequentist_sub.sh'
	echo 'source /afs/cern.ch/sw/lcg/contrib/gcc/4.3/x86_64-slc5/setup.sh' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_frequentist_sub.sh'
	echo 'source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.33.02/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_frequentist_sub.sh'
	echo '' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_frequentist_sub.sh'
	echo 'root -q -l -b runKristenLimit.C'"'"'("'$cut'","'$modelName'",'$m0','$m12',0,0)'"'" >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_frequentist_sub.sh'
	echo 'if [ -e "/tmp/ra2b/ws_expected_frequentist_'$cut'_'$modelName'_'$m0'_'$m12'_kristen.root" ]; then rm /tmp/ra2b/ws_expected_frequentist_'$cut'_'$modelName'_'$m0'_'$m12'_kristen.root; fi' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_frequentist_sub.sh'
	chmod a+x 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_frequentist_sub.sh'
	echo "Submitting job "$modelName'_'$cut'_'$m0'_'$m12"_expected_kristen_frequentist to batch system for processing."
	#bsub -q 1nd -N -o /dev/null -J $modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_frequentist' -Jd $modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_frequentist' $DIR'/submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_frequentist_sub.sh'


	echo '#! /bin/bash -f' > 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_frequentist_sub.sh'
	echo 'export BATCH_DEBUG=true' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_frequentist_sub.sh'
	echo 'if [ ! -d "/tmp/ra2b" ]; then mkdir /tmp/ra2b; fi;' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_frequentist_sub.sh'
	#echo 'if [ -d "/tmp/ra2b" ]; then echo "found directory!"; fi;' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_frequentist_sub.sh'
	echo 'cd '$DIR >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_frequentist_sub.sh'
	echo 'source /afs/cern.ch/sw/lcg/contrib/gcc/4.3/x86_64-slc5/setup.sh' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_frequentist_sub.sh'
	echo 'source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.33.02/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_frequentist_sub.sh'
	echo '' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_frequentist_sub.sh'
	echo 'root -q -l -b runKristenLimit.C'"'"'("'$cut'","'$modelName'",'$m0','$m12',1,0)'"'" >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_frequentist_sub.sh'
	echo 'if [ -e "/tmp/ra2b/ws_expected_frequentist_'$cut'_'$modelName'_'$m0'_'$m12'_kristen.root" ]; then rm /tmp/ra2b/ws_measured_frequentist_'$cut'_'$modelName'_'$m0'_'$m12'_kristen.root; fi' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_frequentist_sub.sh'
	chmod a+x 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_frequentist_sub.sh'
	echo "Submitting job "$modelName'_'$cut'_'$m0'_'$m12"_measured_kristen_frequentist to batch system for processing."
	#bsub -q 1nd -N -o /dev/null -J $modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_frequentist' -Jd $modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_frequentist' $DIR'/submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_frequentist_sub.sh'


	echo '#! /bin/bash -f' > 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_hybrid_sub.sh'
	echo 'export BATCH_DEBUG=true' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_hybrid_sub.sh'
	echo 'if [ ! -d "/tmp/ra2b" ]; then mkdir /tmp/ra2b; fi;' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_hybrid_sub.sh'
	#echo 'if [ -d "/tmp/ra2b" ]; then echo "found directory!"; fi;' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_hybrid_sub.sh'
	echo 'cd '$DIR >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_hybrid_sub.sh'
	echo 'source /afs/cern.ch/sw/lcg/contrib/gcc/4.3/x86_64-slc5/setup.sh' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_hybrid_sub.sh'
	echo 'source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.33.02/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_hybrid_sub.sh'
	echo '' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_hybrid_sub.sh'
	echo 'root -q -l -b runKristenLimit.C'"'"'("'$cut'","'$modelName'",'$m0','$m12',0,1)'"'" >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_hybrid_sub.sh'
	echo 'if [ -e "/tmp/ra2b/ws_expected_hybrid_'$cut'_'$modelName'_'$m0'_'$m12'_kristen.root" ]; then rm /tmp/ra2b/ws_expected_hybrid_'$cut'_'$modelName'_'$m0'_'$m12'_kristen.root; fi' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_hybrid_sub.sh'
	chmod a+x 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_hybrid_sub.sh'
	echo "Submitting job "$modelName'_'$cut'_'$m0'_'$m12"_expected_kristen_hybrid to batch system for processing."
	bsub -q 1nd -N -o /dev/null -J $modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_hybrid' -Jd $modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_hybrid' $DIR'/submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_expected_kristen_hybrid_sub.sh'


	echo '#! /bin/bash -f' > 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_hybrid_sub.sh'
	echo 'export BATCH_DEBUG=true' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_hybrid_sub.sh'
	echo 'if [ ! -d "/tmp/ra2b" ]; then mkdir /tmp/ra2b; fi;' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_hybrid_sub.sh'
	#echo 'if [ -d "/tmp/ra2b" ]; then echo "found directory!"; fi;' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_hybrid_sub.sh'
	echo 'cd '$DIR >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_hybrid_sub.sh'
	echo 'source /afs/cern.ch/sw/lcg/contrib/gcc/4.3/x86_64-slc5/setup.sh' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_hybrid_sub.sh'
	echo 'source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.33.02/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_hybrid_sub.sh'
	echo '' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_hybrid_sub.sh'
	echo 'root -q -l -b runKristenLimit.C'"'"'("'$cut'","'$modelName'",'$m0','$m12',1,1)'"'" >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_hybrid_sub.sh'
	echo 'if [ -e "/tmp/ra2b/ws_expected_hybrid_'$cut'_'$modelName'_'$m0'_'$m12'_kristen.root" ]; then rm /tmp/ra2b/ws_measured_hybrid_'$cut'_'$modelName'_'$m0'_'$m12'_kristen.root; fi' >> 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_hybrid_sub.sh'
	chmod a+x 'submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_hybrid_sub.sh'
	echo "Submitting job "$modelName'_'$cut'_'$m0'_'$m12"_measured_kristen_hybrid to batch system for processing."
	bsub -q 1nd -N -o /dev/null -J $modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_hybrid' -Jd $modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_hybrid' $DIR'/submit_scripts/'$modelName'_'$cut'_'$m0'_'$m12'_measured_kristen_hybrid_sub.sh'


    done
done
exit