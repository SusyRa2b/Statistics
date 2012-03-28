#! /bin/bash 

rm -f ./nominal_results/T1bbbb*
rm -f ./kristen_results/T1bbbb*
cp results/T1bbbb* nominal_results/
mv nominal_results/*kristen* kristen_results/

if [ -e T1bbbb.nominal.expected_frequentist.1BT.dat ] ; then rm T1bbbb.nominal.expected_frequentist.1BT.dat ; fi
cat nominal_results/T1bbbb_ge1bTight_*_expected_frequentist.dat > T1bbbb.nominal.expected_frequentist.1BT.dat
if [ -e T1bbbb.nominal.expected_frequentist.1BL.dat ] ; then rm T1bbbb.nominal.expected_frequentist.1BL.dat ; fi
cat nominal_results/T1bbbb_ge1bLoose_*_expected_frequentist.dat > T1bbbb.nominal.expected_frequentist.1BL.dat
if [ -e T1bbbb.nominal.expected_frequentist.2BT.dat ] ; then rm T1bbbb.nominal.expected_frequentist.2BT.dat ; fi
cat nominal_results/T1bbbb_ge2bTight_*_expected_frequentist.dat > T1bbbb.nominal.expected_frequentist.2BT.dat
if [ -e T1bbbb.nominal.expected_frequentist.2BL.dat ] ; then rm T1bbbb.nominal.expected_frequentist.2BL.dat ; fi
cat nominal_results/T1bbbb_ge2bLoose_*_expected_frequentist.dat > T1bbbb.nominal.expected_frequentist.2BL.dat
if [ -e T1bbbb.nominal.expected_frequentist.3B.dat ] ; then rm T1bbbb.nominal.expected_frequentist.3B.dat ; fi
cat nominal_results/T1bbbb_ge3bLoose_*_expected_frequentist.dat > T1bbbb.nominal.expected_frequentist.3B.dat

if [ -e T1bbbb.kristen.expected_frequentist.1BT.dat ] ; then rm T1bbbb.kristen.expected_frequentist.1BT.dat ; fi
cat kristen_results/T1bbbb_ge1bTight_*_expected_frequentist.dat > T1bbbb.kristen.expected_frequentist.1BT.dat
if [ -e T1bbbb.kristen.expected_frequentist.1BL.dat ] ; then rm T1bbbb.kristen.expected_frequentist.1BL.dat ; fi
cat kristen_results/T1bbbb_ge1bLoose_*_expected_frequentist.dat > T1bbbb.kristen.expected_frequentist.1BL.dat
if [ -e T1bbbb.kristen.expected_frequentist.2BT.dat ] ; then rm T1bbbb.kristen.expected_frequentist.2BT.dat ; fi
cat kristen_results/T1bbbb_ge2bTight_*_expected_frequentist.dat > T1bbbb.kristen.expected_frequentist.2BT.dat
if [ -e T1bbbb.kristen.expected_frequentist.2BL.dat ] ; then rm T1bbbb.kristen.expected_frequentist.2BL.dat ; fi
cat kristen_results/T1bbbb_ge2bLoose_*_expected_frequentist.dat > T1bbbb.kristen.expected_frequentist.2BL.dat
if [ -e T1bbbb.kristen.expected_frequentist.3B.dat ] ; then rm T1bbbb.kristen.expected_frequentist.3B.dat ; fi
cat kristen_results/T1bbbb_ge3bLoose_*_expected_frequentist.dat > T1bbbb.kristen.expected_frequentist.3B.dat

if [ -e T1bbbb.nominal.measured_frequentist.1BT.dat ] ; then rm T1bbbb.nominal.measured_frequentist.1BT.dat ; fi
cat nominal_results/T1bbbb_ge1bTight_*_measured_frequentist.dat > T1bbbb.nominal.measured_frequentist.1BT.dat
if [ -e T1bbbb.nominal.measured_frequentist.1BL.dat ] ; then rm T1bbbb.nominal.measured_frequentist.1BL.dat ; fi
cat nominal_results/T1bbbb_ge1bLoose_*_measured_frequentist.dat > T1bbbb.nominal.measured_frequentist.1BL.dat
if [ -e T1bbbb.nominal.measured_frequentist.2BT.dat ] ; then rm T1bbbb.nominal.measured_frequentist.2BT.dat ; fi
cat nominal_results/T1bbbb_ge2bTight_*_measured_frequentist.dat > T1bbbb.nominal.measured_frequentist.2BT.dat
if [ -e T1bbbb.nominal.measured_frequentist.2BL.dat ] ; then rm T1bbbb.nominal.measured_frequentist.2BL.dat ; fi
cat nominal_results/T1bbbb_ge2bLoose_*_measured_frequentist.dat > T1bbbb.nominal.measured_frequentist.2BL.dat
if [ -e T1bbbb.nominal.measured_frequentist.3B.dat ] ; then rm T1bbbb.nominal.measured_frequentist.3B.dat ; fi
cat nominal_results/T1bbbb_ge3bLoose_*_measured_frequentist.dat > T1bbbb.nominal.measured_frequentist.3B.dat

if [ -e T1bbbb.kristen.measured_frequentist.1BT.dat ] ; then rm T1bbbb.kristen.measured_frequentist.1BT.dat ; fi
cat kristen_results/T1bbbb_ge1bTight_*_measured_frequentist.dat > T1bbbb.kristen.measured_frequentist.1BT.dat
if [ -e T1bbbb.kristen.measured_frequentist.1BL.dat ] ; then rm T1bbbb.kristen.measured_frequentist.1BL.dat ; fi
cat kristen_results/T1bbbb_ge1bLoose_*_measured_frequentist.dat > T1bbbb.kristen.measured_frequentist.1BL.dat
if [ -e T1bbbb.kristen.measured_frequentist.2BT.dat ] ; then rm T1bbbb.kristen.measured_frequentist.2BT.dat ; fi
cat kristen_results/T1bbbb_ge2bTight_*_measured_frequentist.dat > T1bbbb.kristen.measured_frequentist.2BT.dat
if [ -e T1bbbb.kristen.measured_frequentist.2BL.dat ] ; then rm T1bbbb.kristen.measured_frequentist.2BL.dat ; fi
cat kristen_results/T1bbbb_ge2bLoose_*_measured_frequentist.dat > T1bbbb.kristen.measured_frequentist.2BL.dat
if [ -e T1bbbb.kristen.measured_frequentist.3B.dat ] ; then rm T1bbbb.kristen.measured_frequentist.3B.dat ; fi
cat kristen_results/T1bbbb_ge3bLoose_*_measured_frequentist.dat > T1bbbb.kristen.measured_frequentist.3B.dat






if [ -e T1bbbb.nominal.expected_hybrid.1BT.dat ] ; then rm T1bbbb.nominal.expected_hybrid.1BT.dat ; fi
cat nominal_results/T1bbbb_ge1bTight_*_expected_hybrid.dat > T1bbbb.nominal.expected_hybrid.1BT.dat
if [ -e T1bbbb.nominal.expected_hybrid.1BL.dat ] ; then rm T1bbbb.nominal.expected_hybrid.1BL.dat ; fi
cat nominal_results/T1bbbb_ge1bLoose_*_expected_hybrid.dat > T1bbbb.nominal.expected_hybrid.1BL.dat
if [ -e T1bbbb.nominal.expected_hybrid.2BT.dat ] ; then rm T1bbbb.nominal.expected_hybrid.2BT.dat ; fi
cat nominal_results/T1bbbb_ge2bTight_*_expected_hybrid.dat > T1bbbb.nominal.expected_hybrid.2BT.dat
if [ -e T1bbbb.nominal.expected_hybrid.2BL.dat ] ; then rm T1bbbb.nominal.expected_hybrid.2BL.dat ; fi
cat nominal_results/T1bbbb_ge2bLoose_*_expected_hybrid.dat > T1bbbb.nominal.expected_hybrid.2BL.dat
if [ -e T1bbbb.nominal.expected_hybrid.3B.dat ] ; then rm T1bbbb.nominal.expected_hybrid.3B.dat ; fi
cat nominal_results/T1bbbb_ge3bLoose_*_expected_hybrid.dat > T1bbbb.nominal.expected_hybrid.3B.dat

if [ -e T1bbbb.kristen.expected_hybrid.1BT.dat ] ; then rm T1bbbb.kristen.expected_hybrid.1BT.dat ; fi
cat kristen_results/T1bbbb_ge1bTight_*_expected_hybrid.dat > T1bbbb.kristen.expected_hybrid.1BT.dat
if [ -e T1bbbb.kristen.expected_hybrid.1BL.dat ] ; then rm T1bbbb.kristen.expected_hybrid.1BL.dat ; fi
cat kristen_results/T1bbbb_ge1bLoose_*_expected_hybrid.dat > T1bbbb.kristen.expected_hybrid.1BL.dat
if [ -e T1bbbb.kristen.expected_hybrid.2BT.dat ] ; then rm T1bbbb.kristen.expected_hybrid.2BT.dat ; fi
cat kristen_results/T1bbbb_ge2bTight_*_expected_hybrid.dat > T1bbbb.kristen.expected_hybrid.2BT.dat
if [ -e T1bbbb.kristen.expected_hybrid.2BL.dat ] ; then rm T1bbbb.kristen.expected_hybrid.2BL.dat ; fi
cat kristen_results/T1bbbb_ge2bLoose_*_expected_hybrid.dat > T1bbbb.kristen.expected_hybrid.2BL.dat
if [ -e T1bbbb.kristen.expected_hybrid.3B.dat ] ; then rm T1bbbb.kristen.expected_hybrid.3B.dat ; fi
cat kristen_results/T1bbbb_ge3bLoose_*_expected_hybrid.dat > T1bbbb.kristen.expected_hybrid.3B.dat

if [ -e T1bbbb.nominal.measured_hybrid.1BT.dat ] ; then rm T1bbbb.nominal.measured_hybrid.1BT.dat ; fi
cat nominal_results/T1bbbb_ge1bTight_*_measured_hybrid.dat > T1bbbb.nominal.measured_hybrid.1BT.dat
if [ -e T1bbbb.nominal.measured_hybrid.1BL.dat ] ; then rm T1bbbb.nominal.measured_hybrid.1BL.dat ; fi
cat nominal_results/T1bbbb_ge1bLoose_*_measured_hybrid.dat > T1bbbb.nominal.measured_hybrid.1BL.dat
if [ -e T1bbbb.nominal.measured_hybrid.2BT.dat ] ; then rm T1bbbb.nominal.measured_hybrid.2BT.dat ; fi
cat nominal_results/T1bbbb_ge2bTight_*_measured_hybrid.dat > T1bbbb.nominal.measured_hybrid.2BT.dat
if [ -e T1bbbb.nominal.measured_hybrid.2BL.dat ] ; then rm T1bbbb.nominal.measured_hybrid.2BL.dat ; fi
cat nominal_results/T1bbbb_ge2bLoose_*_measured_hybrid.dat > T1bbbb.nominal.measured_hybrid.2BL.dat
if [ -e T1bbbb.nominal.measured_hybrid.3B.dat ] ; then rm T1bbbb.nominal.measured_hybrid.3B.dat ; fi
cat nominal_results/T1bbbb_ge3bLoose_*_measured_hybrid.dat > T1bbbb.nominal.measured_hybrid.3B.dat

if [ -e T1bbbb.kristen.measured_hybrid.1BT.dat ] ; then rm T1bbbb.kristen.measured_hybrid.1BT.dat ; fi
cat kristen_results/T1bbbb_ge1bTight_*_measured_hybrid.dat > T1bbbb.kristen.measured_hybrid.1BT.dat
if [ -e T1bbbb.kristen.measured_hybrid.1BL.dat ] ; then rm T1bbbb.kristen.measured_hybrid.1BL.dat ; fi
cat kristen_results/T1bbbb_ge1bLoose_*_measured_hybrid.dat > T1bbbb.kristen.measured_hybrid.1BL.dat
if [ -e T1bbbb.kristen.measured_hybrid.2BT.dat ] ; then rm T1bbbb.kristen.measured_hybrid.2BT.dat ; fi
cat kristen_results/T1bbbb_ge2bTight_*_measured_hybrid.dat > T1bbbb.kristen.measured_hybrid.2BT.dat
if [ -e T1bbbb.kristen.measured_hybrid.2BL.dat ] ; then rm T1bbbb.kristen.measured_hybrid.2BL.dat ; fi
cat kristen_results/T1bbbb_ge2bLoose_*_measured_hybrid.dat > T1bbbb.kristen.measured_hybrid.2BL.dat
if [ -e T1bbbb.kristen.measured_hybrid.3B.dat ] ; then rm T1bbbb.kristen.measured_hybrid.3B.dat ; fi
cat kristen_results/T1bbbb_ge3bLoose_*_measured_hybrid.dat > T1bbbb.kristen.measured_hybrid.3B.dat


