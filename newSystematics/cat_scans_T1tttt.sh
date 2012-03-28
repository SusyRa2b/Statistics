#! /bin/bash 

rm -f ./nominal_results/T1tttt*

cp results/T1tttt* nominal_results/

if [ -e T1tttt.nominal.expected.1BT.dat ] ; then rm T1tttt.nominal.expected.1BT.dat ; fi
cat nominal_results/T1tttt_ge1bTight_*_expected.dat > T1tttt.nominal.expected.1BT.dat
if [ -e T1tttt.nominal.expected.1BL.dat ] ; then rm T1tttt.nominal.expected.1BL.dat ; fi
cat nominal_results/T1tttt_ge1bLoose_*_expected.dat > T1tttt.nominal.expected.1BL.dat
if [ -e T1tttt.nominal.expected.2BT.dat ] ; then rm T1tttt.nominal.expected.2BT.dat ; fi
cat nominal_results/T1tttt_ge2bTight_*_expected.dat > T1tttt.nominal.expected.2BT.dat
if [ -e T1tttt.nominal.expected.2BL.dat ] ; then rm T1tttt.nominal.expected.2BL.dat ; fi
cat nominal_results/T1tttt_ge2bLoose_*_expected.dat > T1tttt.nominal.expected.2BL.dat
if [ -e T1tttt.nominal.expected.3B.dat ] ; then rm T1tttt.nominal.expected.3B.dat ; fi
cat nominal_results/T1tttt_ge3bLoose_*_expected.dat > T1tttt.nominal.expected.3B.dat

if [ -e T1tttt.nominal.measured.1BT.dat ] ; then rm T1tttt.nominal.measured.1BT.dat ; fi
cat nominal_results/T1tttt_ge1bTight_*_measured.dat > T1tttt.nominal.measured.1BT.dat
if [ -e T1tttt.nominal.measured.1BL.dat ] ; then rm T1tttt.nominal.measured.1BL.dat ; fi
cat nominal_results/T1tttt_ge1bLoose_*_measured.dat > T1tttt.nominal.measured.1BL.dat
if [ -e T1tttt.nominal.measured.2BT.dat ] ; then rm T1tttt.nominal.measured.2BT.dat ; fi
cat nominal_results/T1tttt_ge2bTight_*_measured.dat > T1tttt.nominal.measured.2BT.dat
if [ -e T1tttt.nominal.measured.2BL.dat ] ; then rm T1tttt.nominal.measured.2BL.dat ; fi
cat nominal_results/T1tttt_ge2bLoose_*_measured.dat > T1tttt.nominal.measured.2BL.dat
if [ -e T1tttt.nominal.measured.3B.dat ] ; then rm T1tttt.nominal.measured.3B.dat ; fi
cat nominal_results/T1tttt_ge3bLoose_*_measured.dat > T1tttt.nominal.measured.3B.dat
