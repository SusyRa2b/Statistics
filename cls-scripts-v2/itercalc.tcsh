#!/bin/tcsh -f

   echo processing $1

   if ( ! -f $1 ) then
      echo Not a file.
      exit -1
   endif


   set outfile1 = `echo $1 | sed s/scan-points-50GevBins/iter1-calc/`
   echo output file for iter1 : $outfile1
   awk '{ printf(" %8d %8d %9.2f\n", $1, $2, $5) }' $1 > $outfile1

   set outfile2 = `echo $1 | sed s/scan-points-50GevBins/iter2-calc/`
   echo output file for iter2 : $outfile2
   awk '{ printf(" %8d %8d %9.2f\n", $1, $2, (1.26*$5)) }' $1 > $outfile2

   set outfile3 = `echo $1 | sed s/scan-points-50GevBins/iter3-calc/`
   echo output file for iter3 : $outfile3
   awk '{ printf(" %8d %8d %9.2f\n", $1, $2, (1.59*$5)) }' $1 > $outfile3

   set outfile4 = `echo $1 | sed s/scan-points-50GevBins/iter4-calc/`
   echo output file for iter4 : $outfile4
   awk '{ printf(" %8d %8d %9.2f\n", $1, $2, (2.00*$5)) }' $1 > $outfile4

   set outfile5 = `echo $1 | sed s/scan-points-50GevBins/iter5-calc/`
   echo output file for iter5 : $outfile5
   awk '{ printf(" %8d %8d %9.2f\n", $1, $2, (2.52*$5)) }' $1 > $outfile5

   set outfile6 = `echo $1 | sed s/scan-points-50GevBins/iter6-calc/`
   echo output file for iter6 : $outfile6
   awk '{ printf(" %8d %8d %9.2f\n", $1, $2, (3.18*$5)) }' $1 > $outfile6

   set outfile7 = `echo $1 | sed s/scan-points-50GevBins/iter7-calc/`
   echo output file for iter7 : $outfile7
   awk '{ printf(" %8d %8d %9.2f\n", $1, $2, (4.00*$5)) }' $1 > $outfile7

   set outfile8 = `echo $1 | sed s/scan-points-50GevBins/iter8-calc/`
   echo output file for iter8 : $outfile8
   awk '{ printf(" %8d %8d %9.2f\n", $1, $2, (10.00*$5)) }' $1 > $outfile8




