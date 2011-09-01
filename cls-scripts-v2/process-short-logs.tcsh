#!/bin/tcsh -f

  set selection = $1

  foreach d ( log-files-$1-t1bbbb-iter* )

     echo $d

     set tmp1 = `echo $d | sed s/log-files/output/`
     set outfile = ${tmp1}.txt
     echo $outfile

     cd $d

        zgrep final short*.gz |& tee summary.txt

        awk '{ printf(" %8d %8d %9.3f %9.3f %9.3f %8.2f\n", $4, $5, $6, $7, $8, $9 ) }' summary.txt > $outfile


     cd -

  end

