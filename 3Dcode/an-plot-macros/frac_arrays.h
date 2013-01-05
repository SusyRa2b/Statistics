

   int ncomps(4) ;
   char comp_out_name[4][100] = { "top and $W+$jets",
                                  "QCD",
                                  "\\znunu",
                                  "VV" } ;
   char comp_fracfile_name[4][100] = { "ttwj", "qcd", "znn", "vv" } ;


   float frac_2b_val[4][4] ; // first index is component, second is HT bin (column in table).
   float frac_2b_err[4][4] ; // first index is component, second is HT bin (column in table).

   float frac_3b_val[4][4][5] ; // indices: component, met bin (row), ht bin (column).
   float frac_3b_err[4][4][5] ; // indices: component, met bin (row), ht bin (column).

