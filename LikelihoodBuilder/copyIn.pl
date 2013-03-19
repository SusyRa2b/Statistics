#!/user/bin/perl
use strict;
use warnings;
use File::Path;

my $dir = $ARGV[1];

open( my $flist, "<", "$dir/$dir.dat") or die "Can't open input file!";

while(<$flist>) {
  
  if(/(\S+)(toy)(\d+)/) {
    
    my $string = "$1$2";
    my $toy = $3;
    my $copyto = "$string$toy";
    print "$copyto\n";
    
    my $num = $toy/1;
    print "$num\n";
 
    #copy in binFilesFile and binFilesFileMR
    system("cp copyIn/binFilesFile.dat $copyto/.");
    system("cp copyIn/binFilesFileMR.dat $copyto/.");

    #copy in shape syst files
    system("cp copyIn/singletop-xsec-shapesyst.txt $copyto/.");
    system("cp copyIn/wjets-xsec-shapesyst.txt $copyto/.");

    #copy in sig directory
    system("cp -r copyIn/sig1/ $copyto/.");

    #copy in sig directory
    system("cp -r copyIn/dbdy/ $copyto/.");
    
    #copy in toy counts
    system("cp -r $ARGV[0]/toyMR_$num/ $copyto/countsMR/");
    
  }

}
