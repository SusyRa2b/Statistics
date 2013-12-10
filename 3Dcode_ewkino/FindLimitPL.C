// run a very simple PL scan (using ComputeTestStat.C)
// and find the 1-sided 95% CL upper limit

#include <iostream>
#include <fstream>

void FindLimitPL(TString wsfile, TString outFile) {

  gROOT->LoadMacro("ComputeTestStat.C");

  double startVal(10), curVal(-9);
  double curTestStat(-9.), prevTestStat(-9);
  bool FoundIt(false);
  int iter(0);

  curVal = startVal ;

  while (!FoundIt) {

    prevTestStat = curTestStat ;
    prevVal = curVal ;

    curTestStat = ComputeTestStat(wsfile,curVal);

    if ( curTestStat < 1.0 ) curVal = curVal*1.5 ;
    else if ( curTestStat < 1.8 ) curVal = curVal*1.15 ;
    else if ( curTestStat < 2.3 ) curVal = curVal*1.07 ;
    else if ( curTestStat < 2.5 ) curVal = curVal*1.03 ;
    else curVal = curVal*1.01 ;

    iter++ ;
    if ( curTestStat > 2.70 || iter > 100 ) FoundIt = true ;

  }


  // now get the limit simply by linearly interpolating the two last points

  cout << "\n\nDegugging:" << endl ;
  cout << "prevVal = " << prevVal << endl ;
  cout << "curVal  = " << curVal << endl ;
  cout << "prevTestStat = " << prevTestStat << endl ;
  cout << "curTestStat  = " << curTestStat << endl ;


  double limit = curVal + (( 2.70 - curTestStat )*( prevVal - curVal )/( prevTestStat - curTestStat ));

  cout << "\n\n 95% CL upper limit: = " << limit << "\n\n" << endl ;


  // save the limit into a simple text file:

  TString fname = "PL_results/" ;
  fname += outFile ;

  ofstream outStream ;
  outStream.open(fname,ios::app) ;
  
  outStream << "PL_limit = " << limit << endl ;
  
  outStream.close() ;

  return ;

}
