// initialize model and produce the ws.root file

void make_ws() {

  gROOT->Reset();

  gROOT->LoadMacro("ra2bRoostatsClass3D_1.c+") ;
  ra2bRoostatsClass3D_1 *ra2b = new ra2bRoostatsClass3D_1() ;

  ra2b->initialize("Input.dat",
		   "Susy.dat",
		   300, 200, true, 1.0,
		   "dummy_DeffDbtag.dat");
  
  return;

}
