//
/////////////
// INCLUDES //
/////////////
#include "myheaders.h"
#if defined(__APPLE__)                        // to allow running on mac-os
#include <stdlib.h>
#else
#include <malloc.h>
#endif
/////////////////
//  OPEN STREAMS //
/////////////////
ofstream outfile1("hot.dat");
ofstream outfile2("cold.dat");
ofstream outfile3("bighot.dat");
ofstream outfile4("bigcold.dat");
//
void StartGrace(){
    if (GraceOpen(2048) == -1) {
        fprintf (stderr, "Can't run Grace. \n");
        exit (EXIT_FAILURE);
    }
    /* Send some initialization commands to Grace */
//
    GracePrintf ("page size 700 500");
    GracePrintf ("arrange(1,2,0.2,0.1,0.1)");
//
    GracePrintf( "with g0");                    // g0
    GracePrintf("autoscale onread none");
    GracePrintf ("world xmin 0");
    GracePrintf ("world xmax 100");
    GracePrintf ("world ymin 0");
    GracePrintf ("world ymax 0.8");
    GracePrintf ("autoticks");
    //set symbols   for g0
    GracePrintf ("s0 symbol 1");
    GracePrintf ("s0 symbol size 0.4");
    GracePrintf ("s0 linestyle 0");
    GracePrintf ("s1 symbol 2");
    GracePrintf ("s1 symbol size 0.4");
    GracePrintf ("s1 linestyle 0");
    GracePrintf ("s2 symbol 3");
    GracePrintf ("s2 symbol size 0.4");
    GracePrintf ("s2 linestyle 0");
    GracePrintf ("s3 symbol 4");
    GracePrintf ("s3 symbol size 0.4");
    GracePrintf ("s3 linestyle 0");
//set legend
    GracePrintf ("legend on");
    GracePrintf ("legend loctype view");
    GracePrintf ("legend 0.85,0.95");
    GracePrintf ("s0 legend \"<plaquette> hot start\" ");
    GracePrintf ("s1 legend \"<plaquette> cold start\" ");
    GracePrintf ("s2 legend \"10*<4x4 loop> hot start\" ");
    GracePrintf ("s3 legend \"10*<4x4 loop> cold start\" ");
//
   GracePrintf( "with g1");           //  g1
    GracePrintf("autoscale onread none");
    GracePrintf ("world xmin 0");
    GracePrintf ("world xmax 100");
    GracePrintf ("world ymin 0");
    GracePrintf ("world ymax 0.8");
    GracePrintf ("autoticks");
    //set symbols  same for g1
    GracePrintf ("s0 symbol 1");
    GracePrintf ("s0 symbol size 0.4");
    GracePrintf ("s0 linestyle 0");
    GracePrintf ("s1 symbol 2");
    GracePrintf ("s1 symbol size 0.4");
    GracePrintf ("s1 linestyle 0");
    GracePrintf ("s2 symbol 3");
    GracePrintf ("s2 symbol size 0.4");
    GracePrintf ("s2 linestyle 0");
    GracePrintf ("s3 symbol 4");
    GracePrintf ("s3 symbol size 0.4");
    GracePrintf ("s3 linestyle 0");
}
void EndGrace(){
//   Close GRACE or pipe
    if (GraceIsOpen()) {
        /* Tell Grace to save the data */
        GracePrintf ("saveall \"Thermalise.agr\"");
        /* Flush the output buffer and close Grace */
        GraceClose();
        // Close pipe but keep grace open: ne pas oublier de quitter grace sinon  les process
        //s'accumulent
        //GraceClosePipe();
        /* We are done */
        exit (EXIT_SUCCESS);
    } else {
        exit (EXIT_FAILURE);
    }
}
//// MAIN
int main(int argc, char** argv) {
//////////////////////
//// PARALLEL STUFF //
/////////////////////
  mdp.open_wormholes(argc,argv);
///////////////////////////
//// LATTICE DEFINITION ///
///////////////////////////
  int timesize=8,spacesize=8;
///   int box[]={timesize,spacesize,spacesize,spacesize};
  int box[]={12,12,6,6};
  mdp_lattice lattice(4,box);
  int nc=2;                      //nc colors
  gauge_field U(lattice,nc);
  gauge_field Ucold(lattice,nc);
  coefficients gauge;
  float mybeta=2.5;             //set coupling
  gauge["beta"]=mybeta;
/////////////////////////
//// ITERATION CONTROL //
/////////////////////////
  int ncorel=1,nheat=0;       // decorelation and heating
  int maxiter=10000,niter=300;
  int test, k;
  float rk;
  int autoscale=100,scalex=100;              //    autoscale every 10 iterations
  long int tinitial,tfinal;
//////////////////////////////////
//// LOOP DEFINITION and STORAGE//
//////////////////////////////////
  int mu=0,nu=1,t=1,z=1,bigt=4,bigz=4;  //loop plane and size
  bigt=6;bigz=6;
  float hotloop,hotsum,coldloop,coldsum;
  float bighotloop,bighotsum,bigcoldloop,bigcoldsum,bigloopscale=10;
  float bighotstd,bigcoldstd,hotstd,coldstd;
  float store_hot_loop[maxiter+1],store_cold_loop[maxiter+1];
  float store_big_hot_loop[maxiter+1],store_big_cold_loop[maxiter+1];
/////////////////////////
// THINGS FOR PLOTTING //
/////////////////////////
  char filename[200],makedir[50],graphtitle[200],dirname[30];
  char gracesubtitle[100];
  char gracetitle[100];
  char gracecommand[50];
  sprintf(gracesubtitle,"subtitle \"SU(%d) lattice %dx%dx%dx%d beta=%0.1f\"",nc,box[0],box[1],box[2],box[3],mybeta);
//
  StartGrace();
  GracePrintf("focus g0");
  GracePrintf (gracesubtitle);
  //
////////////////////////////
//CALCULATION STARTS HERE //
////////////////////////////
   tinitial=time(0);
   set_hot(U);
   set_cold(Ucold);
/// Eventually make a heating step to see the difference
   for(k=1; k<=nheat  ; k++) {WilsonGaugeAction::heatbath(U,gauge,1);}
   for(k=1; k<=nheat  ; k++) {WilsonGaugeAction::heatbath(Ucold,gauge,1);}
   hotsum=0.;coldsum=0.;   bighotsum=0.; bigcoldsum=0.;       //initialize the loops
///////////////////////////
// BEGIN  ITERATION LOOP //
///////////////////////////
   for(k=1; k<=niter  ; k++) {
   rk=k;
// MAKE 1 SWEEP OF HEATBATH OVER ALL THE LATTICE
   WilsonGaugeAction::heatbath(U,gauge,ncorel);
   WilsonGaugeAction::heatbath(Ucold,gauge,ncorel);
// COMPUTES THE LOOPS AT TIME k
   	store_hot_loop[k]=my_average_loop(U,mu,t,nu,z);
   	hotsum=hotsum+store_hot_loop[k];
   	store_cold_loop[k]=my_average_loop(Ucold,mu,t,nu,z);
   	coldsum=coldsum+store_cold_loop[k];
	hotloop=hotsum/k;
	coldloop=coldsum/k;
	store_big_hot_loop[k]=my_average_loop(U,mu,bigt,nu,bigz);
	bighotsum=bighotsum+store_big_hot_loop[k];
	store_big_cold_loop[k]=my_average_loop(Ucold,mu,bigt,nu,bigz);
   	bigcoldsum=bigcoldsum+store_big_cold_loop[k];
	bighotloop=bighotsum/k;
	bigcoldloop=bigcoldsum/k;
// COMPUTE STD ERROR ON FLY
         hotstd=0.;coldstd=0.;   bighotstd=0.; bigcoldstd=0.;
	 for(int l=1;l<=k;l++){
	 hotstd=hotstd+pow(store_hot_loop[l]-hotloop,2);
	 coldstd=coldstd+pow(store_cold_loop[l]-coldloop,2);
	 bighotstd=bighotstd+pow(store_big_hot_loop[l]-bighotloop,2);
	 bigcoldstd=bigcoldstd+pow(store_big_cold_loop[l]-bigcoldloop,2);  }
	 hotstd=sqrt(hotstd)/k;
	 coldstd=sqrt(coldstd)/k;
	 bighotstd=sqrt(bighotstd)/k;
	 bigcoldstd=sqrt(bigcoldstd)/k;
// WRITE RESULTS EVERY 10 STEPS
if(k%10==0)
{
outfile1<<k<<" "<<hotloop<<" "<<2*hotstd<<endl;
outfile2<<k<<" "<<coldloop<<" "<<2*coldstd<<endl;
outfile3<<k<<" "<<10*bighotloop<<" "<<10*bighotstd<<endl;
outfile4<<k<<" "<<10*bigcoldloop<<" "<<10*bigcoldstd<<endl;
}
// PLOT IN LINE


   if(k%scalex==0 && k!=niter)
   {sprintf(gracecommand,"world xmax %d",k+scalex);
   GracePrintf (gracecommand);
   GracePrintf ("autoticks");
   }
   tfinal=time(0);
   sprintf(gracetitle,"title \"Thermalisation CPU time: %d secs\"",tfinal-tinitial);
   GracePrintf(gracetitle);

    GracePrintf ("g0.s0 point %0.5g, %0.5g", rk,hotloop);
    GracePrintf ("g0.s1 point %0.5g, %0.5g", rk,coldloop);
    GracePrintf ("g0.s2 point %0.5g, %0.5g", rk,bighotloop*bigloopscale);
    GracePrintf ("g0.s3 point %0.5g, %0.5g", rk,bigcoldloop*bigloopscale);
    if(k%autoscale==0) GracePrintf ("autoscale yaxes");
    GracePrintf ("redraw");

  }
///////////////////////////
// END OF ITERATION LOOP.//
//////////////////////////////////
// CLOSE FILES AND PLOT ERRORS  //
//////////////////////////////////
outfile1.close();outfile2.close();outfile3.close();outfile4.close();
sprintf(gracecommand,"world xmax %d",niter);
GracePrintf ("with g1");
GracePrintf (gracecommand);
GracePrintf ("autoticks");
GracePrintf ("focus g1");
GracePrintf ("read xydy \"hot.dat\" ");
GracePrintf ("read xydy \"cold.dat\" ");
GracePrintf ("read xydy \"bighot.dat\" ");
GracePrintf ("read xydy \"bigcold.dat\" ");
GracePrintf ("redraw");
//
system("xmessage -center  save Grace file Thermalise.agr and then exit?");
//
EndGrace();
//
 mdp.close_wormholes();
//
return 0;
}
