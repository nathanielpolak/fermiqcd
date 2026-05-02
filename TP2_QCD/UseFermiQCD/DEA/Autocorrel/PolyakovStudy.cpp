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
ofstream outfile1("Jacknife.dat");
ofstream outfile2("Autocorrel.dat");
ofstream outfile3("Autotime.dat");
 ofstream outfile4("Zero.dat");
//
void StartGrace(){
    if (GraceOpen(2048) == -1) {
        fprintf (stderr, "Can't run Grace. \n");
        exit (EXIT_FAILURE);
    }
    /* Send some initialization commands to Grace */
//
    GracePrintf ("page size 700 500");
    GracePrintf ("arrange(2,2,0.1,0.2,0.25)");
//g0
    GracePrintf( "with g0");
    GracePrintf("autoscale onread none");
    GracePrintf ("world xmin 0");
    GracePrintf ("world xmax 100");
    GracePrintf ("world ymin 0");
    GracePrintf ("world ymax 0.8");
    GracePrintf ("autoticks");
    GracePrintf ("xaxis label \"iteration time\"");
    //set symbols   for g0
    GracePrintf ("s0 symbol 1");
    GracePrintf ("s0 symbol size 0.4");
    GracePrintf ("s0 linestyle 0");
//set legend for g0
    GracePrintf ("legend on");
    GracePrintf ("legend loctype view");
    GracePrintf ("legend 0.2,0.7");
    GracePrintf ("s0 legend \"<Pol. loop>\" ");
//g1
    GracePrintf( "with g1");
    GracePrintf("autoscale onread none");
    GracePrintf ("world xmin 0");
    GracePrintf ("world xmax 100");
    GracePrintf ("world ymin 0");
    GracePrintf ("world ymax 0.8");
    GracePrintf ("autoticks");
    GracePrintf ("xaxis label \"block size\"");
    //set symbols   for g1
    GracePrintf ("s0 symbol 1");
    GracePrintf ("s0 symbol size 0.4");
    GracePrintf ("s0 linestyle 1");
//set legend for g1
    GracePrintf ("legend on");
    GracePrintf ("legend loctype view");
    GracePrintf ("legend 0.9,0.85");
    GracePrintf ("s0 legend \"jack knife error\" ");
//g2
    GracePrintf( "with g2");
    GracePrintf("autoscale onread none");
    GracePrintf ("world xmin 0");
    GracePrintf ("world xmax 100");
    GracePrintf ("world ymin 0");
    GracePrintf ("world ymax 0.8");
    GracePrintf ("autoticks");
    GracePrintf ("xaxis label \"time shift \"");
//set symbols  for g2
    GracePrintf ("s0 symbol 1");
    GracePrintf ("s0 symbol size 0.4");
    GracePrintf ("s0 linestyle 1");
    GracePrintf ("s1 symbol 0");
    GracePrintf ("s1 linestyle 1");
    GracePrintf ("s1 line color 1");
//set legend for g2
    GracePrintf ("legend on");
    GracePrintf ("legend loctype view");
    GracePrintf ("legend 0.2,0.45");
    GracePrintf ("s0 legend \"Autocorrelation function\" ");
//g3
    GracePrintf( "with g3");
    GracePrintf("autoscale onread none");
    GracePrintf ("world xmin 0");
    GracePrintf ("world xmax 100");
    GracePrintf ("world ymin 0");
    GracePrintf ("world ymax 0.8");
    GracePrintf ("autoticks");
    GracePrintf ("xaxis label \"max. time shift\"");
//set symbols  for g3
    GracePrintf ("s0 symbol 1");
    GracePrintf ("s0 symbol size 0.4");
    GracePrintf ("s0 linestyle 1");
//set legend for g3
    GracePrintf ("legend on");
    GracePrintf ("legend loctype view");
    GracePrintf ("legend 0.9,0.45");
    GracePrintf ("s0 legend \"Autocorrelation time\" ");
}
void EndGrace(){
//   Close GRACE or pipe
    if (GraceIsOpen()) {
        /* Tell Grace to save the data */
        GracePrintf ("saveall \"PolyakovStudy.agr\"");
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
  int timesize=4,spacesize=10;
  int box[]={timesize,spacesize,spacesize,spacesize};
  mdp_lattice lattice(4,box);
  int nc=2;                      //nc colors
  gauge_field U(lattice,nc);
  gauge_field Ucold(lattice,nc);
  coefficients gauge;
  float mybeta=2.31;             //set coupling
  gauge["beta"]=mybeta;
/////////////////////////
//// ITERATION CONTROL //
/////////////////////////
  int decorel=5,nheat=200;       // decorrelation and heating
  int maxiter=1000,niter=500;
  int test, k;
  float rk;
  int autoscale=100,scalex=100;              //    autoscale every 100 iterations
  long int tinitial,tfinal;
//////////////////////////////////
//// LOOP DEFINITION and STORAGE//
//////////////////////////////////
  int mu=0,nu=1,t=4,z=4,bigt=4,bigz=4;  //loop plane and size
  float  loop,sum,store_loop[maxiter+1];
/////////////////////////
// THINGS FOR PLOTTING //
/////////////////////////
  char gracesubtitle[100];
  char gracetitle[100];
  char gracecommand[50];
  sprintf(gracesubtitle,"subtitle \" SU(%d) lattice %dx%dx%dx%d beta=%0.2f, decorel=%d\" ",nc,box[0],box[1],box[2],box[3],mybeta,decorel);
//
  StartGrace();
  GracePrintf("focus g0");
  GracePrintf (gracesubtitle);
  //
////////////////////////////
//CALCULATION STARTS HERE //
////////////////////////////
   tinitial=time(0);
//   set_hot(U);
   set_cold(U);
   sum=0.;    //initialize the loops
   int x1=4,x2=2,x3=6;    //loop location
///////////////////////////
// BEGIN  ITERATION LOOP //
///////////////////////////
   for(k=1; k<=niter  ; k++) {
   rk=k;
// MAKE 1 SWEEP OF HEATBATH OVER ALL THE LATTICE
         WilsonGaugeAction::heatbath(U,gauge,decorel);
// COMPUTES THE LOOPS AT MARKOV TIME k.................works only for su2 (real)
//
   	store_loop[k]=real(my_polyakov_loop(U,x1,x2,x3));
//        store_loop[k]=real(my_polyakov_loop(U));
   	sum=sum+store_loop[k];
	loop=sum/k;
//
// PLOT IN LINE
//
   if(k%scalex==0 && k!=niter)
   {sprintf(gracecommand,"world xmax %d",k+scalex);
   GracePrintf (gracecommand);
   GracePrintf ("autoticks");
   }
   tfinal=time(0);
   sprintf(gracetitle,"title \"    Autocorrelation.   CPU time: %d secs\"",tfinal-tinitial);
   GracePrintf(gracetitle);
   GracePrintf ("g0.s0 point %0.5g, %0.5g", rk,loop);
   GracePrintf ("redraw");
   if(k==5 || k==30 || k%autoscale==0)
   {GracePrintf ("autoscale yaxes");GracePrintf ("autoticks");}
  }
///////////////////////////
// END OF ITERATION LOOP //.
//////////////////////////////////
//JACKKNIFE //
//////////////
//Jack knife parameters
    float hat_loop=0.,jack_err=0.,jack_err1,newloop;
    float bloc_av[maxiter+1],hat_bloc_av[maxiter+1];
    int bloc_size=1,bloc_size_number=50;
    int nconf=niter,nbloc;

    for(bloc_size=1;bloc_size<=nconf/10;bloc_size++){
    if(nconf%bloc_size==0) nbloc=(nconf/bloc_size);
    if(nconf%bloc_size!=0) nbloc=(nconf/bloc_size)+1;
//
//    BLOC AVERAGES
    for (int ibloc=1;ibloc<nbloc;ibloc++)   //average of normal size blocs
    {bloc_av[ibloc]=0.;
    for (int jbloc=(ibloc-1)*bloc_size+1;jbloc<=ibloc*bloc_size;jbloc++)
    {bloc_av[ibloc]=bloc_av[ibloc]+store_loop[jbloc];}
    bloc_av[ibloc]=bloc_av[ibloc]/bloc_size;}
    bloc_av[nbloc]=0.;                      //average last bloc
    for(int jbloc=(nbloc-1)*bloc_size+1;jbloc<=nconf;jbloc++)
    {bloc_av[nbloc]=bloc_av[nbloc]+store_loop[jbloc];}
    bloc_av[nbloc]=bloc_av[nbloc]/(nconf-(nbloc-1)*bloc_size);
    newloop=0.;                       // recomputes loop for locality and check
    for (int ibloc=1;ibloc<=nbloc;ibloc++)
    {newloop=newloop+bloc_av[ibloc];}
    newloop=newloop/nbloc;
//    cout<<" loop: "<<loop<<" newloop: "<<newloop<<endl;
//    BLOC COMPLEMENTARY AVERAGES and check loop  ( THIS USES loop.  not very good)
    hat_loop=0.;
    for (int ibloc=1;ibloc<=nbloc;ibloc++)
    {hat_bloc_av[ibloc]=(nconf*loop-bloc_size*bloc_av[ibloc])/(nconf-bloc_size);
    hat_loop=hat_loop+hat_bloc_av[ibloc];}
    hat_loop=hat_loop/nbloc;
//
    jack_err=0.;
    for (int ibloc=1;ibloc<=nbloc;ibloc++)
    {jack_err=jack_err+pow((hat_bloc_av[ibloc]-hat_loop),2);}
    jack_err=sqrt(jack_err*(nbloc-1)/nbloc);
    if (bloc_size==1)jack_err1=jack_err;
//
   outfile1<<bloc_size<<" "<<jack_err<<endl;

    }
/////////////////////////////////////////////
// CLOSE FILES AND PLOT JACK KNIFE ERRORS  //
/////////////////////////////////////////////
outfile1.close();
sprintf(gracecommand,"world xmax %d",bloc_size);
GracePrintf ("with g1");
GracePrintf (gracecommand);
sprintf(gracecommand,"world ymax %0.5f",2*jack_err1);
GracePrintf (gracecommand);
GracePrintf ("autoticks");
GracePrintf ("focus g1");
GracePrintf ("read xy \"Jacknife.dat\" ");
GracePrintf ("redraw");
/////////////////////////
// AUTOCORRELATION TIME//
/////////////////////////
int nmarkov=niter,ncorel=nmarkov/20;            //ncorel<<nmarkov
float sumshift,chi,chi0,normedchi,int_time,min_normedchi=0.;
if(ncorel*10>nmarkov) cout<<" Warming: nmarkov too small vs ncorrel"<<endl;
int_time=0.5 ;
for(int icorel=0;icorel<=ncorel;icorel++)
{sumshift=0.;
for(k=1;k<=nmarkov-icorel;k++)
{
sumshift=sumshift+store_loop[k]*store_loop[k+icorel];}
chi=sumshift/(nmarkov-icorel)-newloop*newloop;
if(icorel==0) chi0=chi;
normedchi=chi/chi0;
if(normedchi<min_normedchi)min_normedchi=normedchi;
if(icorel>=1) int_time=int_time+ normedchi;
outfile2<<icorel<<" "<<normedchi<<endl;       // write correlation function
outfile3<<icorel<<" "<<int_time<<endl;        // write correlation time
outfile4<<icorel<<" "<<0<<endl;
}
outfile2.close();outfile3.close();outfile4.close();
/////////
// PLOT //
/////////
GracePrintf ("with g2");
sprintf(gracecommand,"world ymax 1");
GracePrintf (gracecommand);
sprintf(gracecommand,"world ymin %1.1f",1.5*min_normedchi);
GracePrintf (gracecommand);
sprintf(gracecommand,"world xmax %d",ncorel);
GracePrintf (gracecommand);
GracePrintf ("autoticks");
GracePrintf ("focus g2");
GracePrintf ("read xy \"Autocorrel.dat\" ");
GracePrintf ("read xy \"Zero.dat\" ");
GracePrintf ("redraw");
//
GracePrintf ("with g3");
sprintf(gracecommand,"world xmax %d",ncorel);
GracePrintf (gracecommand);
GracePrintf ("autoticks");
GracePrintf ("focus g3");
GracePrintf ("read xy \"Autotime.dat\" ");
GracePrintf ("autoscale yaxes");
GracePrintf ("redraw");
//////////
// EXIT //
//////////
system("xmessage -center  save Grace file PolyakovStudy.agr and then exit?");
//
EndGrace();
//
 mdp.close_wormholes();
//
return 0;
}
