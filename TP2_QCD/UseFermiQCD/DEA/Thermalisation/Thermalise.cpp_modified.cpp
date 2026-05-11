//
///                        POUR CHANGER LA TAILLE DU RESEAU ALLER A LIGNE 104
///                        POUR CHANGER LA TAILLE DE LA GRANDE BOUCLE ALLER A LIGNE   127
/////////////
// INCLUDES //
/////////////
#include "MyHeaders.h"
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


//// MAIN
int main(int argc, char** argv) {
//////////////////////
//// PARALLEL STUFF //
/////////////////////
  mdp.open_wormholes(argc,argv);
///////////////////////////
//// LATTICE DEFINITION ///
///////////////////////////
  int box[]={16,12,6,6};                                      ///      define lattice
  mdp_lattice lattice(4,box);
  int nc=2;                                                   ///nc colors
  gauge_field U(lattice,nc);
  gauge_field Ucold(lattice,nc);
  coefficients gauge;
  float mybeta=2.5;             //set coupling
  gauge["beta"]=mybeta;
/////////////////////////
//// ITERATION CONTROL //
/////////////////////////
  int ncorel=1,nheat=0;       // decorelation and heating
  int maxiter=10000,niter=1000;
  int test, k;
  float rk;
  int autoscale=100,scalex=100;              //    autoscale every 10 iterations
  long int tinitial,tfinal;
//////////////////////////////////
//// LOOP DEFINITION and STORAGE//
//////////////////////////////////
  int mu=0,nu=1,t=1,z=1,bigt=4,bigz=4;                      ///loop plane and size
  bigt=5;bigz=5;                                           /// loop of size 5x5
  int   showbigloopscale=20;                                    /// multiply big loop by 20 for plotting
  float hotloop,hotsum,coldloop,coldsum;
  float bighotloop,bighotsum,bigcoldloop,bigcoldsum,bigloopscale;
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
  char graceloopsizehot[50],graceloopsizecold[50];
  bigloopscale=showbigloopscale;
//  sprintf(gracesubtitle,"subtitle \"SU(%d) lattice %dx%dx%dx%d beta=%0.1f\"",nc,box[0],box[1],box[2],box[3],mybeta);
//  sprintf(graceloopsizehot,"s2 legend \"%d*<%dx%d loop> hot\"",showbigloopscale,bigt-t,bigz-z);
// sprintf(graceloopsizecold,"s3 legend \"%d*<%dx%d loop> cold\"",showbigloopscale,bigt-t,bigz-z);
//
 // StartGrace();
 // GracePrintf("focus g0");
 // GracePrintf (gracesubtitle);
 // GracePrintf(graceloopsizehot);
//  GracePrintf(graceloopsizecold);
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
outfile1<<k<<" "<<hotloop<<" "<<hotstd<<endl;
outfile2<<k<<" "<<coldloop<<" "<<coldstd<<endl;
outfile3<<k<<" "<<bigloopscale*bighotloop<<" "<<bigloopscale*bighotstd<<endl;
outfile4<<k<<" "<<bigloopscale*bigcoldloop<<" "<<bigloopscale*bigcoldstd<<endl;
}
// PLOT IN LINE


   if(k%scalex==0 && k!=niter)
  // {sprintf(gracecommand,"world xmax %d",k+scalex);
 //  GracePrintf (gracecommand);
  // GracePrintf ("autoticks");
 //  }
   tfinal=time(0);
 //  sprintf(gracetitle,"title \"Thermalisation CPU time: %d secs\"",tfinal-tinitial);
 //  GracePrintf(gracetitle);

 //   GracePrintf ("g0.s0 point %0.5g, %0.5g", rk,hotloop);
  //  GracePrintf ("g0.s1 point %0.5g, %0.5g", rk,coldloop);
  //  GracePrintf ("g0.s2 point %0.5g, %0.5g", rk,bighotloop*bigloopscale);
  //  GracePrintf ("g0.s3 point %0.5g, %0.5g", rk,bigcoldloop*bigloopscale);
   // if(k%autoscale==0) GracePrintf ("autoscale yaxes");
  //  GracePrintf ("redraw");
//  add iterations interactively
//
	if(k==niter and niter<maxiter)
           {
            system("xmessage -print -center -buttons stop,+10?,+100?,+500?,+1000?  want to add iterations?>message ");
            ifstream tmpmess("message");
            string answer;
            getline(tmpmess,answer);tmpmess.close();
            cout<<answer<<endl;
            int increase;
            if(answer=="stop") increase=0;
            if(answer=="+10?") increase=10;
            if(answer=="+100?") increase=100;
            if(answer=="+500?") increase=500;
            if(answer=="+1000?") increase=1000;
            niter=niter+increase;
            if(niter>maxiter) niter=maxiter;

           }
  }
///////////////////////////
// END OF ITERATION LOOP.//
//////////////////////////////////
// CLOSE FILES AND PLOT ERRORS  //
//////////////////////////////////
outfile1.close();outfile2.close();outfile3.close();outfile4.close();
//sprintf(gracecommand,"world xmax %d",niter);
//GracePrintf ("with g1");
//GracePrintf (gracecommand);
//GracePrintf ("autoticks");
//GracePrintf ("focus g1");
//GracePrintf ("read xydy \"hot.dat\" ");
//GracePrintf ("read xydy \"cold.dat\" ");
//GracePrintf ("read xydy \"bighot.dat\" ");
//GracePrintf ("read xydy \"bigcold.dat\" ");
//GracePrintf ("redraw");
//
//system("xmessage -center  save Grace file Thermalise.agr and then exit?");
//
//EndGrace();
//
 mdp.close_wormholes();
//
return 0;
}
