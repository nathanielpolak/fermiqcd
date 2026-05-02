/// ////////////
/// INCLUDES //
/// /////////
#include "myheaders.h"                        // This includes all the big stuff for computing
                                             // The other includes are just for plotting or fitting
#include "Functions.h"
#include "cppheaders.h"
#if defined(__APPLE__)                        // to allow running on mac-os
#include <stdlib.h>
#else
#include <malloc.h>
#endif
/// /////////
/// MAIN ///
/// ///////
int on_exit();
int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);                // paralelism instruction
/// ///////////////////////
/// LATTICE DEFINITION ///
/// /////////////////////
  int ndim=4;
  int tsize=8,xsize=8,ysize=6,zsize=6;                       ///lattice size   tsize=xsize is compulsory
  if(tsize!=xsize){cout<<"tsize "<<tsize<<" is not equal to xsize "<<xsize<<"\n";abort();}
  int box[]={tsize,xsize,ysize,zsize};
  mdp_lattice lattice(4,box);                 //create lattice in 4 dimensions
  int nc=2;                                   //nc colors     must be 2 if systemchoice==0
  gauge_field U(lattice,nc);                  //create gauge field
  coefficients gauge;
  float mybeta=2.6;
  gauge["beta"]=mybeta;                       ///set coupling
  site x(U.lattice());
  int systemchoice=0;                        /// systemchoice=0   is  faster but implemented only  with nc=2
                                             ///systemchoice=2 use the slow c++ function but works for any nc
/// /////////////////////
/// ITERATION CONTROL //
/// ///////////////////
  int decorel=2,nheat=10;                    /// decorrelation and heating
  int maxiter=10000,niter=100;                /// iterations
  int k=0;                                   // iterations index
  int kupdate=10;                            // update graph every kupdate iterations
  long int tinitial,tfinal;                  // used for CPU info
/// /////////////////////////////////////////////////////////////////////////////////////////////////
///  LOOP DEFINITION and STORAGE
/// for compatibility reasons the loop size is still called (tvalue,zvalue) eventhough the loop plane is t,x now
/// /////////////////////////////////////////////////////////////////////////////////////////////////
//
  int maxloopsize=4;                          /// maximum loop size in t,x directions, not larger than 10
  int mu=0,nu=1;                              //loop plane
  int zmax=maxloopsize,tmax=maxloopsize;                        // number of zvalue and tvalue used
///  DONT CHANGE THIS
  int tvalue[]={1,2,3,4,5,6,7,8,9,10};              //end of loops in time direction This is kept for historical reasons. DONT CHANGE THIS
  int zvalue[]={1,2,3,4,5,6,7,8,9,10};              //end of loops in z direction
  if(zmax>10||tmax>10){cout<<"tmax,zmax: "<<tmax<<" "<<zmax<<"  are too large!\n";abort();}
//
/// ///////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////
/// NOTHING TO CHANGE  BEYOND THAT POINT //////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////
//
///  allocation for storing full lattice field (separately real and imag) and receive loops
//
   int fullsize=box[0]*box[1]*box[2]*box[3]*nc*nc;
   float * RfullCCfield_t,* IfullCCfield_t,* RfullCCfield_x,* IfullCCfield_x;
   RfullCCfield_t=(float*)malloc(fullsize*sizeof(float));
   IfullCCfield_t=(float*)malloc(fullsize*sizeof(float));
   RfullCCfield_x=(float*)malloc(fullsize*sizeof(float));
   IfullCCfield_x=(float*)malloc(fullsize*sizeof(float));
//
   float* CPU_average_loop;
   CPU_average_loop=(float*)malloc(maxloopsize*maxloopsize*sizeof(float));
//
  ofstream logfile("log.txt");
  if( tvalue[tmax-1]>=box[0]|| zvalue[zmax-1]>=box[1] )   //loop must be smaller than lattice!
  {system("xmessage -center  maximum loop larger than lattice! &");abort();}
  float ftvalue[tmax],fzvalue[zmax];
  for(int t=0;t<tmax;t++)ftvalue[t]=tvalue[t];        //casting
  for(int z=0;z<zmax;z++)fzvalue[z]=zvalue[z];        //casting
  float store_loop[maxiter+1][zmax][tmax];            //storing all the loops
  float averaged_loop[zmax][tmax],jack_error[zmax][tmax];   //storing the averaged loop and error
                                                            //The naming is incorrect because
	                                                    //now this  is log[<W(t,z)>/<W(t-1,z)>]
	                                                    //rather than <W(t,z)>
  #include "ThingsForPlotting.h"
  #include "ThingsForFitting.h"                         //obvious meaning
//  #include "GraceTitle_1.h"
/// ///////////////////////////
/// CALCULATION STARTS HERE //
/// /////////////////////////
    tinitial=time(0);
/// /////////////
///  HEATING ///
/// ///////////
    set_hot(U);
    for(k=1;k<=nheat;k++)
    {WilsonGaugeAction::heatbath(U,gauge,1);        //make nheat heatbath sweeps
//     #include "GraceTitle_2.h"
    }
//     #include "GraceTitle_3.h"
/// //////////////////////////
/// BEGIN  ITERATION LOOP ///
/// ////////////////////////
//
/// these char arrays are used to avoid reading a file when plotting
char replace_tmp[10][50];
char replace_tmp_x[10][50],replace_tmp_y[10][50],replace_tmp_dy[10][50];
char replace_tmpfit[100][50],replace_tmpfit_as[100][50];
//
float t_update,t_loop;
tinitial=time(0);
for(k=1; k<=niter  ; k++)                          //big iteration loop
{
 t_update=(float)clock();
 WilsonGaugeAction::heatbath(U,gauge,decorel);    // make decorel  sweeps of heatbath
 t_update=((float)clock()-t_update)/1000000;      //CPU time for heatbath
/// fill the arrays with the U just generated
for(int it=0;it<box[0];it++)  //loop on site
  for(int ix=0;ix<box[1];ix++)
     for(int iy=0;iy<box[2];iy++)
       for(int iz=0;iz<box[2];iz++)
    {x.set(it,ix,iy,iz);
       for(int ic=0;ic<nc;ic++)
       for(int jc=0;jc<nc;jc++)
       {int fadress;
       fadress=full_adress(it,box[0],ix,box[1],iy,box[2],iz,box[3],ic,jc,nc);
       RfullCCfield_t[fadress]=real(U(x,0,ic,jc));
       IfullCCfield_t[fadress]=imag(U(x,0,ic,jc));
       RfullCCfield_x[fadress]=real(U(x,1,ic,jc));
       IfullCCfield_x[fadress]=imag(U(x,1,ic,jc));
       }
    }
/// compute loops
/// if systemchoice==0 all loops sizes from 1 to maxloopsize are computed in a single shot and the stored in store_loop
//
   if(systemchoice==0)
   {t_loop=(float)clock();
    make_aver_loops_with_CPU(RfullCCfield_t,IfullCCfield_t,RfullCCfield_x,IfullCCfield_x,box,nc,maxloopsize,CPU_average_loop);
//
      for(int zloop=0;zloop<zmax;zloop++)              //computes and store the loops
      for(int tloop=0;tloop<tmax;tloop++)
         {
         int iloop=(tvalue[tloop]-1)+maxloopsize*(zvalue[zloop]-1);
         store_loop[k][zloop][tloop]=CPU_average_loop[iloop];
         }
     t_loop=((float)clock()-t_loop)/1000000;        //CPU time for loop calculation
   }
//
   if(systemchoice==2)                                  ///slow option, normally used only for checks
      {t_loop=(float)clock();
      for(int zloop=0;zloop<zmax;zloop++)              //computes and store the loops
      for(int tloop=0;tloop<tmax;tloop++)
	 {store_loop[k][zloop][tloop]=my_average_loop(U,mu,tvalue[tloop],nu,zvalue[zloop]);}
      t_loop=((float)clock()-t_loop)/1000000;
     }
//
if(k%kupdate==0)                                   // update plots every kupdate
    {
//    #include "GraceTitle_4.h"
    for(int zloop=0;zloop<zmax;zloop++)
	{ofstream tmp("tmp.dat");                  //stream for pollting v(t,z)
	int data_index=-1;                         // array index for good zvalues
	for(int tloop=0;tloop<tmax;tloop++)
	      {mdp_jackboot jack(k,2);                //create container jack
	      for(int kp=1;kp<=k;kp++)                //fill jack with data from 1 to current k
	           {jack(kp-1,1)=store_loop[kp][zloop][tloop];
	           jack(kp-1,0)=1.;
	           if(tloop>0){jack(kp-1,0)=store_loop[kp][zloop][tloop-1]; }
	           }
	      float delta_t=tvalue[tloop];
	      if(tloop>0) delta_t=tvalue[tloop]-tvalue[tloop-1];
	      jack.handle=(void*) &delta_t;                 //transmit  delta_t to jack
	      jack.f=f1;                                    //tell jack the fonction to use
	      averaged_loop[zloop][tloop]=jack.mean();
	      jack_error[zloop][tloop]=jack.j_err();
	      flagnumber[tloop]=0;
              #include "FlagBad.h"	                     //flag  bad points
              #include "PrepareFit.h"                        //prepare fit and plot of v(t,z)
	      }
	tmp.flush();tmp.close();
//
        #include "GetPot.h"              // fit v(t,z) to get the potential at z
//        #include "Plot_v_t_z.h"          // plot v(t,z) and the fit for several z
//
        if(zloop==zmax-1)     // Fit and plot potential when all z are calculated
            {
            #include "NoGraceFitPlot_Pot.h"      //fit the potential and plot it
            #include "NoGracePlot_a.h"           //plot the deduced value of a
            }
	}                     //end of z loop
//    #include "CPUtime.h"            //plot CPU time for heatbath and loop
    }                           //end of update block
// #include "AddIterations.h"      //ask for more iterations
}                               //end of iteration loop
/// ///////////
/// FINISH ///
/// /////////
logfile.close();
#include "WriteResults.h"                // write results on file Result_n, n incremented
//#include "Finish.h"                      // save file.agr with overwrite warning and close GRACE
/// /////////
/// EXIT ///
/// ///////
 mdp.close_wormholes();                         //parallelism instruction
}
