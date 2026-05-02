// this is called by Confine.cpp   compile with svcc
//*******************************************************************************************************************************
/// includes and define
// ******************************************************************************************************************************
#include <stdio.h>
#include "parameters.h"
typedef struct {
    float re,im;
} mycomplex;
typedef struct {
   mycomplex m00,m01,m10,m11;} su2;
//
int field_adress(int t,int nt,int x,int nx,int i,int j,int nc)        // used for 2d calculation. periodicity is implemented
{int ct=t,cx=x;
   if(t>=nt)ct=t-nt;
   if(x>=nx)cx=x-nx;
   if(i>=nc||j>=nc)return(-1);
return(ct+nt*(cx+nx*(i+nc*j)));}
//
int full_adress(int t,int nt,int x,int nx,int y,int ny,int z,int nz,int i,int j,int nc)
{
   return(t+nt*(x+nx*(y+ny*(z+nz*(i+nc*j)))));               //no periodicity is necessary. On stores only the t,x  directions
}
int CCfull_adress(int t,int nt,int x,int nx,int y,int ny,int z,int nz,int i,int j,int nc)
{
   return(t+nt*(x+nx*(y+ny*(z+nz*(i+nc*j)))));               //no periodicity is necessary. On stores only the t,x  directions
}

int CCfield_adress(int t,int nt,int x,int nx,int i,int j,int nc)
{int ct=t,cx=x;
   if(t>=nt)ct=t-nt;
   if(x>=nx)cx=x-nx;
   if(i>=nc||j>=nc)return(-1);
return(ct+nt*(cx+nx*(i+nc*j)));}
mycomplex CCc_mult(mycomplex a,mycomplex b)
   {mycomplex c;
   c.re=a.re*b.re-a.im*b.im;
   c.im=a.re*b.im+a.im*b.re;
   return(c);}
su2 CCsu2_set(mycomplex a,mycomplex b,mycomplex c,mycomplex d)
   {su2 matrix;
   matrix.m00.re=a.re;matrix.m01.re=b.re;matrix.m10.re=c.re;matrix.m11.re=d.re;
   matrix.m00.im=a.im;matrix.m01.im=b.im;matrix.m10.im=c.im;matrix.m11.im=d.im;
   return(matrix);
   }
su2 CCsu2_equal(su2 left_hand_side)
   {su2 matrix;
   matrix.m00.re=left_hand_side.m00.re;
   matrix.m01.re=left_hand_side.m01.re;
   matrix.m10.re=left_hand_side.m10.re;
   matrix.m11.re=left_hand_side.m11.re;
   matrix.m00.im=left_hand_side.m00.im;
   matrix.m01.im=left_hand_side.m01.im;
   matrix.m10.im=left_hand_side.m10.im;
   matrix.m11.im=left_hand_side.m11.im;
   return(matrix);}
su2 CCsu2_adjoint(su2 left_hand_side)
   {su2 matrix;
   matrix.m00.re=left_hand_side.m00.re;
   matrix.m01.re=left_hand_side.m10.re;
   matrix.m10.re=left_hand_side.m01.re;
   matrix.m11.re=left_hand_side.m11.re;
   //
   matrix.m00.im=-left_hand_side.m00.im;
   matrix.m01.im=-left_hand_side.m10.im;
   matrix.m10.im=-left_hand_side.m01.im;
   matrix.m11.im=-left_hand_side.m11.im;
   return(matrix);}
// this is the product
su2 CCsu2_prod(su2 mat1,su2 mat2)
   {su2 matrix;
   matrix.m00.re=CCc_mult(mat1.m00,mat2.m00).re+CCc_mult(mat1.m01,mat2.m10).re;
   matrix.m10.re=CCc_mult(mat1.m10,mat2.m00).re+CCc_mult(mat1.m11,mat2.m10).re;
   matrix.m01.re=CCc_mult(mat1.m00,mat2.m01).re+CCc_mult(mat1.m01,mat2.m11).re;
   matrix.m11.re=CCc_mult(mat1.m10,mat2.m01).re+CCc_mult(mat1.m11,mat2.m11).re;
   matrix.m00.im=CCc_mult(mat1.m00,mat2.m00).im+CCc_mult(mat1.m01,mat2.m10).im;
   matrix.m10.im=CCc_mult(mat1.m10,mat2.m00).im+CCc_mult(mat1.m11,mat2.m10).im;
   matrix.m01.im=CCc_mult(mat1.m00,mat2.m01).im+CCc_mult(mat1.m01,mat2.m11).im;
   matrix.m11.im=CCc_mult(mat1.m10,mat2.m01).im+CCc_mult(mat1.m11,mat2.m11).im;
   return(matrix);
   }
// this compute the local loop at t0,x0  of length tloop x xloop  : tested
float betterCC_local_loop(mycomplex* Ut,mycomplex* Ux,int nt,int nx,int nc,int tloop,int xloop,int t0,int x0)
{int k;
su2 downloop,downfac,uploop,upfac,adjuploop,loop;
//tloop is the number of links along the t direction
if(tloop==0||xloop==0)return(1.);
//downloop
// initialize downloop  with first link
downloop=CCsu2_set(Ut[CCfield_adress(t0,nt,x0,nx,0,0,nc)],
       Ut[CCfield_adress(t0,nt,x0,nx,0,1,nc)],
       Ut[CCfield_adress(t0,nt,x0,nx,1,0,nc)],
       Ut[CCfield_adress(t0,nt,x0,nx,1,1,nc)]);
// make  products for  down horizontal link
for( k=1;k<tloop;k++)
   {
   downfac=CCsu2_set(Ut[CCfield_adress(t0+k,nt,x0,nx,0,0,nc)],
         Ut[CCfield_adress(t0+k,nt,x0,nx,0,1,nc)],
         Ut[CCfield_adress(t0+k,nt,x0,nx,1,0,nc)],
         Ut[CCfield_adress(t0+k,nt,x0,nx,1,1,nc)]);
   downloop=CCsu2_prod(downloop,downfac);
   }
// make  products for vertical right link
for( k=0;k<xloop;k++)
   {downfac=CCsu2_set(Ux[CCfield_adress(t0+tloop,nt,x0+k,nx,0,0,nc)],
      Ux[CCfield_adress(t0+tloop,nt,x0+k,nx,0,1,nc)],
      Ux[CCfield_adress(t0+tloop,nt,x0+k,nx,1,0,nc)],
      Ux[CCfield_adress(t0+tloop,nt,x0+k,nx,1,1,nc)]);
   downloop=CCsu2_prod(downloop,downfac);}
//uploop
//initialize uploop  with first vertical link
   uploop=CCsu2_set(Ux[CCfield_adress(t0,nt,x0,nx,0,0,nc)],
        Ux[CCfield_adress(t0,nt,x0,nx,0,1,nc)],
        Ux[CCfield_adress(t0,nt,x0,nx,1,0,nc)],
        Ux[CCfield_adress(t0,nt,x0,nx,1,1,nc)]);
// make  products for left vertical link
for( k=1;k<xloop;k++)
   {
   upfac=CCsu2_set(Ux[CCfield_adress(t0,nt,x0+k,nx,0,0,nc)],
       Ux[CCfield_adress(t0,nt,x0+k,nx,0,1,nc)],
       Ux[CCfield_adress(t0,nt,x0+k,nx,1,0,nc)],
       Ux[CCfield_adress(t0,nt,x0+k,nx,1,1,nc)]);
   uploop=CCsu2_prod(uploop,upfac);
   }
// make  products for up  horizontal link
for( k=0;k<tloop;k++)
{// put current U in X
   upfac=CCsu2_set(Ut[CCfield_adress(t0+k,nt,x0+xloop,nx,0,0,nc)],
         Ut[CCfield_adress(t0+k,nt,x0+xloop,nx,0,1,nc)],
         Ut[CCfield_adress(t0+k,nt,x0+xloop,nx,1,0,nc)],
         Ut[CCfield_adress(t0+k,nt,x0+xloop,nx,1,1,nc)]);
   uploop=CCsu2_prod(uploop,upfac);
}
//
// compute Tr(downlink*uplink^+)/nc
adjuploop=CCsu2_adjoint(uploop);
loop=CCsu2_prod(downloop,adjuploop);
float resR=0.,resI=0.;
resR=(loop.m00.re+loop.m11.re)/nc;
resI=(loop.m00.im+loop.m11.im)/nc;
return(resR);
}
// ********************************************************************************************************************************
/// Make average loops
// ********************************************************************************************************************************
int make_aver_loops_with_CPU(float *RfullCCfield_t,float *IfullCCfield_t,float *RfullCCfield_x,float *IfullCCfield_x,int *box,int nc,float* CCstore_aver_loops)
{
   if(box[0]!=BLOCK_SIZE || box[1]!=BLOCK_SIZE) {printf("t,x box size is not =BLOCK_SIZE !\n");return(1);}
   if(nc!=2){printf("nc is not=2!");return(1);}
   printf("lattice size: %d %d %d %d \n",box[0],box[1],box[2],box[3]);
///allocate memory to receive local  CCstore_loops
   float* CCstore_loops;
   CCstore_loops=(float*)malloc(box[0]*box[1]*box[2]*box[3]*LOOPSIZE*LOOPSIZE*sizeof(float));
/// load field in t,x plane for all lattice; this is used only for the cpu calculation now
   int size=box[0]*box[1]*nc*nc*sizeof(mycomplex);
   mycomplex *CCfield_t,*CCfield_x;
   CCfield_t=(mycomplex*)malloc(size);
   CCfield_x=(mycomplex*)malloc(size);
//
/// compute averaged loops starting from 1
   if(LOOPSIZE_MIN!=1){printf("LOOPSIZE_MIN must be 1 %d",LOOPSIZE_MIN);return(-1);}
   int adr,av_adr;
   int tloop,xloop,tloc,xloc,yloc,zloc;
   float volume=(box[0]*box[1]*box[2]*box[3]);
/// compute local loops
   for(yloc=0;yloc<box[2];yloc++)
   for(zloc=0;zloc<box[2];zloc++)
   {          int it,ix;
              for(it=0;it<box[0];it++)  //loop on site
              for(ix=0;ix<box[1];ix++)
              {int ic,jc;
              for( ic=0;ic<nc;ic++)
              for( jc=0;jc<nc;jc++)
              {int adress,fadress;
              adress=field_adress(it,box[0],ix,box[1],ic,jc,nc);
              fadress=CCfull_adress(it,box[0],ix,box[1],yloc,box[2],zloc,box[3],ic,jc,nc);
              CCfield_t[adress].re=RfullCCfield_t[fadress];
              CCfield_t[adress].im=IfullCCfield_t[fadress];
              CCfield_x[adress].re=RfullCCfield_x[fadress];
              CCfield_x[adress].im=IfullCCfield_x[fadress];
              }
              }
       for(tloop=1;tloop<=LOOPSIZE;tloop++)
       for(xloop=1;xloop<=LOOPSIZE;xloop++)
       for(tloc=0;tloc<box[0];tloc++)
       for(xloc=0;xloc<box[1];xloc++)
         {float betterloop;
         betterloop=betterCC_local_loop(CCfield_t,CCfield_x,box[0],box[1],nc,tloop,xloop,tloc,xloc);
         adr=tloc+box[0]*(xloc+box[1]*(yloc+box[2]*(zloc+box[3]*((tloop-1)+LOOPSIZE*(xloop-1)))));
         CCstore_loops[adr]=betterloop;
         }
   }
///compute average loops
   if(LOOPSIZE_MIN!=1){printf("LOOPSIZE_MIN must be 1 %d",LOOPSIZE_MIN);return(-1);}
   for(tloop=1;tloop<=LOOPSIZE;tloop++)
   for(xloop=1;xloop<=LOOPSIZE;xloop++)
      {float average=0;
      for(tloc=0;tloc<box[0];tloc++)                            //loop on site
      for(xloc=0;xloc<box[1];xloc++)
      for(yloc=0;yloc<box[2];yloc++)
      for(zloc=0;zloc<box[2];zloc++)
         {adr=tloc+box[0]*(xloc+box[1]*(yloc+box[2]*(zloc+box[3]*((tloop-1)+LOOPSIZE*(xloop-1)))));
         average=average+CCstore_loops[adr];
         }
      av_adr=(tloop-1)+LOOPSIZE*(xloop-1);
      CCstore_aver_loops[av_adr]=average/volume;
      }
//
free(CCfield_t);free(CCfield_x);free(CCstore_loops);
return(0);
}
//