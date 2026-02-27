// Program3: Integration in D dimensions without bounds
/////////////////////////////////////////////////////////
#include "MyHeaders.h"
#include <cstdlib>          // pour random()
#include <cmath>            // pour pow()
#include <iostream>
#include <fstream>
using namespace std;
//
ofstream outfile1("calcul_Nd.dat");
ofstream outfile2("sample.dat");
ofstream outfile3("walk_2d.dat");
//FUNCTIONS:
// myrandom: random number in [a,b]//
float myrandom(float a, float b)
{
    float range = pow(2.0, 31);     // 2^31
    return a + (b - a) * random() / range;
}
/////////////////////////////////////////////////////////
// prob: unnormalized probability distribution in D/
/////////////////////////////////////////////////////////
float prob(int d,float *x){
float arg=0.;
for(int i=0;i<d;i++){arg=arg+x[i]*x[i];}
return exp(-arg);}
/////////////////////////////////////////////////////////
// func: function multiplying the proba//
/////////////////////////////////////////////////////////
float func(int d,float *x){
float arg=0.;
for(int i=0;i<d;i++){arg=arg+x[i]*x[i];}
return arg;}
////////
//MAIN/
////////
int main(){
int j,dim=20; // set dimension
float epsi=0.3; // set parameter for update

int Nmin=1;
int Nmax=100000;
int Nstep=1000;

float exact=dim*0.5;

// write a small header (comment style)
outfile1 << "# N mean sigma accept_rate exact" << endl;

for(int nmarkov=Nstep; nmarkov<=Nmax; nmarkov+=Nstep)
{
float x[dim],y[dim];
float accept=0.;

// first point of markov chain
for(j=0;j<dim;j++)
{
x[j]=0.;
}

// running sums (no need to store all values)
float sum=0., mean=0., sum2=0., mean2=0., sigma=0.;

for(int i=1;i<=nmarkov;i++)
{
for(j=0;j<dim;j++)
{
y[j]=x[j]+myrandom(-epsi,epsi);
}
if( prob(dim,y)/prob(dim,x) > myrandom(0.,1.) )
{
for(j=0;j<dim;j++) x[j]=y[j];accept=accept+1.;
}

// store sample output (optional): keep only last run files readable
if(nmarkov==Nmax)
{
if(dim==1) // write 1d chain for histogram
{
outfile2<<i<<" "<<x[0]<<endl;
}
if(dim==2) // write 2d chain to watch walk
{
outfile3<<x[0]<<" "<<x[1]<<endl;
}
}

// update running mean and standard error (same formula as before)
float fx = func(dim,x);
sum=sum+fx; mean=sum/i;
sum2=sum2+fx*fx; mean2=sum2/i;
sigma=sqrt((mean2-mean*mean)/i);
}

// write one line per chain length
outfile1 << nmarkov << " " << mean << " " << sigma << " " << accept/nmarkov << " " << exact << endl;
}

cout<<" calcul en dimension: "<<dim<<endl;
cout<<" exact: "<<exact<<endl;

outfile1.close();outfile2.close();outfile3.close();
}                                   //End of Program3