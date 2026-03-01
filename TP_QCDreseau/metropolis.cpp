// Program3: Integration in D dimensions without bounds
////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// prob: unnormalized probability distribution in D/
////////////////////////////////////////////////////////////

// Fonction qui génère un nombre aléatoire dans [a,b]
float myrandom(float a, float b)
{
    float range = pow(2.0, 31);     // 2^31
    return a + (b - a) * random() / range;
}

float prob(int d,float *x){
    float arg=0.;
    for(int i=0;i<d;i++){arg=arg+x[i]*x[i];}
    return exp(-arg);
}
////////////////////////////////////////////////////////////
// func: function multiplying the proba//
////////////////////////////////////////////////////////////
float func(int d,float *x){
    float arg=0.;
    for(int i=0;i<d;i++){arg=arg+x[i]*x[i];}
    return arg;
}
/////////
//MAIN//
/////////
int main(){
    int j,dim=2; // set dimension
    float epsi=1.2; // set parameter for update
////////////////////////////////////////////////////////////
//2D MARKOV CHAIN. Store function//
////////////////////////////////////////////////////////////
    int i,nmarkov=100000;
    float x[dim],y[dim];
    float store_func[nmarkov+1];
    float accept=0.;
    for(j=0;j<dim;j++)
    {
        x[j]=0.;
    } // first point of markov chain

    for(i=1;i<=nmarkov;i++)
    {
        for(j=0;j<dim;j++)
        {
            y[j]=x[j]+myrandom(-epsi,epsi);
        }
        if( prob(dim,y)/prob(dim,x) > myrandom(0.,1.) )
        {
            for(j=0;j<dim;j++) x[j]=y[j];accept=accept+1.;
        }
        store_func[i]=func(dim,x);
        if(dim==1)        // write 1d chain for histogram
        {
            outfile2<<i<<" "<<x[0]<<endl;
        }
        if(dim==2)        // write 2d chain to watch walk
        {
            outfile3<<x[0]<<" "<<x[1]<<endl;
        }
    }
///////////////////////////////
//MEAN AND STANDARD ERROR//
///////////////////////////////
    float sum=0.,mean,sum2=0.,mean2,sigma,exact=dim*0.5;
    for(i=1;i<=nmarkov;i++)
    {
        sum=sum+store_func[i];mean=sum/i;
        sum2=sum2+store_func[i]*store_func[i];mean2=sum2/i;
        sigma=sqrt((mean2-mean*mean)/i);
        if(i%1000==0)
        {
            outfile1<<i<<" "          //write result
                    <<mean<<" "<<sigma<<endl;
        }
    }
    //
    cout<<" calcul en dimension: "<<dim<<endl;
    cout<<" accept rate: "<<accept/nmarkov<<endl;
    cout<<"mean: "<<mean<<" std error: "
        <<sigma<<" exact: "<<exact<<endl;

///////////////////
//AUTOCORRELATION//
///////////////////
ofstream outfile("autocorel_Nd.dat");
float sumshift,chi,chi0,normedchi,rough_time;
int icorel,ncorel=100;
for(icorel=0;icorel<=ncorel;icorel++)
{
    sumshift=0.;
    for(i=1;i<=nmarkov-icorel;i++)
    {
        sumshift=sumshift+store_func[i]*store_func[i+icorel];
    }
    chi=sumshift/(nmarkov-icorel)-mean*mean;
    if(icorel==0) chi0=chi;
    normedchi=chi/chi0;
    if(icorel==1)
    {
        rough_time=-1/log(normedchi);
    } // rough estimate of correlation time
    outfile<<icorel<<" "
           <<normedchi<<endl; // write correl. function
}
cout<<" rough estimate of exponential time: "
    <<rough_time<<endl;
//
outfile.close();
    //
    outfile1.close();outfile2.close();outfile3.close();
}
                                                        //End of Program3