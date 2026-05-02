// Program3: Integration in D dimensions without bounds
////////////////////////////////////////////////////////////
#include "MyHeaders.h"
#include <cstdlib>          // pour random()
#include <cmath>            // pour pow(), sqrt(), exp(), log()
#include <iostream>
#include <fstream>
using namespace std;

ofstream outfile1("calcul_Nd_pas_corr.dat");        // erreur standard naive
ofstream outfile_corr("calcul_Nd_corr.dat"); // erreur corrigée autocorr
ofstream outfile2("sample.dat");
ofstream outfile3("walk_2d.dat");

// Fonction qui génère un nombre aléatoire dans [a,b]
float myrandom(float a, float b)
{
    float range = pow(2.0, 31);     // 2^31
    return a + (b - a) * random() / range;
}

float prob(int d, float *x){
    float arg=0.f;
    for(int i=0;i<d;i++){ arg += x[i]*x[i]; }
    return exp(-arg);
}

// f(x) = ||x||^2
float func(int d, float *x){
    float arg=0.f;
    for(int i=0;i<d;i++){ arg += x[i]*x[i]; }
    return arg;
}

int main(){
    int j, dim=2;         // set dimension
    float epsi=1.2f;      // set parameter for update

    int i, nmarkov=100000;
    float x[dim], y[dim];
    float store_func[nmarkov+1];

    float accept=0.f;
    for(j=0;j<dim;j++) x[j]=0.f; // first point

    // --- Markov chain
    for(i=1;i<=nmarkov;i++)
    {
        for(j=0;j<dim;j++) y[j] = x[j] + myrandom(-epsi, epsi);

        if( prob(dim,y)/prob(dim,x) > myrandom(0.f, 1.f) )
        {
            for(j=0;j<dim;j++) x[j]=y[j];
            accept += 1.f;
        }

        store_func[i]=func(dim,x);

        if(dim==1) outfile2 << i << " " << x[0] << "\n";
        if(dim==2) outfile3 << x[0] << " " << x[1] << "\n";
    }

    // --- Mean (global) needed for autocorrelation
    double mean_global = 0.0;
    for(i=1;i<=nmarkov;i++) mean_global += store_func[i];
    mean_global /= (double)nmarkov;

    // --- AUTOCORRELATION + tau_int (tronqué quand |A(tau)|<1e-3)
    ofstream outfile_ac("autocorel_Nd.dat");
    const int max_lag = 5000;      // assez large, mais pas énorme
    const double cut = 1e-3;       // comme dans ton CR
    double chi0 = 0.0;
    double tau_int = 0.5;          // 1/2 + somme_{tau>=1} A(tau)

    // chi(0)
    {
        double sum0 = 0.0;
        for(i=1;i<=nmarkov;i++) sum0 += (double)store_func[i]*(double)store_func[i];
        chi0 = sum0/(double)nmarkov - mean_global*mean_global;
        if(chi0 <= 0.0) chi0 = 1e-30; // sécurité numérique
        outfile_ac << 0 << " " << 1.0 << "\n";
    }

    int last_used_lag = 0;
    for(int lag=1; lag<=max_lag && lag < nmarkov; lag++)
    {
        double sumshift = 0.0;
        const int n = nmarkov - lag;
        for(i=1;i<=n;i++)
            sumshift += (double)store_func[i]*(double)store_func[i+lag];

        double chi = sumshift/(double)n - mean_global*mean_global;
        double A = chi/chi0;

        outfile_ac << lag << " " << A << "\n";

        // tronquage quand l'autocorr devient négligeable
        if (fabs(A) < cut) {
            last_used_lag = lag;
            break;
        }

        tau_int += A;
        last_used_lag = lag;
    }
    outfile_ac.close();

    const double corr_factor = sqrt(2.0 * tau_int);

    cout << " calcul en dimension: " << dim << endl;
    cout << " accept rate: " << accept/nmarkov << endl;
    cout << " mean_global: " << mean_global << endl;
    cout << " tau_int (cut |A|<" << cut << "): " << tau_int
         << "  (last lag used: " << last_used_lag << ")" << endl;
    cout << " correction factor sqrt(2*tau_int): " << corr_factor << endl;

    // --- MEAN + ERREURS (naive + corrigée) au fil de N, mêmes points que ton plot (tous les 1000)
    double sum=0.0, sum2=0.0;
    double mean=0.0, mean2=0.0, sigma=0.0;

    const double exact = dim * 0.5;

    for(i=1;i<=nmarkov;i++)
    {
        sum  += store_func[i];
        sum2 += (double)store_func[i]*(double)store_func[i];

        mean  = sum  / (double)i;
        mean2 = sum2 / (double)i;

        double var = mean2 - mean*mean;
        if(var < 0.0) var = 0.0;
        sigma = sqrt(var / (double)i);  // erreur standard naive

        if(i%1000==0)
        {
            outfile1      << i << " " << mean << " " << sigma << "\n";
            outfile_corr  << i << " " << mean << " " << (sigma*corr_factor) << "\n";
        }
    }

    cout << "mean(final): " << mean
         << " std error naive: " << sigma
         << " std error corrected: " << sigma*corr_factor
         << " exact: " << exact << endl;

    outfile1.close();
    outfile_corr.close();
    outfile2.close();
    outfile3.close();
}