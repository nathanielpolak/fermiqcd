// Program5: Program4 + Jackknife error versus block size
////////////////////////////////////////////////////////////
#include "MyHeaders.h"
#include <cstdlib>          // random()
#include <cmath>            // exp(), log(), pow(), sqrt()
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

// Fichiers de sortie
ofstream outfile1("calcul_Nd.dat");
ofstream outfile2("sample.dat");
ofstream outfile3("walk_2d.dat");

// Fonction qui génère un nombre aléatoire dans [a,b]
float myrandom(float a, float b)
{
    float range = pow(2.0, 31);     // 2^31
    return a + (b - a) * random() / range;
}

// prob: unnormalized probability distribution in D
float prob(int d, float *x)
{
    float arg = 0.0f;
    for(int i = 0; i < d; i++) arg += x[i] * x[i];
    return exp(-arg);
}

// func: function multiplying the proba (ici: x^2)
float func(int d, float *x)
{
    float arg = 0.0f;
    for(int i = 0; i < d; i++) arg += x[i] * x[i];
    return arg;
}

int main()
{
    int j, dim = 2;          // set dimension
    float epsi = 1.2f;       // set parameter for update

    // 2D MARKOV CHAIN. Store function
    int i, nmarkov = 100000, nstep = 5;   // nombre d'itérations à vide;
    float x[dim], y[dim];

    vector<float> store_func(nmarkov + 1, 0.0f);

    float accept = 0.0f;
    for(j = 0; j < dim; j++) x[j] = 0.0f; // first point of markov chain

    for(i = 1; i <= nmarkov; i++)
{
    // on fait nstep itérations à vide
    for(int istep = 0; istep < nstep; istep++)
    {
        for(j = 0; j < dim; j++)
        {
            y[j] = x[j] + myrandom(-epsi, epsi);
        }

        if(prob(dim, y) / prob(dim, x) > myrandom(0.0f, 1.0f))
        {
            for(j = 0; j < dim; j++) x[j] = y[j];
            accept = accept + 1.0f;
        }
    }

    // on mesure seulement après nstep updates
    store_func[i] = func(dim, x);

    if(dim == 1)
    {
        outfile2 << i << " " << x[0] << endl;
    }
    if(dim == 2)
    {
        outfile3 << x[0] << " " << x[1] << endl;
    }
}

    ///////////////////////////////
    // MEAN AND STANDARD ERROR   //
    ///////////////////////////////
    float sum = 0.0f, mean = 0.0f;
    float sum2 = 0.0f, mean2 = 0.0f;
    float sigma = 0.0f, exact = dim * 0.5f;

    for(i = 1; i <= nmarkov; i++)
    {
        sum += store_func[i];
        mean = sum / i;

        sum2 += store_func[i] * store_func[i];
        mean2 = sum2 / i;

        sigma = sqrt((mean2 - mean * mean) / i);

        if(i % 1000 == 0)
        {
            outfile1 << i << " " << mean << " " << sigma << endl;
        }
    }

    cout << " calcul en dimension: " << dim << endl;
    cout << " accept rate: " << accept / nmarkov << endl;
    cout << " mean: " << mean << " std error: " << sigma << " exact: " << exact << endl;

        ///////////////////
    //AUTOCORRELATION//
    ///////////////////
    ofstream outfileA("autocorel_Nd.dat");
    float sumshift, chi, chi0, normedchi, rough_time;
    int icorel, ncorel = 100;

    // Pour tau_int (TP: 1/2 + somme A(tau) tronquée avant fluctuations)
    float tau_int = 0.5f;
    int tau_cut = 0;

    for(icorel = 0; icorel <= ncorel; icorel++)
    {
        sumshift = 0.0f;

        for(i = 1; i <= nmarkov - icorel; i++)
        {
            sumshift = sumshift + store_func[i] * store_func[i + icorel];
        }

        chi = sumshift / (nmarkov - icorel) - mean * mean;

        if(icorel == 0) chi0 = chi;

        normedchi = chi / chi0;

        if(icorel == 1)
        {
            rough_time = -1.0f / log(normedchi);
        } // rough estimate of correlation time

        outfileA << icorel << " " << normedchi << endl;

        // --- intégration tronquée: on somme tant que A(tau) est "non négligeable" ---
        // On ne somme pas tau=0 (déjà pris en compte par le +1/2)
        if(icorel >= 1)
        {
            if(normedchi > 0.001f)
            {
                tau_int += normedchi;
                tau_cut = icorel;
            }
            else
            {
                break; // dès que ça passe <=0, on considère qu'on est dans les fluctuations
            }
        }
    }

    outfileA.close();

    float sigma_true = sigma * sqrt(2.0f * tau_int);

    cout << " rough estimate of exponential time: " << rough_time << endl;
    cout << " tau_cut = " << tau_cut << endl;
    cout << " tau_int (truncated) = " << tau_int << endl;
    cout << " sigma * sqrt(2 tau_int) = " << sigma_true << endl;

    //////////////////////////
    // JACKKNIFE ERROR (b)  //
    //////////////////////////
    // σ_JK(b) en fonction de la taille de bloc b
    ofstream outfileJK("jackknife_Nd.dat");

    int N = nmarkov;

    // Choix d'une plage de tailles de blocs (à adapter si besoin)
    // Ici on calcule pour b = 1..1000 (comme l'exemple du TP)
    for(int b = 1; b <= 1000; b++)
    {
        // Nombre de blocs (dernier bloc éventuellement plus petit)
        int Nb = (N + b - 1) / b; // ceil(N/b)

        // On calcule les moyennes de blocs fbar_i
        vector<float> fbar(Nb, 0.0f);
        vector<int>   bsize(Nb, 0);

        for(int ib = 0; ib < Nb; ib++)
        {
            int start = ib * b + 1;          // indices store_func: 1..N
            int end   = (ib + 1) * b;        // inclus
            if(end > N) end = N;

            int bi = end - start + 1;
            bsize[ib] = bi;

            float s = 0.0f;
            for(int k = start; k <= end; k++)
            {
                s += store_func[k];
            }
            fbar[ib] = s / bi;
        }

        // Jackknife: hat_f_i = (N * fbar_global - b_i * fbar_i) / (N - b_i)
        // et σ_JK^2 = (Nb-1)/Nb * Σ_i (hat_f_i - fbar_global)^2
        float sumsq = 0.0f;

        for(int ib = 0; ib < Nb; ib++)
        {
            int bi = bsize[ib];
            float hat_f = (N * mean - bi * fbar[ib]) / (N - bi);
            float diff  = hat_f - mean;
            sumsq += diff * diff;
        }

        float sigmaJK2 = ((float)(Nb - 1) / (float)Nb) * sumsq;
        float sigmaJK  = sqrt(sigmaJK2);

        outfileJK << b << " " << sigmaJK << endl;
    }

    outfileJK.close();

    outfile1.close();
    outfile2.close();
    outfile3.close();

    return 0;
}