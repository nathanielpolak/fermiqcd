// Program5_common_chain: one Markov chain, then JK for nstep in {1,5,10,20}
////////////////////////////////////////////////////////////
#include "MyHeaders.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

// random in [a,b]
float myrandom(float a, float b)
{
    float range = pow(2.0, 31);
    return a + (b - a) * random() / range;
}

// prob ~ exp(-x^2)
float prob(int d, float *x)
{
    float arg = 0.0f;
    for(int i = 0; i < d; i++) arg += x[i]*x[i];
    return exp(-arg);
}

// f(x) = x^2
float func(int d, float *x)
{
    float arg = 0.0f;
    for(int i = 0; i < d; i++) arg += x[i]*x[i];
    return arg;
}

// standard error on a series f[1..N]
void mean_and_sigma(const vector<float>& f, int N, float &mean, float &sigma)
{
    float sum = 0.0f, sum2 = 0.0f;
    mean = 0.0f; sigma = 0.0f;

    for(int i = 1; i <= N; i++)
    {
        sum += f[i];
        mean = sum / i;
        sum2 += f[i]*f[i];
        float mean2 = sum2 / i;
        sigma = sqrt((mean2 - mean*mean) / i);
    }
}

int main()
{
    int dim = 2;
    float epsi = 1.2f;

    // ---- paramètres ----
    int Nmeas = 100000;             // nb de mesures (comme le TP)
    int nstep_list[4] = {1,5,10,20};
    int nstep_max = 20;
    int ncorel = 100;
    int bmax = 1000;
    // --------------------

    // 1) Générer UNE SEULE chaîne "fine" (mesure à chaque update)
    int Nfine = Nmeas * nstep_max;  // nombre total d'updates
    vector<float> f_fine(Nfine + 1, 0.0f);

    float x[dim], y[dim];
    for(int j=0;j<dim;j++) x[j]=0.0f;

    float accept = 0.0f;

    for(int i = 1; i <= Nfine; i++)
    {
        for(int j = 0; j < dim; j++)
            y[j] = x[j] + myrandom(-epsi, epsi);

        if(prob(dim, y)/prob(dim, x) > myrandom(0.0f, 1.0f))
        {
            for(int j = 0; j < dim; j++) x[j] = y[j];
            accept += 1.0f;
        }

        // on "mesure" f à chaque update : chaîne fine
        f_fine[i] = func(dim, x);
    }

    // 2) sigma commun (calculé sur toute la chaîne fine)
    float mean_fine, sigma_common;
    mean_and_sigma(f_fine, Nfine, mean_fine, sigma_common);

    cout << "==============================" << endl;
    cout << "COMMON chain (fine) length Nfine = " << Nfine << endl;
    cout << "accept rate (per update): " << accept / Nfine << endl;
    cout << "mean_fine: " << mean_fine << endl;
    cout << "sigma_common (std error on full fine chain): " << sigma_common << endl;
    cout << "==============================" << endl;

    // 3) Pour chaque nstep, sous-échantillonnage + autocorr + jackknife
    for(int idx=0; idx<4; idx++)
    {
        int nstep = nstep_list[idx];

        // série mesurée: f_meas[k] = f_fine[k*nstep], k=1..Nmeas
        vector<float> f_meas(Nmeas + 1, 0.0f);
        for(int k=1; k<=Nmeas; k++)
            f_meas[k] = f_fine[k*nstep];

        // moyenne et sigma de CETTE série (si tu en as besoin),
        // mais pour ton graphe tu veux sigma_common identique pour tous.
        float mean_meas, sigma_meas;
        mean_and_sigma(f_meas, Nmeas, mean_meas, sigma_meas);

        // fichiers
        string f_auto = "autocorel_nstep" + to_string(nstep) + ".dat";
        string f_jk   = "jackknife_nstep" + to_string(nstep) + ".dat";
        ofstream outfileA(f_auto.c_str());
        ofstream outfileJK(f_jk.c_str());

        ///////////////////
        // AUTOCORRELATION (TP) + tau_int tronqué
        ///////////////////
        float chi0 = 0.0f;
        float rough_time = 0.0f;

        float tau_int = 0.5f;
        int tau_cut = 0;

        for(int tau = 0; tau <= ncorel; tau++)
        {
            float sumshift = 0.0f;
            for(int i=1; i<=Nmeas - tau; i++)
                sumshift += f_meas[i] * f_meas[i+tau];

            float chi = sumshift / (Nmeas - tau) - mean_meas*mean_meas;
            if(tau==0) chi0 = chi;

            float A = chi / chi0;

            if(tau==1) rough_time = -1.0f / log(A);

            outfileA << tau << " " << A << endl;

            if(tau >= 1)
            {
                if(A > 0.0f)
                {
                    tau_int += A;
                    tau_cut = tau;
                }
                else break; // tronquage avant fluctuations
            }
        }
        outfileA.close();

        float sigma_corr = sigma_common * sqrt(2.0f * tau_int);

        //////////////////////////
        // JACKKNIFE σ_JK(b) (TP)
        //////////////////////////
        int N = Nmeas;

        for(int b=1; b<=bmax; b++)
        {
            int Nb = (N + b - 1) / b; // ceil

            vector<float> fbar(Nb, 0.0f);
            vector<int> bsize(Nb, 0);

            for(int ib=0; ib<Nb; ib++)
            {
                int start = ib*b + 1;
                int end = (ib+1)*b;
                if(end > N) end = N;

                int bi = end - start + 1;
                bsize[ib] = bi;

                float s = 0.0f;
                for(int k=start; k<=end; k++) s += f_meas[k];
                fbar[ib] = s / bi;
            }

            float sumsq = 0.0f;
            for(int ib=0; ib<Nb; ib++)
            {
                int bi = bsize[ib];
                float hat_f = (N*mean_meas - bi*fbar[ib]) / (N - bi);
                float diff = hat_f - mean_meas;
                sumsq += diff*diff;
            }

            float sigmaJK2 = ((float)(Nb - 1) / (float)Nb) * sumsq;
            float sigmaJK  = sqrt(sigmaJK2);

            outfileJK << b << " " << sigmaJK << endl;
        }
        outfileJK.close();

        // résumé console
        cout << "nstep = " << nstep << endl;
        cout << "mean_meas = " << mean_meas << "  sigma_meas = " << sigma_meas << endl;
        cout << "tau_cut = " << tau_cut << "  tau_int = " << tau_int << endl;
        cout << "sigma_common * sqrt(2 tau_int) = " << sigma_corr << endl;
        cout << "--------------------------------" << endl;
    }

    return 0;
}