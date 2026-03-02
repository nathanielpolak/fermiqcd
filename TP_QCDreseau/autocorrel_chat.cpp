// Program: Metropolis + autocorrelations for dimensions 1..10
////////////////////////////////////////////////////////////
#include "MyHeaders.h"
#include <cstdlib>   // random(), srandom()
#include <cmath>     // exp(), log(), pow(), sqrt()
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

// Fonction qui génère un nombre aléatoire dans [a,b]
float myrandom(float a, float b)
{
    float range = pow(2.0, 31); // 2^31
    return a + (b - a) * random() / range;
}

// prob: unnormalized probability distribution in D
float prob(int d, const vector<float>& x)
{
    float arg = 0.0f;
    for (int i = 0; i < d; i++) arg += x[i] * x[i];
    return exp(-arg);
}

// func: function multiplying the proba (ici: x^2)
float func(int d, const vector<float>& x)
{
    float arg = 0.0f;
    for (int i = 0; i < d; i++) arg += x[i] * x[i];
    return arg;
}

int main()
{
    int nmarkov = 100000;
    int ncorel  = 100;   // autocorrelation up to tau=100
    float epsi  = 1.2f;  // set parameter for update

    // (Optionnel) initialisation du générateur si besoin :
    // srandom(1);

    for (int dim = 1; dim <= 10; dim++)
    {
        // Fichiers de sortie (un par dimension)
        string f_calcul = "calcul_Nd_D" + to_string(dim) + ".dat";
        string f_auto   = "autocorel_Nd_D" + to_string(dim) + ".dat";

        ofstream outfile1(f_calcul.c_str());
        ofstream outfileA(f_auto.c_str());

        // Chaîne de Markov
        vector<float> x(dim, 0.0f), y(dim, 0.0f);
        vector<float> store_func(nmarkov + 1, 0.0f);
        float accept = 0.0f;

        for (int i = 1; i <= nmarkov; i++)
        {
            for (int j = 0; j < dim; j++)
                y[j] = x[j] + myrandom(-epsi, epsi);

            if (prob(dim, y) / prob(dim, x) > myrandom(0.0f, 1.0f))
            {
                for (int j = 0; j < dim; j++) x[j] = y[j];
                accept = accept + 1.0f;
            }

            store_func[i] = func(dim, x);
        }

        // Moyenne et erreur standard (comme dans Program3)
        float sum = 0.0f, mean = 0.0f;
        float sum2 = 0.0f, mean2 = 0.0f;
        float sigma = 0.0f, exact = dim * 0.5f;

        for (int i = 1; i <= nmarkov; i++)
        {
            sum += store_func[i];
            mean = sum / i;

            sum2 += store_func[i] * store_func[i];
            mean2 = sum2 / i;

            sigma = sqrt((mean2 - mean * mean) / i);

            if (i % 1000 == 0)
                outfile1 << i << " " << mean << " " << sigma << endl;
        }

        cout << "Dimension: " << dim << endl;
        cout << "accept rate: " << accept / nmarkov << endl;
        cout << "mean: " << mean << " std error: " << sigma << " exact: " << exact << endl;

        // AUTOCORRELATION (comme dans le TP, Program4)
        float sumshift, chi, chi0, normedchi, rough_time;
        int icorel;

        for (icorel = 0; icorel <= ncorel; icorel++)
        {
            sumshift = 0.0f;

            for (int i = 1; i <= nmarkov - icorel; i++)
                sumshift = sumshift + store_func[i] * store_func[i + icorel];

            chi = sumshift / (nmarkov - icorel) - mean * mean;
            if (icorel == 0) chi0 = chi;

            normedchi = chi / chi0;

            if (icorel == 1)
                rough_time = -1.0f / log(normedchi); // rough estimate of correlation time

            outfileA << icorel << " " << normedchi << endl; // write correl. function
        }

        cout << "rough estimate of exponential time: " << rough_time << endl;
        cout << "----------------------------------------" << endl;

        outfile1.close();
        outfileA.close();
    }

    return 0;
}