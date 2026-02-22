// PROGRAM1: Use and test of random generator

#include "MyHeaders.h"      // Voir Préliminaires
#include <cstdlib>          // pour random()
#include <cmath>            // pour pow()
#include <iostream>
#include <fstream>

using namespace std;

ofstream outfile("random.dat");   // ouvre un fichier pour écrire

// Fonction qui génère un nombre aléatoire dans [a,b]
float myrandom(float a, float b)
{
    float range = pow(2.0, 31);     // 2^31
    return a + (b - a) * random() / range;
}

int main()
{
    float x, sum = 0.0, mean;
    int i, hitnb = 100000;   // nombre de tirages

    cout << "create " << hitnb 
         << " rand. numb. in [0,1] and compute mean" 
         << endl;

    for (i = 1; i <= hitnb; i++)
    {
        x = myrandom(0.0, 1.0);   // tirage aléatoire entre 0 et 1

        outfile << i << " " << x << endl;   // écrit dans le fichier

        sum = sum + x;     // somme cumulée
        mean = sum / i;    // moyenne

        if (i % (hitnb / 10) == 0)
        {
            cout << "i = " << i 
                 << " mean = " << mean 
                 << endl;
        }
    }

    outfile.close();   // ferme le fichier
    
     ifstream infile("random.dat");

    if (!infile) {
        cout << "Erreur ouverture fichier !" << endl;
        return 1;
    }

    double y;
    int index;
    double somme = 0.0;
    int N = 0;

    while (infile >> index >> y)
    {
        somme += exp(-pow(y,2.0)) * pow(y,2.0);   // f(y) à intégrer
        N++;
    }

    infile.close();

    double integral = somme / N;

    cout << "Nombre de points = " << N << endl;
    cout << "Intégrale ≈ " << integral << endl;

    return 0;
}


   
