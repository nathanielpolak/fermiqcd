#include "MyHeaders.h"
#include "Functions.h"
#include "cppheaders.h"

// 1. DÉFINITIONS CONSTANTES POUR LES HEADERS (Fixe les erreurs de scope)
#define tmax 10
#define zmax 10
#define npar 3
#define nparpot 3

// 2. VARIABLES GLOBALES (Protège contre malloc/fastbin corruption)
char replace_tmp[11][50], replace_tmp_x[11][50], replace_tmp_y[11][50], replace_tmp_dy[11][50];
char replace_tmpfit[101][50], replace_tmpfit_as[101][50];
float averaged_loop[11][11], jack_error[11][11];
float ftvalue[11], fzvalue[11];
static float store_loop[501 * 10 * 10]; 

// Variables de fitting requises par ThingsForFitting.h
float potential[zmax], poterror[zmax];
Measure ydy[tmax];
Measure potdpot[zmax];

void run_simulation(int argc, char** argv) {
    // Configuration Lattice
    int box[] = {12, 12, 6, 6};
    mdp_lattice lattice(4, box);
    gauge_field U(lattice, 2); 
    coefficients gauge;
    gauge["beta"] = 2.5;

    // Paramètres de boucle
    int niter = 100;
    int nheat = 200;
    int tvalue[] = {1,2,3,4,5,6,7,8,9,10};
    int zvalue[] = {1,2,3,4,5,6,7,8,9,10};
    
    for(int i=0; i<10; i++) { 
        ftvalue[i] = (float)tvalue[i]; 
        fzvalue[i] = (float)zvalue[i]; 
    }

    ofstream logfile("log.txt");

    // INCLUSIONS CRITIQUES (Maintenant que tmax/zmax sont définis)
    #include "ThingsForPlotting.h"
    #include "ThingsForFitting.h"

    set_hot(U);

    for(int k=1; k<=(nheat + niter); k++) {
        if(k <= nheat) {
            WilsonGaugeAction::heatbath(U, gauge, 1);
            continue; 
        }

        WilsonGaugeAction::heatbath(U, gauge, 5);
        int cur_iter = k - nheat;

        // Mesure des boucles
        for(int zl=0; zl<zmax; zl++) {
            for(int tl=0; tl<tmax; tl++) {
                store_loop[cur_iter*100 + zl*10 + tl] = my_average_loop(U, 0, tvalue[tl], 1, zvalue[zl]);
            }
        }

        // Analyse Statistique
        if(cur_iter % 10 == 0) {
            for(int zloop=0; zloop<zmax; zloop++) {
                ofstream tmp("tmp.dat"); // Défini ici pour PrepareFit.h
                int data_index = -1;     // Défini ici pour PrepareFit.h

                for(int tloop=0; tloop<tmax; tloop++) {
                    mdp_jackboot jack(cur_iter, 2);
                    for(int kp=1; kp<=cur_iter; kp++) {
                        jack(kp-1, 1) = store_loop[kp*100 + zloop*10 + tloop];
                        jack(kp-1, 0) = (tloop == 0) ? 1.0 : store_loop[kp*100 + zloop*10 + (tloop-1)];
                    }
                    
                    float delta_t = (tloop == 0) ? (float)tvalue[0] : (float)(tvalue[tloop]-tvalue[tloop-1]);
                    jack.handle = (void*)&delta_t;
                    jack.f = f1;

                    averaged_loop[zloop][tloop] = jack.mean();
                    jack_error[zloop][tloop] = jack.j_err();

                    #include "FlagBad.h"
                    #include "PrepareFit.h"
                }
                tmp.close();

                #include "GetPot.h"
                if(zloop == zmax-1) {
                    #include "NoGraceFitPlot_Pot.h"
                    #include "NoGracePlot_a.h"
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    mdp.open_wormholes(argc, argv);
    run_simulation(argc, argv);
    mdp.close_wormholes();
    return 0;
}