#!/bin/bash

for beta in 1.5 1.7 1.9 2.1 ; do

    echo "======================================"
    echo "=== Calcul pour beta=$beta ==="
    echo "======================================"

    sed -i "s/float mybeta=.*/float mybeta=$beta;/" NoGraceConfine.cpp

    g++ -O3 -fpermissive -o NoGraceConfine NoGraceConfine.cpp CPULoop.cpp \
    -I ../../Include \
    -I /fermiqcd/Libraries \
    -I include

    if [ $? -ne 0 ]; then
        echo "ERREUR de compilation pour beta=$beta — on passe au suivant"
        continue
    fi

    ./NoGraceConfine

    mkdir -p results_beta_bis$beta
    cp Result_bis_* results_beta_bis$beta/ 2>/dev/null
    cp log.txt results_beta_bis$beta/log_beta$beta.txt 2>/dev/null
    cp tmppot.dat results_beta_bis$beta/tmppot_beta$beta.dat 2>/dev/null
    cp tmpfitpot.dat results_beta_bis$beta/tmpfitpot_beta$beta.dat 2>/dev/null

    echo "=== beta=$beta terminé ==="

done

echo "======================================"
echo "=== Tous les calculs sont terminés ==="
echo "======================================"