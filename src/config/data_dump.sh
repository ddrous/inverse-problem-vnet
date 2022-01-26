#!/bin/bash

# Script pour generer les donnees en variant le nom du fichier pour 'rho'
#   Ces variables ainsi sont ecrites dans 'src/config/io_names*.txt'
#   Les constantes pour ces simulation sont ecrites dans 'src/config/template_*.txt'
#   Les fichiers 'template_*' et 'ionames' forment le fichier de configuaration 'src/config/final_*.cfg'


# Nombre d'iterations totales
debut=3004
fin=3006
nbsimu=$(($fin-$debut))


# Simulation function
simulate() {
        echo "simulation numéro $h jusqu'à $(($fin-1)) en cours ..."

        OMP_NUM_THREADS=1 build/transfer src/config/final_u.cfg > /dev/null
        retVal=$?
        if [ $retVal -ne 0 ]; then
            error=1
        fi

        OMP_NUM_THREADS=1 build/transfer src/config/final_d.cfg > /dev/null
        retVal=$?
        if [ $retVal -ne 0 ]; then
            error=1
        fi

        OMP_NUM_THREADS=1 build/transfer src/config/final_l.cfg > /dev/null
        retVal=$?
        if [ $retVal -ne 0 ]; then
            error=1
        fi

        OMP_NUM_THREADS=1 build/transfer src/config/final_r.cfg > /dev/null
        retVal=$?
        if [ $retVal -ne 0 ]; then
            error=1
        fi
}


# Boucle des simulations
for (( h = $debut ; h < $fin; h++ )); do
    ## Buld file names
    input_name="rho numpy:data/inputs/"$h".npy\n"
    output_name="export_file data/outputs/"$h".sds\n"
    echo -e $input_name$output_name > src/config/io_names.txt
    cat src/config/io_names.txt src/config/template_u.txt > src/config/final_u.cfg
    cat src/config/io_names.txt src/config/template_d.txt > src/config/final_d.cfg
    cat src/config/io_names.txt src/config/template_l.txt > src/config/final_l.cfg
    cat src/config/io_names.txt src/config/template_r.txt > src/config/final_r.cfg
    
    ## Simulate
    error=0
    simulate    # can change the value of error to 1
    while [ $error -eq 1 ]
    do
        error=0
        simulate    
        # echo $error
    done

done