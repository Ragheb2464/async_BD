#!/bin/bash -l

flags=-O
#flags=-g

mpic++ -DNDEBUG -DIL_STD -c main_async.cpp ${flags}
mpic++ ${flags} -o main_async main_async.o -lilocplex -lconcert -lcplex -lm -lpthread

instcateg=S
nb_proc=5

for nodes in 16 30
do
    for origin in 14
    do
	insname="instances-${instcateg}/SND${origin}C${nodes}N.dat"

	for scenarios in 1000 #10 20
	do
	    scenname="instances-${instcateg}/scenario${origin}C${scenarios}S.dat"
	    echo qsub -o res_${nodes}_${origin}_${scenarios}.out -e res_${nodes}_${origin}_${scenarios}.err -pe smp $nb_proc ./script.sh $instcateg $nodes $origin $scenarios $insname $scenname $nb_proc
	    qsub -o res_${nodes}_${origin}_${scenarios}.out -e res_${nodes}_${origin}_${scenarios}.err -pe smp $nb_proc ./script.sh $instcateg $nodes $origin $scenarios $insname $scenname $nb_proc
	done
    done
    for origin in 40 80
    do
	insname="instances-${instcateg}/SND${origin}C${nodes}N.dat"

	for scenarios in 1000 #20 60 90
	do
	    scenname="instances-${instcateg}/scenario${origin}C${scenarios}S.dat"
	    echo qsub -o res_${nodes}_${origin}_${scenarios}.out -e res_${nodes}_${origin}_${scenarios}.err -pe smp $nb_proc ./script.sh $instcateg $nodes $origin $scenarios $insname $scenname
	    qsub -o res_${nodes}_${origin}_${scenarios}.out -e res_${nodes}_${origin}_${scenarios}.err -pe smp $nb_proc ./script.sh $instcateg $nodes $origin $scenarios $insname $scenname $nb_proc
	done
 
    done
done
