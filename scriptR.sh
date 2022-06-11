#!/bin/bash -l

flags=-O
#flags=-g

mpic++ -DNDEBUG -DIL_STD -c main_async.cpp ${flags}
mpic++ ${flags} -o main_async main_async.o -lilocplex -lconcert -lcplex -lm -lpthread

instcateg=R
nb_proc=5

for Ins_R_name in 4 5 6 7 8 9 10 #11 ,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30
do

		if [ $Ins_R_name -eq 4 ]; then
			  nodes=10
			  n_arcs=60
			  origin=10
fi
		if [ $Ins_R_name -eq 5 ]; then
			  nodes=10
			  n_arcs=60
			  origin=25
fi
		if [ $Ins_R_name -eq 6 ]; then
			  nodes=10
			  n_arcs=60
			  origin=50
fi
		if [ $Ins_R_name -eq 7 ]; then
			  nodes=10
			  n_arcs=82
			  origin=10
fi
		if [ $Ins_R_name -eq 8 ]; then
			  nodes=10
			  n_arcs=83
			  origin=25
fi
		if [ $Ins_R_name -eq 9 ]; then
			  nodes=10
			  n_arcs=83
			  origin=50
fi
		if [ $Ins_R_name -eq 10 ]; then
			  nodes=20
			  n_arcs=120
			  origin=40
fi
		if [ $Ins_R_name -eq 11 ]; then
			  nodes=20
			  n_arcs=120
			  origin=100
fi
		if [ $Ins_R_name -eq 12 ]; then
			  nodes=20
			  n_arcs=120
			  origin=200
fi
		if [ $Ins_R_name -eq 13 ]; then
			  nodes=20
			  n_arcs=220
			  origin=40
fi
		if [ $Ins_R_name -eq 14 ]; then
			  nodes=20
			  n_arcs=220
			  origin=100
fi
		if [ $Ins_R_name -eq 15 ]; then
			  nodes=20
			  n_arcs=220
			  origin=200
fi
		if [ $Ins_R_name -eq 16 ]; then
			  nodes=20
			  n_arcs=314
			  origin=40
fi
		if [ $Ins_R_name -eq 17 ]; then
			  nodes=20
			  n_arcs=318
			  origin=100
fi
		if [ $Ins_R_name -eq 18 ]; then
			  nodes=20
			  n_arcs=315
			  origin=200
fi
		if [ $Ins_R_name -eq 19 ]; then
			  nodes=25
			  n_arcs=100
			  origin=10
fi
		if [ $Ins_R_name -eq 20 ]; then
			  nodes=25
			  n_arcs=400
			  origin=10
fi
		if [ $Ins_R_name -eq 21 ]; then
			  nodes=25
			  n_arcs=100
			  origin=30
fi
		if [ $Ins_R_name -eq 22 ]; then
			  nodes=100
			  n_arcs=400
			  origin=30
fi
		if [ $Ins_R_name -eq 23 ]; then
			  nodes=20
			  n_arcs=230
			  origin=40
fi
		if [ $Ins_R_name -eq 24 ]; then
			  nodes=20
			  n_arcs=294
			  origin=40
fi
		if [ $Ins_R_name -eq 25 ]; then
			  nodes=30
			  n_arcs=519
			  origin=100
fi
		if [ $Ins_R_name -eq 26 ]; then
			  nodes=30
			  n_arcs=680
			  origin=100
fi
		if [ $Ins_R_name -eq 27 ]; then
			  nodes=20
			  n_arcs=230
			  origin=200
fi
		if [ $Ins_R_name -eq 28 ]; then
			  nodes=20
			  n_arcs=294
			  origin=200
fi
		if [ $Ins_R_name -eq 29 ]; then
			  nodes=30
			  n_arcs=520
			  origin=300
fi
		if [ $Ins_R_name -eq 30 ]; then
			  nodes=30
			  n_arcs=685
			  origin=300
fi

    for FIXCOST_CAP_RATIOS in 1 3 9
    do
	insname="instances-${instcateg}/r0${Ins_R_name}.${FIXCOST_CAP_RATIOS}.dow"
	if [ $Ins_R_name -ge 10 ]; then
	    insname="instances-${instcateg}/r${Ins_R_name}.${FIXCOST_CAP_RATIOS}.dow"
	fi

	for CORRS in 0.2 #0 0.2 0.8 #[0,0.2,0.4,0.6,0.8]
	do
	     for SCENARIOS in 1000 #16 32 64
	     do
		scenname="instances-${instcateg}/scenarios/r0${Ins_R_name}-${CORRS}-${SCENARIOS}"
		if  [ $Ins_R_name -ge 10 ]; then
		     scenname="instances-${instcateg}/scenarios/r${Ins_R_name}-${CORRS}-${SCENARIOS}"
		fi
		echo qsub -o res_${Ins_R_name}_${FIXCOST_CAP_RATIOS}_${CORRS}_${SCENARIOS}.out -e res_${Ins_R_name}_${FIXCOST_CAP_RATIOS}_${CORRS}_${SCENARIOS}.err -pe smp $nb_proc ./script.sh $instcateg $nodes $origin $SCENARIOS $insname $scenname $nb_proc
		qsub -o res_${Ins_R_name}_${FIXCOST_CAP_RATIOS}_${CORRS}_${SCENARIOS}.out -e res_${Ins_R_name}_${FIXCOST_CAP_RATIOS}_${CORRS}_${SCENARIOS}.err -pe smp $nb_proc ./script.sh $instcateg $nodes $origin $SCENARIOS $insname $scenname $nb_proc

	     done
	done
    done    
done
