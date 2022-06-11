#!/bin/bash -l
echo $"InsName" $"Time_I" $"Obj_I" $"Time_II" $"LB_II" $"UB_II" $"LB" $"UB" $"TotalTime(s)" $"gap(%)" $"nbIter"
for inst in $(ls *.out)
do
    time_root=`grep "the time spend on the relaxation" $inst  |cut -d " " -f 8`
    obj_root=`grep "the objective of the relaxation" $inst  |cut -d " " -f 7`
    time_interm=`grep "the time spend on the intermediary phase" $inst  |cut -d " " -f 8`
    lb_interm=`grep "The lower bound after intermediary phase" $inst  |cut -d " " -f 7`
    ub_interm=`grep "The upper bound after intermediary phase" $inst  |cut -d " " -f 7`
    glob_lb=`grep "The global lower bound" $inst  |cut -d " " -f 5`
    glob_ub=`grep "The global upper bound" $inst  |cut -d " " -f 5`
    run_time=`grep "The whole run time is" $inst  |cut -d " " -f 6`
    gap=`grep "Optimality gap in % is" $inst  |cut -d " " -f 6`
    nb_iter=`grep "The whole run time is" $inst  |cut -d " " -f 12`
    echo $inst $time_root $obj_root $time_interm $lb_interm $ub_interm $glob_lb $glob_ub $run_time $gap $nb_iter
done
