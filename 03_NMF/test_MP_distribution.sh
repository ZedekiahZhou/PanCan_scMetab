#! /bin/bash
paras=("p50_.7_.3" "p30_.7_.2" "p30_.7_.3" "p30_.7_.4" "p30_.8_.4" "p50_.7_.2" "p50_.7_.4" "p50_.8_.4")
MinScores=(0.5 0.8 1.0)

for para in ${paras[@]};
do
    for MinScore in ${MinScores[@]};
    do
        echo -e "______________________________\nBegein ${para} ${MinScore} at `date`\n------------------------------\n"
        Rscript src/2.5_NMF/MP_distribution.R ${para} ${MinScore}
        echo -e "------------------------------\nOver ${para} ${MinScore} at `date`\n______________________________\n"
    done
done