export LOOPITERSKIP_PFUPDATE='1'
export LOOPCHANGE_PFUPDATE='1'

export LOOPITERSKIP_IEINSIDE_OUTER='1'
export LOOPCHANGE_IEINSIDE_OUTER='2'

export LOOPITERSKIP_IEINSIDE_INNER='1'
export LOOPCHANGE_IEINSIDE_INNER='1'

export LOOPITERSKIP_IEEDGE_OUTER='1'
export LOOPCHANGE_IEEDGE_OUTER='1'

export LOOPITERSKIP_IEEDGE_INNER='1'
export LOOPCHANGE_IEEDGE_INNER='1'

export LOOPCHANGE_INSIDEERROR_OUTER='1'
export LOOPITERSKIP_INSIDEERROR_OUTER='1'

export LOOPCHANGE_INSIDEERROR_INNER='1'
export LOOPITERSKIP_INSIDEERROR_INNER='1'

export LOOPCHANGE_EDGEERROR_1='1'
export LOOPITERSKIP_EDGEERROR_1='1'

export LOOPCHANGE_EDGEERROR_2='1'
export LOOPITERSKIP_EDGEERROR_2='1'






export LOOPCHANGE_PFUPDATE='2'

export LOOPCHANGE_IEINSIDE_OUTER='2'

export LOOPCHANGE_IEINSIDE_INNER='2'

export LOOPCHANGE_IEEDGE_OUTER='2'

export LOOPCHANGE_IEEDGE_INNER='2'

export LOOPCHANGE_INSIDEERROR_OUTER='2'

export LOOPCHANGE_INSIDEERROR_INNER='2'

export LOOPCHANGE_EDGEERROR_1='2'

export LOOPCHANGE_EDGEERROR_2='2'



# parsecmgmt -a build -p parsec.bodytrack   -c gcc-serial

# ../obj/amd64-linux.gcc/TrackingBenchmark/bodytrack ../inputs/sequenceB_1 4 1 1000 5 0 1 > log5.txt
# ./obj/amd64-linux.gcc-serial/TrackingBenchmark/bodytrack ../inputs/sequenceB_1 4 1 1000 5 0 1 > log5.txt