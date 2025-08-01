#!/bin/bash

LISTA="test_influeza_filelist.txt"
#LISTA="test.txt"
THRESHOLD="0.9"
REPS=1
THREADS_ARR=(1 8)
#MH_SIZE_ARR=(16 32 64 128 256 512 1024)
MH_SIZE_ARR=(512)

CPU_BINARY="./build/time_smh"
GPU_BINARY="./build/time_smh_cuda"

LOG="experimento_smh_comparativo.csv"
echo "impl,threads,mh_size,rep,criterio,tiempo" > $LOG

# ----------- CPU -----------
for T in "${THREADS_ARR[@]}"; do
  for M in "${MH_SIZE_ARR[@]}"; do
    for REP in $(seq 1 $REPS); do
      OUTPUT=$($CPU_BINARY -l $LISTA -t $T -h $THRESHOLD -m $M)
      echo "$OUTPUT" | grep ';build_smh;'   | awk -F';' -v t=$T -v m=$M -v r=$REP '{print "cpu,"t","m","r",build_smh,"$4}'   >> $LOG
      echo "$OUTPUT" | grep ';smh_a;'       | awk -F';' -v t=$T -v m=$M -v r=$REP '{print "cpu,"t","m","r",smh_a,"$4}'       >> $LOG
      echo "$OUTPUT" | grep ';CB+smh_a;'    | awk -F';' -v t=$T -v m=$M -v r=$REP '{print "cpu,"t","m","r",CB+smh_a,"$4}'    >> $LOG
    done
  done
done

# ----------- GPU (CUDA) -----------
for M in "${MH_SIZE_ARR[@]}"; do
    for REP in $(seq 1 $REPS); do
        OUTPUT=$($GPU_BINARY -l $LISTA -h $THRESHOLD -m $M)
        echo "$OUTPUT" | grep ';build_smh;'   | awk -F';' -v m=$M -v r=$REP '{print "gpu,"t","m","r",build_smh,"$4}'   >> $LOG
        echo "$OUTPUT" | grep ';smh_a;'       | awk -F';' -v m=$M -v r=$REP '{print "gpu,"t","m","r",smh_a,"$4}'       >> $LOG
        echo "$OUTPUT" | grep ';CB+smh_a;'    | awk -F';' -v m=$M -v r=$REP '{print "gpu,"t","m","r",CB+smh_a,"$4}'    >> $LOG
    done
done

echo "Listo, resultados en $LOG"
