#!/usr/bin/env bash
set -euo pipefail

########################  PARÁMETROS AJUSTABLES  ########################
LISTA="text.txt"   
THRESHOLD="0.9"    
REPS=1             

THREADS=8          
BLOCK_SIZE=128     

MH_SIZE_ARR=(512) 
EPS=1e-6          
#########################################################################

CPU_BINARY="./build/selection"
GPU_BINARY="./build/selection_cuda"

OUT="comparacion_cpu_gpu.csv"
echo "cfg,card1,card2,sim_cpu,sim_gpu,diff" > "$OUT"

#########################################################################
# f_run  → corre un binario, añade llave «cardA_cardB» y ordena          #
#########################################################################
f_run () {            # $1 = comando completo + args
  local tmp
  tmp=$(mktemp)
  "$@" | awk '{key=$1"_"$2; print key,$0}' | \
        LC_ALL=C sort -k1,1 -S1G --parallel=2 > "$tmp"
  echo "$tmp"
}

#########################################################################
#                        BUCLE PRINCIPAL                                #
#########################################################################
for m in "${MH_SIZE_ARR[@]}"; do
  for r in $(seq 1 "$REPS");    do

    cpu_out=$(f_run "$CPU_BINARY" -l "$LISTA" -t "$THREADS" -h "$THRESHOLD" -a "$m" -c "smh_a")
    gpu_out=$(f_run "$GPU_BINARY" -l "$LISTA" -b "$BLOCK_SIZE" -h "$THRESHOLD" -a "$m" )

    # Une ambas salidas por la clave y calcula la diferencia
    join -1 1 -2 1 "$cpu_out" "$gpu_out" |                        \
    awk -v cfg="t${THREADS}_b${BLOCK_SIZE}_m${m}_r${r}" -v eps="$EPS" '{
          sim_cpu=$4; sim_gpu=$7;
          diff=(sim_cpu>sim_gpu?sim_cpu-sim_gpu:sim_gpu-sim_cpu);
          if (diff<eps) diff=0;
          print cfg","$2","$3","sim_cpu","sim_gpu","diff
    }' >> "$OUT"

    rm -f "$cpu_out" "$gpu_out"
  done
done

echo "Comparación completada: $(wc -l < "$OUT") líneas en '$OUT'"
