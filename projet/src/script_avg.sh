#!/bin/bash

CSV_FILE="results.csv"
ENV_VEAR_NAME="OMP_NUM_THREADS"
EXECUTABLE_PATH="./simulation.exe"
ENV_VALUES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)
N_RUNS=5 # Número de execuções para média

# Cabeçalho do CSV
echo "ENV_VALUE,n_iterations,avg_global_time,avg_time_update,avg_time_affichage,avg_time_for" > "$CSV_FILE"

for val in "${ENV_VALUES[@]}"; do
  # Variáveis para acumular resultados
  sum_global=0
  sum_update=0
  sum_affichage=0
  sum_for=0
  n_iter=0

  for ((run=1; run<=N_RUNS; run++)); do
    # Executa o programa e captura a saída
    full_output=$(export $ENV_VEAR_NAME=$val && $EXECUTABLE_PATH -n 300 -v $val)
    echo "$full_output"
    # Extrai valores (usando awk para garantir números inteiros)
    curr_iter=$(echo "$full_output" | awk '/n_iterations/ {print $2}')
    curr_global=$(echo "$full_output" | awk '/global_time/ {print $2}')
    curr_update=$(echo "$full_output" | awk '/avg_time_update/ {print $2}')
    curr_for=$(echo "$full_output" | awk '/avg_time_for/ {print $2}')
    curr_affichage=$(echo "$full_output" | awk '/avg_time_affichage/ {print $2}')

    # Validação básica dos resultados
    if [ -z "$curr_global" ] || [ -z "$curr_update" ] || [ -z "$curr_affichage" ]; then
      echo "Erro na execução $run com $ENV_VEAR_NAME=$val. Resultados ignorados."
      continue
    fi

    # Acumula para cálculo da média
    sum_global=$(echo "$sum_global + $curr_global" | bc)
    sum_update=$(echo "$sum_update + $curr_update" | bc)
    sum_affichage=$(echo "$sum_affichage + $curr_affichage" | bc)
    sum_for=$(echo "$sum_for + $curr_for" | bc)
    
    # Guarda o n_iterations da primeira execução válida
    if [ "$n_iter" -eq 0 ]; then
      n_iter=$curr_iter
    fi
  done

  # Calcula médias (usando awk para precisão decimal)
  avg_global=$(awk -v s=$sum_global -v n=$N_RUNS 'BEGIN {printf "%.10f", s/n}')
  avg_update=$(awk -v s=$sum_update -v n=$N_RUNS 'BEGIN {printf "%.10f", s/n}')
  avg_affichage=$(awk -v s=$sum_affichage -v n=$N_RUNS 'BEGIN {printf "%.10f", s/n}')
  avg_for=$(awk -v s=$sum_for -v n=$N_RUNS 'BEGIN {printf "%.10f", s/n}')

  # Adiciona ao CSV
  echo "$val,$n_iter,$avg_global,$avg_update,$avg_affichage,$avg_for" >> "$CSV_FILE"
done

echo "CSV gerado: $CSV_FILE"
