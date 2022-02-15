# par_nbody

## Correr individualmente
1. Los dos archivos .py (mainintegrator y parallelsketch) deben estar en el mismo directorio que la carpeta que contenga los sistemas
2. Correr como parallelsketch.py no_pert_nb nombre_ejemplo

## Correr en cluster

1. qsub task.sh 
2. times.sh nombre_ejemplo y zipatcluster.sh nombre_ejemplo se usan después de haber terminado la corrida (diagnóstico de tiempo total que alcanzó a evolucionar cada sistema / zipear lo importante para sacarlo con scp
3. sc_plots.py hace plots con info sobre conservación de energía de planetas ligados, energy_plots.py hace esos plots más bonitos pero no corre bien en el cluster
