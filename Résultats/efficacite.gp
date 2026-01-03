# ==========================================
# CONFIGURATION
# ==========================================
set terminal pngcairo size 1200,800 enhanced font 'Verdana,12'
set output 'Résultats/comparaison_efficacite.png'

set title "Diagramme d'Efficacité : Quel schéma est le plus rentable ?" font ",16"
set xlabel "Temps de Calcul [s] (log)" font ",14"
set ylabel "Erreur L2 (log)" font ",14"

set grid lc rgb "#DDDDDD" lt 1 
set logscale xy
set format x "%.1e"; set format y "%.0e"
set key top right box opaque
set offset graph 0.1, graph 0.1, graph 0.1, graph 0.1

# ==========================================
# REGRESSIONS
# ==========================================
# Modèles : Erreur = a * Temps^b
f_euler(x) = a1 * x**b1; a1=1; b1=-0.25;
f_rk2(x)   = a2 * x**b2; a2=1; b2=-0.25;
f_rk4(x)   = a3 * x**b3; a3=1; b3=-0.25;

fit f_euler(x) "results/convergence_ExplicitEuler.dat" using 4:2 via a1, b1
fit f_rk2(x)   "results/convergence_RungeKutta2.dat" using 4:2 via a2, b2
fit f_rk4(x)   "results/convergence_RungeKutta4.dat" using 4:2 via a3, b3

t_euler = sprintf("Euler (Pente %.2f)", b1)
t_rk2   = sprintf("RK2   (Pente %.2f)", b2)
t_rk4   = sprintf("RK4   (Pente %.2f)", b3)

# ==========================================
# TRACÉ SUPERPOSÉ
# ==========================================
plot \
  "results/convergence_ExplicitEuler.dat" u 4:2 w p pt 7 ps 1.5 lc rgb "red" t "Euler", \
  f_euler(x) w l lw 2 lc rgb "red" dt 2 t t_euler, \
  "results/convergence_RungeKutta2.dat" u 4:2 w p pt 9 ps 1.5 lc rgb "blue" t "RK2", \
  f_rk2(x)   w l lw 2 lc rgb "blue" dt 2 t t_rk2, \
  "results/convergence_RungeKutta4.dat" u 4:2 w p pt 5 ps 1.5 lc rgb "green" t "RK4", \
  f_rk4(x)   w l lw 2 lc rgb "green" dt 2 t t_rk4