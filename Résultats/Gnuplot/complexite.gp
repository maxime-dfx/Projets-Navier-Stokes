# ==========================================
# CONFIGURATION
# ==========================================
set terminal pngcairo size 1200,800 enhanced font 'Verdana,12'
set output 'Résultats/Courbes/comparaison_complexite.png'

set title "Analyse de Complexité Algorithmique (Coût vs N)" font ",16"
set xlabel "Résolution N (Nb points par axe)" font ",14"
set ylabel "Temps de Calcul [s] (log)" font ",14"

set grid lc rgb "#DDDDDD" lt 1 
set logscale xy
set format y "%.1e"
set key top left box opaque
set offset graph 0.1, graph 0.1, graph 0.1, graph 0.1

# ==========================================
# REGRESSIONS
# ==========================================
# Modèles : Temps = a * N^b (N = 1/h)
f_euler(x) = a1 * x**b1; a1=1e-6; b1=4;
f_rk2(x)   = a2 * x**b2; a2=1e-6; b2=4;
f_rk4(x)   = a3 * x**b3; a3=1e-6; b3=4;

fit f_euler(x) "results/Convergence/convergence_ExplicitEuler.dat" using (1.0/$1):4 via a1, b1
fit f_rk2(x)   "results/Convergence/convergence_RungeKutta2.dat" using (1.0/$1):4 via a2, b2
fit f_rk4(x)   "results/Convergence/convergence_RungeKutta4.dat" using (1.0/$1):4 via a3, b3

t_euler = sprintf("Euler : O(N^{%.2f})", b1)
t_rk2   = sprintf("RK2   : O(N^{%.2f})", b2)
t_rk4   = sprintf("RK4   : O(N^{%.2f})", b3)

# ==========================================
# TRACÉ SUPERPOSÉ
# ==========================================
plot \
  "results/Convergence/convergence_ExplicitEuler.dat" u (1.0/$1):4 w p pt 7 ps 1.5 lc rgb "red" t "Euler", \
  f_euler(x) w l lw 2 lc rgb "red" dt 2 t t_euler, \
  "results/Convergence/convergence_RungeKutta2.dat" u (1.0/$1):4 w p pt 9 ps 1.5 lc rgb "blue" t "RK2", \
  f_rk2(x)   w l lw 2 lc rgb "blue" dt 2 t t_rk2, \
  "results/Convergence/convergence_RungeKutta4.dat" u (1.0/$1):4 w p pt 5 ps 1.5 lc rgb "green" t "RK4", \
  f_rk4(x)   w l lw 2 lc rgb "green" dt 2 t t_rk4