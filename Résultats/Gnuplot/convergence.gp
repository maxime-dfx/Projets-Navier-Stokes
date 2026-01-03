# ==========================================
# CONFIGURATION
# ==========================================
set terminal pngcairo size 1200,800 enhanced font 'Verdana,12'
set output 'Résultats/Courbes/comparaison_convergence.png'

set title "Comparaison des Ordres de Convergence (L2)" font ",16"
set xlabel "Pas de maillage h (log)" font ",14"
set ylabel "Erreur L2 (log)" font ",14"

# Grille et Echelles
set grid lc rgb "#DDDDDD" lt 1 
set logscale xy
set format x "%.0e"; set format y "%.0e"
set key bottom right box opaque
set offset graph 0.1, graph 0.1, graph 0.1, graph 0.1

# ==========================================
# REGRESSIONS (Calcul des pentes)
# ==========================================
# Modèles : Erreur = a * h^b
f_euler(x) = a1 * x**b1; a1=1; b1=1;
f_rk2(x)   = a2 * x**b2; a2=1; b2=1;
f_rk4(x)   = a3 * x**b3; a3=1; b3=2;

# Fits (On ignore les erreurs si un fichier manque)
fit f_euler(x) "results/Convergence/convergence_ExplicitEuler.dat" using 1:2 via a1, b1
fit f_rk2(x)   "results/Convergence/convergence_RungeKutta2.dat" using 1:2 via a2, b2
fit f_rk4(x)   "results/Convergence/convergence_RungeKutta4.dat" using 1:2 via a3, b3

# Titres dynamiques
t_euler = sprintf("Euler (Ordre %.2f)", b1)
t_rk2   = sprintf("RK2   (Ordre %.2f)", b2)
t_rk4   = sprintf("RK4   (Ordre %.2f)", b3)

# ==========================================
# TRACÉ SUPERPOSÉ
# ==========================================
plot \
  "results/Convergence/convergence_ExplicitEuler.dat" u 1:2 w p pt 7 ps 1.5 lc rgb "red" t "Euler Data", \
  f_euler(x) w l lw 2 lc rgb "red" dt 2 t t_euler, \
  "results/Convergence/convergence_RungeKutta2.dat" u 1:2 w p pt 9 ps 1.5 lc rgb "blue" t "RK2 Data", \
  f_rk2(x)   w l lw 2 lc rgb "blue" dt 2 t t_rk2, \
  "results/Convergence/convergence_RungeKutta4.dat" u 1:2 w p pt 5 ps 1.5 lc rgb "green" t "RK4 Data", \
  f_rk4(x)   w l lw 2 lc rgb "green" dt 2 t t_rk4