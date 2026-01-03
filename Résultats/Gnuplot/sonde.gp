# ==============================================================================
#  VISUALISATION SONDE - NAVIER STOKES 2D
# ==============================================================================

# 1. Configuration de la sortie (Image PNG haute qualité)
set terminal pngcairo size 1200,800 enhanced font 'Verdana,10'
set output 'Résultats/Courbes/resultats_sonde.png'

# 2. Configuration générale
set multiplot layout 2,1 title "Analyse de la Sonde (Von Kármán)" font ",14"
set grid
set key top right box opaque

# Définition des styles de ligne
set style line 1 lc rgb '#0060ad' lt 1 lw 1.5 pt 7 ps 0.5 # Bleu pro
set style line 2 lc rgb '#dd181f' lt 1 lw 2.0 # Rouge vif

# Nom du fichier de données (ADAPTE LE CHEMIN SI NECESSAIRE)
# Si ton fichier s'appelle autrement (ex: results/Karman_probe.dat), change-le ici.
FILES = system("ls results/Karman/Karman_Re120_Explicit_Euler_probe.dat | head -1")

# ==============================================================================
# GRAPHIQUE 1 : Signal Complet
# ==============================================================================
set title "1. Signal Temporel Complet (Transitoire + Établi)"
set xlabel "Temps [s]"
set ylabel "Vitesse Verticale V [m/s]"

plot FILES using 1:2 with lines ls 1 title "Vitesse V"

# ==============================================================================
# GRAPHIQUE 2 : Zoom sur les 20% derniers pourcents (Régime Établi)
# ==============================================================================
# On récupère le temps final (Tmax) via une commande système stats
stats FILES using 1 nooutput
T_max = STATS_max
T_start_zoom = T_max * 0.8

set title "2. Zoom sur le Régime Établi (Oscillations)"
set xrange [T_start_zoom : T_max]
set xlabel "Temps [s]"
set ylabel "Vitesse V [m/s]"

# On trace une ligne à 0 pour repère
set arrow from graph 0,first 0 to graph 1,first 0 nohead lc rgb "black" dt 2

plot FILES using 1:2 with lines ls 2 title "Oscillations Périodiques"

unset multiplot
print ">> Graphique généré : resultats_sonde.png"