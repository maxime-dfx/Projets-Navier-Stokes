# Solveur Navier-Stokes 2D Incompressible

Ce projet implémente un solveur de mécanique des fluides (Navier-Stokes) sur une grille MAC (Marker-and-Cell) en C++. Il utilise la méthode de projection de Chorin pour la pression et gère des obstacles arbitraires.

## 🚀 Fonctionnalités
- **Schéma numérique :** Différences finies sur grille décalée (Staggered Grid).
- **Temps :** Euler Explicite.
- **Solveur Linéaire :** Eigen (SimplicialLLT pour Poisson).
- **Visualisation :** Génération automatique de scripts Gnuplot et GIFs animés.
- **Traceurs :** Système de particules lagrangiennes pour visualiser l'écoulement.

## 📦 Prérequis
- Compilateur C++ (g++) supportant C++17.
- **Eigen** (Bibliothèque d'algèbre linéaire).
- **Gnuplot** (Pour la visualisation).
- **Make** (Pour la compilation).

## 🛠️ Compilation et Utilisation
Le projet dispose d'un `makefile` automatisé.

### 1. Compilation simple
```bash
make