# 🌊 Solveur Navier-Stokes 2D Incompressible

Ce projet est un code de calcul scientifique performant (CFD) résolvant les équations de Navier-Stokes pour un fluide incompressible en 2D. Il est développé en **C++17 Moderne** et met l'accent sur la modularité, la performance et la rigueur mathématique.

L'architecture repose sur une grille décalée **MAC (Marker-and-Cell)** et utilise la méthode de **Projection de Chorin** pour le couplage vitesse-pression.

---

## 🚀 Fonctionnalités Clés

### 🧠 Cœur Numérique
* **Discrétisation Spatiale :** Différences Finies sur grille décalée (Staggered Grid).
* **Intégration Temporelle :**
    * Euler Explicite (Ordre 1).
    * **Runge-Kutta 2** (Ordre 2).
    * **Runge-Kutta 4** (Ordre 4, haute stabilité).
* **Solveur de Pression :** Solveur direct Cholesky (`SimplicialLLT`) via la bibliothèque **Eigen** pour une résolution rapide de l'équation de Poisson.
* **Transport (Advection) :** Schéma **Upwind** (décentré amont) pour la stabilité numérique.

### ⚙️ Configuration & I/O
* **Configuration :** Fichiers `.toml` pour paramétrer la simulation sans recompiler (Géométrie, Temps, Physique, Conditions limites).
* **Conditions aux Limites :** Dirichlet, Neumann, Profil de Poiseuille, Obstacles cylindriques arbitraires.
* **Sorties :**
    * Format **VTK** (`.vtk`) pour visualisation 3D sous Paraview (Pression, Vitesse, Vorticité).
    * Sondes ponctuelles (`.dat`) pour analyse temporelle.

### 📊 Analyse & Validation
* **Validation :** Calcul automatique de l'erreur $L^2$ par rapport à la solution analytique de Poiseuille.
* **Convergence :** Script d'étude de convergence en temps et en espace.
* **Physique :** Simulation de l'allée de tourbillons de Von Kármán et analyse spectrale (Nombre de Strouhal).

---

## 📦 Prérequis

Pour compiler et lancer ce projet, vous avez besoin de :

* **Compilateur C++ :** Compatible C++17 (ex: `g++` >= 8, `clang++`).
* **Eigen 3 :** Bibliothèque d'algèbre linéaire (Header-only).
* **Gnuplot :** Pour visualiser les courbes de convergence et les sondes.
* **Paraview :** Pour visualiser les champs de vitesse/pression en 2D/3D.

---

## 🛠️ Installation & Compilation

1.  **Cloner le dépôt :**
    ```bash
    git clone [https://github.com/votre-repo/navier-stokes-2d.git](https://github.com/votre-repo/navier-stokes-2d.git)
    cd navier-stokes-2d
    ```

2.  **Configuration des chemins :**
    Ouvrez le fichier `makefile` et ajustez les variables si nécessaire (si vos librairies ne sont pas dans les chemins standards) :
    ```makefile
    EIGEN_INCLUDE = ${HOME}/libraries/eigen/
    TOML_INCLUDE  = ${HOME}/libraries/
    ```

3.  **Compilation :**
    Le projet utilise un Makefile optimisé (`-O3 -march=native`).
    ```bash
    make        # Compile l'exécutable
    ```

---

## 🖥️ Utilisation (Scénarios)

Le `makefile` propose des cibles prêtes à l'emploi pour les cas tests classiques.

### 1. Cas Von Kármán (Allée tourbillonnaire) 🌪️
Simule un écoulement autour d'un cylindre. C'est le cas le plus visuel.
```bash
make karman