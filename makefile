# ==============================================================================
#  SOLVEUR NAVIER-STOKES 2D - MAKEFILE (MODERNE & VTK)
# ==============================================================================

CXX = g++

# --- DOSSIERS ---
SRC_DIR = src
LIB_DIR = lib
BUILD_DIR = build
BIN_DIR = bin
DATA_DIR = results

# --- FICHIERS ---
EXEC = $(BIN_DIR)/ns_solver_2d
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# --- CONFIGURATION ---
# Ajuste ces chemins selon ton environnement
EIGEN_INCLUDE = ${HOME}/libraries/eigen/
TOML_INCLUDE  = ${HOME}/libraries/

# Options de compilation :
# -std=c++17 : Standard moderne requis
# -O3 : Optimisation maximale (mieux que -O2 pour la CFD)
# -march=native : Optimise pour ton processeur spécifique
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -I$(LIB_DIR) -I$(EIGEN_INCLUDE) -I$(TOML_INCLUDE)

# ==============================================================================
# RÈGLES GÉNÉRALES
# ==============================================================================

.PHONY: all clean clean_all help mkdirs poiseuille stokes karman turbulent check

all: $(EXEC)

help:
	@echo "--- SOLVEUR NS 2D (Sortie VTK) ---"
	@echo "  make poiseuille  : Lance le cas Poiseuille"
	@echo "  make karman      : Lance le cas Von Karman"
	@echo "  make stokes      : Lance le cas Stokes"
	@echo "  make turbulent   : Lance le cas Turbulent"
	@echo "  make check       : Lance l'étude de convergence Python"
	@echo "----------------------------------"
	@echo "  make clean       : Nettoie les objets"
	@echo "  make clean_all   : Nettoie tout (exécutable + résultats)"

# Édition des liens
$(EXEC): $(OBJS) | mkdirs
	@echo "Liaison de l'exécutable..."
	$(CXX) $(CXXFLAGS) $^ -o $@

# Compilation des sources
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | mkdirs
	@echo "Compilation de $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

mkdirs:
	@mkdir -p $(BUILD_DIR) $(BIN_DIR) $(DATA_DIR)

# ==============================================================================
# SCENARIOS (SIMULATION UNIQUEMENT)
# ==============================================================================

poiseuille: $(EXEC)
	@echo ">>> Lancement Simulation : Poiseuille"
	./$(EXEC) input/poiseuille.toml

stokes: $(EXEC)
	@echo ">>> Lancement Simulation : Stokes"
	./$(EXEC) input/stokes.toml

karman: $(EXEC)
	@echo ">>> Lancement Simulation : Von Karman"
	./$(EXEC) input/karman.toml

turbulent: $(EXEC)
	@echo ">>> Lancement Simulation : Turbulent"
	./$(EXEC) input/turbulent.toml

# ==============================================================================
# OUTILS & NETTOYAGE
# ==============================================================================

# Lance l'étude de convergence (Python + Matplotlib)
check: $(EXEC)
	@echo ">>> Lancement Étude de Convergence..."
	./$(EXEC) input/poiseuille.toml check
	gnuplot résultats/validation.gp

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

clean_all: clean
	rm -rf $(DATA_DIR) *.png *.dat

clean_data: 
	rm -rf $(DATA_DIR) 