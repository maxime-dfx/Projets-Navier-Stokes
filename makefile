# ==============================================================================
#  SOLVEUR NAVIER-STOKES 2D - MAKEFILE (MODERNE & VTK)
# ==============================================================================

CXX = g++

# --- DOSSIERS ---
SRC_DIR   = src
LIB_DIR   = lib
BUILD_DIR = build
BIN_DIR   = bin
DATA_DIR  = results

# --- FICHIERS ---
EXEC = $(BIN_DIR)/ns_solver_2d
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# --- CONFIGURATION BIBLIOTHEQUES ---
# Assurez-vous que le dossier 'eigen' contient bien le sous-dossier 'Eigen'
# ou les headers directement.
EIGEN_INCLUDE = ${HOME}/librairies/eigen/
TOML_INCLUDE  = ${HOME}/librairies/toml11/include/

INCLUDES      = -I$(LIB_DIR) -I$(EIGEN_INCLUDE) -I$(TOML_INCLUDE)

# --- FLAGS DE COMPILATION ---
# -std=c++17    : Requis pour std::filesystem et std::make_unique
# -O3           : Optimisation maximale
# -march=native : Optimise pour ton CPU (AVX, etc.)
# -Wall -Wextra : Avertissements complets
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra $(INCLUDES)

# --- FLAGS D'ÉDITION DE LIENS ---
# -lstdc++fs : Nécessaire sur GCC < 9 pour std::filesystem
LDFLAGS = -lstdc++fs

# ==============================================================================
# RÈGLES GÉNÉRALES
# ==============================================================================

.PHONY: all clean clean_all help mkdirs poiseuille stokes karman turbulent validation

all: $(EXEC)

help:
	@echo "--- SOLVEUR NS 2D (Sortie VTK) ---"
	@echo "  make poiseuille  : Lance le cas Poiseuille"
	@echo "  make karman      : Lance le cas Von Karman"
	@echo "  make stokes      : Lance le cas Stokes"
	@echo "  make turbulent   : Lance le cas Turbulent"
	@echo "  make validation  : Lance l'étude de convergence"
	@echo "----------------------------------"
	@echo "  make clean       : Nettoie les objets (.o)"
	@echo "  make clean_all   : Nettoie tout (exécutable + résultats)"

# Édition des liens
$(EXEC): $(OBJS) | mkdirs
	@echo ">> [LINK] Liaison de l'exécutable..."
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Compilation des sources
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | mkdirs
	@echo ">> [CXX]  Compilation de $<..."
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
	@echo ">>> Lancement Simulation : St	okes"
	./$(EXEC) input/stokes.toml

karman: $(EXEC)
	@echo ">>> Lancement Simulation : Von Karman"
	./$(EXEC) input/karman.toml

turbulent: $(EXEC)
	@echo ">>> Lancement Simulation : Turbulent"
	./$(EXEC) input/turbulent.toml

validation: $(EXEC)
	@echo ">>> Lancement Validation : Poisseuille"
	./$(EXEC) input/convergence.toml

# ==============================================================================
# OUTILS & NETTOYAGE
# ==============================================================================

plots:
	@echo ">>> Génération des graphiques Gnuplot..."
	@mkdir -p Résultats/Gnuplot
	@for file in Résultats/Gnuplot/*.gp; do \
		if [ -f "$$file" ]; then \
			echo "   >> Exécution de $$file"; \
			gnuplot "$$file"; \
		else \
			echo "   >> Aucun fichier .gp trouvé dans Résultats/Gnuplot"; \
		fi; \
	done
	
clean:
	@echo ">> Nettoyage des objets..."
	rm -rf $(BUILD_DIR)

clean_all: clean
	@echo ">> Nettoyage complet..."
	rm -rf $(BIN_DIR) $(DATA_DIR) *.png *.dat