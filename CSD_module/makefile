#*********************************************************************
#*  Author: Rimple Sandhu
#*  Email: rimple_sandhu@outlook.com
#*********************************************************************

SRC_DIR = src
DATA_DIR = data
MESH_OUT_DIR = output/mesh
RESTART_OUT_DIR = output/restart
EIGEN_OUT_DIR = output/eigen
TRANSIENT_OUT_DIR = output/transient

GEO_FILE = WING
EIGEN_SIF_FILE = EIGEN
TRANSEINT_SIF_FILE = TRANSIENT
STATIC_SIF_FILE = STATIC
MY_FORCING_FILE = AERO_PRESSURE
MY_FUNCTIONS_FILE = common

mesh1:
	gmsh -3 $(SRC_DIR)/$(GEO_FILE).GEO -optimize -o $(MESH_OUT_DIR)/$(GEO_FILE).msh -order 1

mesh2:
	gmsh -3 $(SRC_DIR)/$(GEO_FILE).GEO -optimize -o $(MESH_OUT_DIR)/$(GEO_FILE).msh -order 2

clean:
	rm -r $(TRANSIENT_OUT_DIR)/AERO_PRESSURE.dSYM;
	rm -f *.out fort* $(EIGEN_OUT_DIR)/* $(TRANSIENT_OUT_DIR)/*  $(RESTART_OUT_DIR)/*t.str $(RESTART_OUT_DIR)/initial.str debug* $(EIGEN_OUT_DIR)/../.DS_Store $(DATA_DIR)/X*.dat $(DATA_DIR)/PARAM*.dat $(DATA_DIR)/F_*.dat ../*.mod;

cleanSPLINE:
	rm -f $(DATA_DIR)/MAT_H.dat 

cleanSPLINEB:
	rm -f $(DATA_DIR)/MAT_HB.dat 

show:
	gmsh $(MESH_OUT_DIR)/$(GEO_FILE).msh &

show_g:
	gmsh $(SRC_DIR)/$(GEO_FILE).geo &

to_elmer:
	ElmerGrid 14 2 $(MESH_OUT_DIR)/$(GEO_FILE).msh -autoclean -out $(MESH_OUT_DIR)

eigen:
	ElmerSolver $(SRC_DIR)/$(EIGEN_SIF_FILE).SIF > $(EIGEN_OUT_DIR)/eigen.log 2>&1

static:
	ElmerSolver $(SRC_DIR)/$(STATIC_SIF_FILE).SIF > $(RESTART_OUT_DIR)/static.log 2>&1

transient:
	ElmerSolver $(SRC_DIR)/$(TRANSEINT_SIF_FILE).SIF >> $(TRANSIENT_OUT_DIR)/transient.log 2>> $(TRANSIENT_OUT_DIR)/error.log

transient_zero:
	ElmerSolver $(SRC_DIR)/$(TRANSEINT_SIF_FILE)_ZERO.SIF >> $(TRANSIENT_OUT_DIR)/transient_zero.log 2>> $(TRANSIENT_OUT_DIR)/error_zero.log

compile:
	elmerf90 $(SRC_DIR)/$(MY_FUNCTIONS_FILE).F $(SRC_DIR)/$(MY_FORCING_FILE).F -o $(TRANSIENT_OUT_DIR)/$(MY_FORCING_FILE) > $(TRANSIENT_OUT_DIR)/compile.log 2>&1;
	rm $(MY_FUNCTIONS_FILE).mod;


