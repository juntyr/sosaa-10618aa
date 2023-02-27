#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Check if CHEM and CASE are defined
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# define \n
#
#
# endef
#
# ifndef CHEM
# $(error Usage:${\n}$$ make CHEM=chemistry_folder_name CASE=case_name${\n})
# endif
#
# ifndef CASE
# $(error Usage:${\n}$$ make CHEM=chemistry_folder_name CASE=case_name${\n})
# endif

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Include M_DEF
# You need to create your own M_DEF to fit for your system.
# An example of M_DEF is shown in M_DEF_template.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
IN = M_DEF
include $(IN)

# Export variables to all sub make commands, these variables are not accessable outside make command
export IN
export CASE  # case name
export CHEM  # chemistry scheme folder

CHEMISTRY_LIB   = -lchemistry
AEROSOL_LIB     = -laerosol
ACDC_LIB        = -lacdc
DATA_LIB        = -ldata
MEGAN_LIB       = -lmegan
DEPOSITION_LIB  = -ldeposition
METEOROLOGY_LIB = -lmeteorology

OBJ = $(NETCDF_LIBS) -L$(LIB) $(CHEMISTRY_LIB) $(ACDC_LIB) $(DATA_LIB) $(AEROSOL_LIB) $(MEGAN_LIB) $(DEPOSITION_LIB) $(METEOROLOGY_LIB)

OBJ_BOX = -L$(LIB) $(CHEMISTRY_LIB)

sosaa_1d: sosaa.exe

sosaa_box: sosaa_box.exe

sosaa.exe: preprocess chemistry_lib acdc_lib aerosol_lib data_lib megan_lib deposition_lib meteorology_lib
	$(LINK) $(OPT) $(SYS_INCL) $(LIB_INCL) $(NETCDF_INCL) -o $@ MAIN/Sosa.f90 $(OP) $(OBJ) $(CL)
	@mv -v sosaa.exe $(CASE_DIR)/$(CASE)/$(ALT_NAME)  # move the executable to the case folder
	@cp -nv $(INIT_FILE) $(CASE_DIR)/$(CASE)/  # copy the init file from the template if not existed

sosaa_box.exe: preprocess chemistry_lib
	$(LINK) $(OPT) $(SYS_INCL) $(LIB_INCL) -o $@ MAIN/sosaa_box.f90 $(OP) $(OBJ_BOX) $(CL)
	@mv -v sosaa_box.exe $(CASE_DIR)/$(CASE)/  # move the executable to the case folder

preprocess:
	@mkdir -p $(LIB)  # Create the folder $(LIB) if it does not exist.
	@mkdir -p $(CASE_DIR)/$(CASE)  # Create the case folder if it does not exist.

# Notice: Every time you link a new chemistry scheme, you need to 'make clean'.
chemistry_lib:
	@cd CHEMISTRY; make -f Make_src
	@echo "====================================================================================================="

aerosol_lib:
	@cd AEROSOL; make -f Make_src
	@echo "====================================================================================================="

data_lib:
	@cd DATA; make -f Make_src
	@echo "====================================================================================================="

megan_lib:
	@cd MEGAN; make -f Make_src
	@echo "====================================================================================================="

deposition_lib:
	@cd DEPOSITION; make -f Make_src
	@echo "====================================================================================================="

meteorology_lib:
	@cd METEOROLOGY; make -f Make_src
	@echo "====================================================================================================="

acdc_lib:
	@cd ACDC_NEW; make -f Make_src
	@echo "====================================================================================================="
.PHONY: kpp clean cleanchem cleanall

# Use kpp to generate second files for a chemistry scheme
# Use $$ to let the shell see one $
# Use \ to execute the commands in one child process
KPP_DIR=$(CODE_DIR)/tools/kpp
SEC_DIR=$(CHEMALL_DIR)/$(CHEM)
kpp:
	@cd $(SEC_DIR) ; \
	export KPP_HOME=$(KPP_DIR) ; \
	echo 'Now KPP_HOME is temporarily set to ' $$KPP_HOME; \
	$(KPP_DIR)/bin/kpp second.kpp

clean:
	@rm -fv $(LIB)/*.a
	@rm -fv $(LIB)/*.mod
	@rm -fv *.exe
	@rm -fv CHEMISTRY/*.o CHEMISTRY/*.mod
	@rm -fv AEROSOL/*.o AEROSOL/*.mod
	@rm -fv DATA/*.o DATA/*.mod
	@rm -fv MEGAN/*.o MEGAN/*.mod
	@rm -fv DEPOSITION/*.o DEPOSITION/*.mod
	@rm -fv METEOROLOGY/*.o METEOROLOGY/*.mod
	@rm -fv MAIN/*.o MAIN/*.mod
	@rm -fv ACDC_NEW/*.o ACDC_NEW/*.mod

cleanchem:
	@rm -fv $(CHEMALL_DIR)/$(CHEM)/LIB/*.a
	@rm -fv $(CHEMALL_DIR)/$(CHEM)/LIB/*.o
	@rm -fv $(CHEMALL_DIR)/$(CHEM)/LIB/*.mod

cleanall: clean cleanchem
