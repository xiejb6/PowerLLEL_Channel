#==========================================#
#    Makefile for PowerLLEL_Postprocess    #
#==========================================#
# Choose compilation options as needed
#   -D_SINGLE_PREC - use single-precision

HDF5_PATH = /WORK/app/hdf5/1.10.4/01-CF-14
DECOMP2D_PATH = ../decomp2d

#FC = mpiifort
FC = $(HDF5_PATH)/bin/h5pfc
FFLAGS = -O3 -xHost -ipo -fpp

OBJ = main_postproc.o \
      mod_param_postproc.o \
	  mod_bc_postproc.o \
	  mod_dataio_postproc.o \
	  mod_mesh_postproc.o \
	  mod_variables_postproc.o \
	  mod_vortex_postproc.o \
      mod_type.o \
	  mod_utils.o \
      mod_mpi.o \
      mod_hdf5.o \
      decomp_2d.o
TARGET = PowerLLEL_Postprocess PowerLLEL_Convert_Precision

all: $(TARGET)

PowerLLEL_Postprocess: $(OBJ)
	$(FC) $(FFLAGS) $^ -o $@

PowerLLEL_Convert_Precision: convert_precision.f90 mod_dataio_postproc.o mod_hdf5.o mod_mpi.o decomp_2d.o
	$(FC) $(FFLAGS) $^ -o $@

main_postproc.o: main_postproc.f90 mod_type.o mod_param_postproc.o mod_variables_postproc.o \
                 mod_mesh_postproc.o mod_dataio_postproc.o mod_bc_postproc.o mod_vortex_postproc.o \
				 mod_mpi.o mod_hdf5.o
	$(FC) $(FFLAGS) -c $<
mod_param_postproc.o: mod_param_postproc.f90 mod_type.o mod_mpi.o
	$(FC) $(FFLAGS) -c $<
mod_bc_postproc.o: mod_bc_postproc.f90 mod_type.o mod_mpi.o
	$(FC) $(FFLAGS) -c $<
mod_dataio_postproc.o: mod_dataio_postproc.f90 mod_mpi.o mod_hdf5.o
	$(FC) $(FFLAGS) -c $<
mod_mesh_postproc.o: mod_mesh_postproc.f90 mod_type.o mod_param_postproc.o mod_mpi.o
	$(FC) $(FFLAGS) -c $<
mod_variables_postproc.o: mod_variables_postproc.f90 mod_type.o mod_utils.o
	$(FC) $(FFLAGS) -c $<
mod_vortex_postproc.o: mod_vortex_postproc.f90 mod_type.o mod_param_postproc.o mod_variables_postproc.o \
                       mod_mesh_postproc.o mod_dataio_postproc.o mod_mpi.o
	$(FC) $(FFLAGS) -c $<
mod_type.o: ../mod_type.f90
	$(FC) $(FFLAGS) -c $<
mod_utils.o: ../mod_utils.f90 mod_type.o mod_mpi.o
	$(FC) $(FFLAGS) -c $<
mod_mpi.o: ../mod_mpi.f90 decomp_2d.o
	$(FC) $(FFLAGS) -c $<
# mod_hdf5.o: ../mod_hdf5.f90 mod_type.o
mod_hdf5.o: mod_hdf5_postproc.f90
	$(FC) $(FFLAGS) -c $< -o $@
decomp_2d.o: $(DECOMP2D_PATH)/decomp_2d.f90
	$(FC) $(FFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf *.o *.mod $(TARGET)