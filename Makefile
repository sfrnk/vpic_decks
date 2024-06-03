# Makefile for building final exe and cleaning directory

# path to installed compiler wrapper
VPIC = /anvil/projects/x-phy230179/wham/build/20240131.nsmeb/bin/vpic
# name of input deck file
DECK = deck/WHAM-phase2-3d.cxx
# name of platform to build for
PLAT = Linux

# compile final exe
all:
	$(VPIC) $(DECK)

clean:
	-$(RM) $(basename $(notdir $(DECK))).$(PLAT)
	-$(RM) -r global.vpc
	-$(RM) -r info.dat info.bin info.hdf5
	-$(RM) -r energies.dat particles_count.dat
	-$(RM) -r restart*
	-$(RM) -r particles
	-$(RM) -r fields
	-$(RM) -r hydro
	-$(RM) -r data

# example mpirun commands. --tpp is number of vpic pipelines. Each
# pipeline has domain duplication of the accumulator array.

# openmpi examples
# one rank on a node, bound to all available cores.
# mpirun --bind-to none --report-bindings -n 1 ./harris.Linux 1 1 --tpp 16

# two ranks on a node with two sockets. ranks bound to sockets.
# mpirun --map-by socket --bind-to socket --report-bindings -n 2 ./harris.Linux  1 1 --tpp 16

# 4 ranks on a node with 4 numa nodes. ranks bound to numa nodes.
# mpirun --map-by numa --bind-to numa --report-bindings -n 4 ./harris.Linux 1 1 --tpp 8

# aprun examples on four nodes. see man aprun for flags.
# 1 rank/core = 128 ranks
# aprun -n 128 ./reconnection.Linux
# 1 rank/core, 2 threads/core
# aprun -n 128 -d 2 -j 2 ./reconnection.Linux --tpp 2
# 8 ranks, 2 ranks/node, 1 rank/numa, 32 CPUs/rank, 2 CPUs/compute unit
# aprun -n 8 -N 2 -S 1 -d 32 -j 2 ./reconnection.Linux --tpp 32
