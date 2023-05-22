#!/bin/bash
 
# On an example HPC server:
# . ~/.bashrc
# module purge
# module load toolchain/intel/2018.5.274
# module load numlib/ScaLAPACK/2.1.0-gompi-2021b-fb

make clean
make all

# Make a checkpoints folder if one does not exist
if [ -d checkpoints/ ]; then
  mkdir checkpoints/
fi

# Get the most recent checkpoint if there is one
#checkpointfile=$(ls -rt checkpoints/net.1* | tail -1)
checkpointfile=$(ls checkpoints/net.1* | sort --version-sort | grep -v 'net.1-0' | tail -1)
if [ ! -z "$checkpointfile" ]; then
  cp $checkpointfile ./net.1
  cp $checkpointfile ./net.1.restart
fi
rm checkpoints/*

# Get the upper and lower bounds of the training
# set given this neural network
./lookAtBounds.x

# Train the NN
./training.x

# Measure the error in both energies and forces
# for the most recent (non-rejected) model
checkpointfile=$(ls checkpoints/net.1* | sort --version-sort | grep -v 'net.1-0' | tail -1)
cp $checkpointfile ./net.1.restart
./testEnergyAndOrForces.x

