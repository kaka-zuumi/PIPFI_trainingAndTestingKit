# Simplified PIP Neural Network Training and Testing

A set of fortran code that trains a neural network given (1) a training set, (2) a set of inputs, and (3) a neural network architecture. To the current date, only networks with exactly two hidden layers and one output layer (i.e., the energy) are supported. Energy and force errors can subsequently be tested.

## Getting Started

See https://doi.org/10.1063/5.0010104 for information on how PIP/FI neural networks can be used to model potential energy surfaces. The following repositories have fortran code that may be used as *input* for the training and testing: https://github.com/kjshao/FI, https://github.com/Chenrongjun/FI. Two examples are given here: BrClH2 (molecular formula ABC2) and BrCH5 (ABC5).

A neural network is defined by its (1) architecture, (2) inputs, and (3) parameters, where the parameters is usually a large number of adjustable weights and biases. Nearly all of this information is specified in a single file: "net.1". The header describes the architecture: the number of hidden layers, the number of neurons per layer, and the total number of parameters. The second section describes the inputs' bounds: the minimum and maximum value of each input (and one extra line for the output). The last section describes the parameters: The connection between each layer and the next neurons is given as a set of four indexes, followed by the value of each parameter in a single column. Finally, the way the inputs are defined (with respect to the coordinates) are specified in the PIPfile. 

### Prerequisites

Intel fortran with LAPACK is required.

```
blah blah blah
```

### Installing

After getting the preqrequisites ready, the code can be quickly downloaded, compiled, and executed.

First, download the github source code:
```
git clone https://github.com/kaka-zuumi/PIPFI_trainingAndTestingKit.git
```

Second, prepare the header of the GNUmakefile so as to point to the correct fortran package.
```
# On an example Linux desktop, these LIBs and INLCUDEs work:
MKLROOT = /opt/intel/oneapi/mkl/2023.1.0
LIB = -L${MKLROOT}/lib/intel64 -qmkl -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64
INCLUDE = -I${MKLROOT}/include/intel64/ilp64 -i8  -I"${MKLROOT}/include"
```

Third, make.
```
make clean; make all
```

## Running the tests

There are three executables for three different actions.

To figure out the minimum and maximum bounds for the training set given, execute:
```
./lookAtBounds.x
```

The input neural network can be adjusted in the file "net.1". Here, the bounds of the inputs (as seen above) can be adjusted. Following this, the network may be trained with:
```
./training.x
```

If there is no restart file "net.1.restart", then the neural network will use the architecture and bounds specified in "net.1" but initialize the weights from random values. Otherwise, the weights are read in from "net.1.restart"

At any time, a neural network may be tested with:
```
./testEnergyAndOrForces.x
```
where the predicted energies and forces are compared with the actual energies and forces in the training set (force predictions may be turned off if needed).
