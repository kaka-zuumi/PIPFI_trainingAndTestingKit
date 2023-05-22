# Simplified PIP Neural Network Training and Testing

A set of fortran code that trains a neural network given (1) a training set and (2) neural network architecture. To the current date, only networks with exactly two hidden layers and one output layer (i.e., the energy) are supported. Energy and force errors can subsequently be tested.

## Getting Started

See XXX for information on how PIP/FI neural networks can be used to model potential energy surfaces. The following repositories have fortran code that may be used as *input* for the training and testing. Two examples are given here: BrClH2 (molecular formula ABC2) and BrCH5 (ABC5).

### Prerequisites

Intel fortran with LAPACK is required.

```
blah blah blah
```

### Installing

After getting the preqrequisites ready, the code can be quickly downloaded, compiled, and executed.

First, download the github source code:
```
blah blah blah
```

Second, prepare the header of the GNUmakefile so as to point to the correct fortran package.
```
blah blah blah
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
