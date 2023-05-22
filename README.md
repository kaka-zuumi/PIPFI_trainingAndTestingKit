# Simplified PIP Neural Network Training and Testing

A set of fortran code that trains a neural network given (1) a training set, (2) a set of inputs, and (3) a neural network architecture. To the current date, only networks with exactly two hidden layers and one output layer (i.e., the energy) are supported. Energy and force errors can subsequently be tested.

## Getting Started

See https://doi.org/10.1063/5.0010104 for information on how PIP/FI neural networks can be used to model potential energy surfaces. The following repositories have fortran code that may be used as *input* for the training and testing: https://github.com/kjshao/FI, https://github.com/Chenrongjun/FI. Two examples are given here: BrClH2 (molecular formula ABC2) and BrCH5 (ABC5).

### Prerequisites

Intel fortran with LAPACK is required. This may be installed from: https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html

### Installing

After getting the preqrequisites ready, the code can be quickly downloaded, compiled, and executed.

First, download the github source code:
```
git clone https://github.com/kaka-zuumi/PIPFI_trainingAndTestingKit.git
```

Second, prepare the header of the "GNUmakefile" so as to point to the correct fortran package.
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

### Running the tests

There are three executables for three different actions.

To figure out the minimum and maximum bounds for the training set given, execute:
```
./lookAtBounds.x
```

The resulting bounds are stored in the header section of "checkpoints/net.1-0".

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

## How is the neural network stored?

![Alt text](arch1.png?raw=true "Neural Network Architecture")

A neural network is defined by its (1) architecture, (2) inputs, and (3) parameters, where the parameters is usually a large number of adjustable weights and biases. Nearly all of this information is specified in a single network file: "net.1". The header describes the architecture: the number of hidden layers, the number of neurons per layer, and the total number of parameters. The second section describes the inputs' bounds: the minimum and maximum value of each input (and one extra line for the output). The last section describes the parameters: The connection between each layer and the next neurons is given as a set of four indexes, followed by the value of each parameter in a single column. Finally, the way the inputs are defined (with respect to the coordinates) are specified in the PIPfile. 

### The network file

For example, in "net.1.BrClH2", an example neural network is specified for the BrClH2 system:
```
           2
           7
          40
          49
        2379
  7.375507464055091E-071  0.666554315818088     
  1.930174954872492E-027  0.579718363531905     
  5.216195185551709E-071  0.479075556337930     
  1.667706802377780E-027  0.410968326894610     
  4.426156027365823E-036  0.302648560752083     
  3.109719783312207E-038  0.545271533681930     
  1.889943433156477E-071  0.153914592695499     
  0.000000000000000E+000   5.31782066167973     
           1           8           1           8
           1           8           9          16
           1           8          17          24
           1           8          25          32
           1           8          33          40
           1           8          41          48
           1           8          49          56
...
           9          49        2166        2206
           9          49        2207        2247
           9          49        2248        2288
           9          49        2289        2329
          50          99        2330        2379
 -1.403441280868372E-002
 -5.593424987956366E-002
 -4.792931618446141E-002
 -7.161430277231182E-002
...
```

The first single-column section is the neural network architecture. The first line indicates there are two hidden layers. The next line indicates there are 7 input neurons (this must corespond to the true number of inputs in the corresponding PIPfile). The next two numbers are the number of neurons in the two hidden layers (40 and 49). The last single-column number is the number of parameters.

The two-column section is the minimum and maximum values of each of the 7 input neurons and the minimum and maximum value of the single output neuron (the eneergy), for a total of 8 lines. If a training set is given, these can automatically be generated 

The four-column section details how each neuron from one layer connects to a neuron in the next layer; this can be automatically generated given an awk command in "README.txt".

The final single-column section is the value of each parameter in the neural network. If the header is correct, there should be 2379 lines in this section.

### The PIP file

The PIP file specifies how the inputs are made. For the BrClH2 system, an example is given in "obj/BrClH2.f90":
```
module pipvariables
implicit none

! Parameters required for interfacing with
! the "pes" and "training" software...
! must match the subroutines below

integer,parameter :: Natoms = 4  ! Number of atoms
integer,parameter :: Np     = 7  ! Number of polynomials

end module pipvariables

! The PIPs/FIs are defined here... switch
! the content of this subroutine with
! whatever PIPs/FIs you want to use
subroutine get_p(r,p)
implicit none
real(kind=8),intent(in) :: r(6)
real(kind=8),intent(out) :: p(7)
p(1)=r(5)+r(4)
p(2)=r(3)+r(2)
p(3)=r(5)**2+r(4)**2
p(4)=r(3)**2+r(2)**2
p(5)=r(5)*r(3)+r(4)*r(2)
p(6)=r(6)
p(7)=r(1)
end subroutine get_p

...
```

The number of polynomials described here (7) must match the number of input neurons specified in the network file above (7). The polynomials are specified in terms of the pairwise atomic distances, which are a function of the atomic coordinates. So for example, the system above with 4 atoms will have 4x3/2 = 6 atomic distances. Symmetrizing the polynomials is not a trivial task and so these should be deliberately chosen. Whichever PIP/FI file is used, the fortran code must be compiled into the executable. An example is given in "GNUmakefile":

```
# Specify the file that has the input PIPs/FIs:
PIPfile = obj/BrClH2.f90
dPIPfile = obj/BrClH2-derivatives.f90
```

If forces are wanted, the derivatives of the PIPs/FIs are also necessary. The derivatives of a new PIP/FI file can automatically be generated acccording to the command in "README.txt". In this example, they are then stored in "obj/BrClH2-derivatives.f90".

