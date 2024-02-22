# learning-based-homothetic-tube-mpc
This is the MATLAB code for the article
```
@article{gao2024learninghomothetic,
  title={Learning-based homothetic tube {MPC}},
  author={Yulong Gao, Shuhao Yan, Jian Zhou, and Mark Cannon,},
  year={2024},
  pages={},
  doi={ }
} 
```

Gao is with the Department of Electrical and Electronic Engineering, Imperial College London, UK. Yan is with the Department of Mathematics, University of Stuttgart, Germany. Zhou is with the  Department of Electrical Engineering, Linköping University, Sweden, and Cannon is with the  Department of Engineering Science, University of Oxford, UK.
## Packages for running the code
To run the code you need to install:

**MPT3**: https://www.mpt3.org/Main/Installation

**CasADi**: https://web.casadi.org/

**HSL Solver**: https://licences.stfc.ac.uk/product/coin-hsl

Note: Installing the HSL package can be a bit comprehensive, but the solvers just speed up the solutions. You can comment out the places where the HSL solver is used, i.e., options.ipopt.linear_solver = 'ma57', and just use the default linear solver of CasADi.

## Introduction to the files
The folder `Functions` contains basic functions like calculating the robust invariant set and so on.

The files `HomotheticUQMPC.m` and `RigidUQMPC.m` define the functions for learning-based homothetic MPC and learning-based rigid MPC, respectively. Here 'UQ' stands for 'uncertainty quantification', i.e., quantifying the uncertainties by learning-based approaches.


The file `systemmodel2V.m` defines all the basic parameters for the platooning case with 2 vehicles, and `main_2V.m` executes the controller in the 2-vehicle platooning scenario. The definition is similar for the files `systemmodel3V.m`--`systemmodel5V.m` and `main_3V.m`--`main_5V.m`.

To run the controller on MATLAB, you need to first install the needed packages, then for the 2-vehicle case, you just need to run the files `main_2V.m`--`main_5V.m` to see the results for different cases.




