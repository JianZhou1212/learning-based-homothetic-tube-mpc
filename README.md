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
The article includes three case studies. For the convenience of the reviewers during the peer-review process, we provide the code for the first case study, which involves adjusting the homotopy iterations $n$ in the optimization problem. The full implementation will be published if the paper is accepted.

The file `MinTimeCoGDoubleTrackInitialization.m` defines the optimization problem for initializing the homotopic optimization.

The file `MinTimeCoGDoubleTrack.m` defines the optimization problem during the homotopy iterations.

The file `parameters_define.m` defines the static parameters.

The files `main_1.m` to `main_11.m` contain the cases for implementing the homotopic optimization with different iterations $n$. For example, in `main_1.m`, $n=120$. The only difference in these files is just different values of $n$.

To implement the method, first run `parameters_define.m`, then you can choose any case from `main_1.m` to `main_11.m` to get the results with different $n$ values.




