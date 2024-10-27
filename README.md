# Learning-based Rigid Tube Model Predictive Control
This repository contains code for the article
```
@article{gao2024learninghomo,
  title={Learning-based Homothetic Tube MPC},
  author={Yulong Gao, Shuhao Yan, Jian Zhou, Mark Cannon},
  year={2024},
  pages={}
} 
```
## Packages for running the code
To run the code you need to install:

**CasADi**: https://web.casadi.org/;

**MPT**: https://www.mpt3.org/Main/Installation; (The installation will automatically install Yalmip, which is also necessary for running the code.)

## Introduction to the files
See the introductions to the two folders below.

## CarFollowingSimulation
This folder contains a subfolder **Function**, which include the basic functions for running the algorithm, and several scripts. To run the method, you need to first run **systemmodel2V.m**, which define the parameters offline. After this, you can run **main_simulation.m**, which compares three method: the proposed learning-based homothetic tube MPC, the learning-based rigid tube MPC [1], and the traditional homothetic tube MPC [2], where these three MPC methods are defined in the scripts **HomotheticMPC.m**, **RigidMPC.m**, and **TradHomotheicMPC.m**, respectively.

## RobustSatisfaction
The files in this folder are similar to those in the above. The script **main_robustness.m** exams the statistical satisfaction of the method in stochastic simulations.

## Reference
[1] Gao, Yulong, Shuhao Yan, Jian Zhou, Mark Cannon, Alessandro Abate, and Karl Henrik Johansson. "Learning-based rigid tube model predictive control." In 6th Annual Learning for Dynamics & Control Conference, pp. 492-503. PMLR, 2024.

[2] Kouvaritakis, Basil, and Mark Cannon. "Model predictive control." Switzerland: Springer International Publishing 38 (2016).