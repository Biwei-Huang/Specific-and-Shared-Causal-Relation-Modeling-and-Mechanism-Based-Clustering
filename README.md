# Specific-and-Shared-Causal-Relation-Modeling-and-Mechanism-Based-Clustering

A unified framework for causal discovery and mechanism-based group identification.

### EXAMPLE: 
see example.m

### IMPORTANT FUNCTIONS:
#### Main function for causal discovery: 
function [W_hat, thetaW_hat,thetaE_hat, W_save,thetaW_save,thetaE_save,Q_save] = SAEM_ins(Data,thetaW0,thetaE0,np,Mask)
* INPUT: 
  *  Data: data from each subject are saved in a cell
  *  thetaW0: initial values of the parameters related to W (W = I-B and B is the causal adjacency matrix), see equation (7) in the paper
  *  thetaE0: initial values of the parameters related to E
  *  np: number of particles that need to be sampled
  *  Mask: use it to fix some entries of B to zero, where B = I-W

* OUTPUT:
  *   W_hat: estimated W for each individual
  *   thetaW_hat: estimated parameters related to W
  *   thetaE_hat: estimated parameters related to E
  *   W_save: sampled W's in each iteration
  *   thetaW_save: estimated theta_W in each iteration
  *   thetaE_save: estimated theta_E in each iteration
  *   Q_save: estimated Q value in each iteration


#### Main function for mechanism-based clustering:
function Pz = clustering_ins(Data,thetaW,thetaE,Mask,nZ)
* INPUT: 
  *   Data: data from each subject are saved in a cell
  *   thetaW: related parameters of W
  *   thetaE: related parameters of E
  *   Mask: use it to fix some entries of B to zero, where B = I-W
  *   nZ: number of groups

* OUTPUT:
  *   Pz: Pz(i,j) the posterior probability that individual i is in group j



### Citation:
B. Huang, K. Zhang, P. Xie, M. Gong, E. Xing, C. Glymour. Specific and Shared Causal Relation Modeling and Mechanism-based Clustering. NeurIPS'19

If there are any questons, please send emails to biweih@andrew.cmu.edu
