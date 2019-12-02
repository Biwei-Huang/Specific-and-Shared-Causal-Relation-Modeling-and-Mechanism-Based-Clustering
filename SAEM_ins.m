% estimate the W matrix for each individual, where W = I-B and B is the
%   causal adjacency matrix
% estimate the parameters related to W and E
% using stochastic approximated EM
function [W_hat, thetaW_hat,thetaE_hat, W_save,thetaW_save,thetaE_save,Q_save] = SAEM_ins(Data,thetaW0,thetaE0,np,Mask)
% INPUT: 
%     Data: data from each subject are saved in a cell
%     thetaW0: initial values of the parameters related to W (W = I-B and B is the causal 
%              adjacency matrix), see equation (7) in the paper
%     thetaE0: initial values of the parameters related to E
%     np: number of particles that need to be sampled
%     Mask: use it to fix some entries of B to zero, where B = I-W

% OUTPUT:
%     W_hat: estimated W for each individual
%     thetaW_hat: estimated parameters related to W
%     thetaE_hat: estimated parameters related to E
%     W_save: sampled W's in each iteration
%     thetaW_save: estimated theta_W in each iteration
%     thetaE_save: estimated theta_E in each iteration
%     Q_save: estimated Q value in each iteration

nS = length(Data); % number of subjects
m = size(Data{1},1); % number of variables
% np: number of sampled particles
maxIter = 25; % maximum number of iterations
kappa = 1;
gamma = zeros(1,maxIter);
gamma(1:3) = 1;
for i = 4:maxIter
    gamma(i) = 0.98;
end
% gamma(55:end) = gamma(54)*(((0:maxIter-55)+kappa)/kappa).^(-0.4);

% initialize
for s = 1:nS
    W_old(:,:,s) = zeros(m,m);
    for i = 1:m
        for j = 1:m
            if(i==j)
                W_old(i,j,s) = 1;
            else
                if(i>j)
                    W_old(i,j,s) = 0.5;
                end
            end
        end
    end
end
thetaW_old = thetaW0; thetaE_old = thetaE0;

% SAEM iteration
for iter = 1:maxIter
    iter
    % E step
    W_new = gibbs_saem_ins2_high2(Data,W_old,thetaW_old,thetaE_old,Mask,np);
    
    % M step (inner EM to update parameters)
    W_save{iter} = W_new;
    [thetaW_new,thetaE_new,Q_new] = EM_inner_new2_high2(Data,thetaW_old,thetaE_old,gamma,W_save,Mask);
    
    thetaW_save{iter} = thetaW_new;
    thetaE_save{iter} = thetaE_new;
    Q_save{iter} = Q_new;
    
    %     if(abs((Q_new-Q_old)/Q_old)<threshold)
    %         break;
    %     end
    thetaW_old = thetaW_new;
    thetaE_old = thetaE_new;
    Q_old = Q_new;
    W_old = W_new{end};
    
end

W_hat = W_new;
thetaW_hat = thetaW_new;
thetaE_hat = thetaE_new;



