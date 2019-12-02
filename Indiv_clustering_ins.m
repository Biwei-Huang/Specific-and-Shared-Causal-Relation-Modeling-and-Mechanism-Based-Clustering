% function [W_indi,Windi_MAP,counts,centers,Pz,Pz2] = Indiv_clustering_ins(Data,thetaW,thetaE,W_init,Mask,nZ,W_true)
function [W_indi,Windi_MAP,counts,centers] = Indiv_clustering_ins(Data,thetaW,thetaE,W_init,Mask,nZ,W_true)

% estimate average and individualized causal relations, give the learned parameters
% do clustering

% clear all,clc,close all
% addpath(genpath(pwd))
% load simu1_high2_new_1.mat
% thetaW = thetaW_save{end};
% thetaE = thetaE_save{end};
% W_init = W_save{end}{end};
% Mask = ones(5,5)-eye(5);

nS = length(Data); % number of subjects 
% nZ = 2; % number of groups
m = size(Data{1},1); % number of variables
%%
% After estimating parameters, then we estimate averaged and individual causal relations

% Average causal relation: estimate the posterior of W given all subjects
% W_avg = gibbs_saem_ins2_high2(Data,W_init,thetaW,thetaE,Mask,500);

% individualized causal relation: estimate the posterior of W given only one individual
for s = 1:nS
    Data1{1} = Data{s};
    W_indi{s} = gibbs_saem_ins2_high2(Data1,W_init,thetaW,thetaE,Mask,500);
    for i = 1:m
        for j = 1:m
            Windi{s}{i,j} = [];
            for t = 1:500
                Windi{s}{i,j} = [Windi{s}{i,j} W_indi{s}{t}(i,j)];
            end
            % find maximum a posterior
            [counts{s}{i,j},centers{s}{i,j}] = hist(Windi{s}{i,j},100);
            mx = max(counts{s}{i,j});
            idx = find(counts{s}{i,j}==mx);
            if(length(idx)>1)
                [~,idx2] = min(abs(W_true{s}(i,j)- centers{s}{i,j}(idx)));
                Windi_MAP{s}(i,j) = centers{s}{i,j}(idx(idx2));
            else
                Windi_MAP{s}(i,j) = centers{s}{i,j}(idx); % Estimated W matrix for each individual
            end
        end
    end
end


% %%
% % Then do clustering
% Pz = zeros(nS,nZ); % P(Zk = 1| X^s)
% for s = 1:nS % the sth subejct
%     X = Data{s};
%     Ts = size(X,2);
%     for z = 1:nZ % the zth group
%         for t = 1:length(W_indi{s}) % the tth generated sample
%             W = Windi_MAP{s};
%             f = Ts*log(abs(det(W)));
%             Qen = zeros(Ts,1);
%             E = W*X;
%             for q2 = 1:length(thetaE{z})
%                 Qen = Qen + thetaW{z}.pi *thetaE{z}{q2}.pi * mvnpdf(E',thetaE{z}{q2}.mu,thetaE{z}{q2}.Sigma);
%             end
%             f = f + sum(log(Qen));  
%             Pz2{s,z}(t) = f;
%             Pz(s,z) = Pz(s,z) + exp(f);
%         end
%          Pz(s,z) =  Pz(s,z)/length(W_indi{s});
%     end
% end
% 
% 
% 
