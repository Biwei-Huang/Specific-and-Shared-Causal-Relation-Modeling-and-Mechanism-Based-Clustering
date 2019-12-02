% mechanism-based clustering
function Pz = clustering_ins(Data,thetaW,thetaE,Mask,nZ)
% INPUT: 
%     Data: data from each subject are saved in a cell
%     thetaW: related parameters of W
%     thetaE: related parameters of E
%     Mask: use it to fix some entries of B to zero, where B = I-W
%     nZ: number of groups

% OUTPUT:
%     Pz: Pz(i,j) the posterior probability that individual i is in group j


nS = length(Data); % number of subjects 
m = size(Data{1},1); % number of variables

%%
% do mechanism-based clustering
Pz = zeros(nS,nZ); % P(Zk = 1| X^s)
for s = 1:nS % the sth subejct
    X = Data{s};
    Ts = size(X,2);
    for z = 1:nZ % the zth group
        % samlping 
        W_vec_random = mvnrnd(thetaW{z}.mu,thetaW{z}.Sigma,80000);
        for t = 1:size(W_vec_random,1) % the tth generated sample
            W = vec2matrix(W_vec_random(t,:),Mask) + eye(m);
            f = Ts*log(abs(det(W)));
            Qen = zeros(Ts,1);
            E = W*X;
            for q2 = 1:length(thetaE{z})
                Qen = Qen + thetaW{z}.pi *thetaE{z}{q2}.pi * mvnpdf(E',thetaE{z}{q2}.mu,thetaE{z}{q2}.Sigma);
            end
            f = f + sum(log(Qen));
            Pz(s,z) = Pz(s,z) + exp(f);
        end
        Pz(s,z) =  Pz(s,z)/size(W_vec_random,1);
    end
    Pz(s,:) = Pz(s,:)/sum(Pz(s,:));
end



