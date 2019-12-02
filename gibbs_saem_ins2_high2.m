% gibbs sampling (for SAEM)
function W_all_new = gibbs_saem_ins2_high2(Data,W_all,thetaW,thetaE,Mask,np)
% np: number of particles that need to be sampled
maxIter = 10;
nS = length(Data); % number of groups
m = size(Data{1},1); % number of variables

for s = 1:nS
    X = Data{s};
    W = squeeze(W_all(:,:,s));
    for iter = 1:maxIter % burn-in period
        for i = 1:m
            for j = 1:m
                if(Mask(i,j))
                    W(i,j) = IA2RMS_ins2_high2(X,W,thetaW,thetaE,Mask,i,j,2,0);  % adpatative rejection metropolis sampling
                end
            end
        end
    end
    
    W_new = W;
    for iter = 1:np % sampling period
        for i = 1:m
            for j = 1:m
                if(Mask(i,j))
                    W_new(i,j) = IA2RMS_ins2_high2(X,W_new,thetaW,thetaE,Mask,i,j,2,0);  % adpatative rejection metropolis sampling
                end
            end
        end
        W_all_new{iter}(:,:,s) = W_new;
    end
end

