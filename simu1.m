clear all,clc,close all
% In this simulation, we add prior on the whole W and E
% W: mixture of Gaussian
% E: mixture of mixture of Gaussian

load Simu_structure % G:graph structures
m = 5; % number of variables
nZ = 2; % number of groups
nS = 100; % number of subjects
nT = 20; % number of data points per subject or group
T = nS*nT;
np = 30; % number of particles that need to be sampled

for trial = 1:20
    rng(trial)
    
    Mask = ones(m,m)-eye(m);
    B0 = G{trial}';
    Mask0 = B0;
    W0 = eye(m)-G{trial}';
    
    % parameters for W
    W_mu{1} = 0*ones(1,sum(B0(:)))';
    W_mu{2} = 1*ones(1,sum(B0(:)))';
    tmp = unique(randi(sum(B0(:)),[1,2]));
    W_mu{1}(tmp) = 1;
    W_sigma{1} = diag(0.01*ones(1,sum(B0(:))));
    W_sigma{2} = diag(0.01*ones(1,sum(B0(:))));
    W_pi(1) = 0.5;
    W_pi(2) = 0.5;
    
    thetaW_true.mu = W_mu;
    thetaW_true.sigma = W_sigma;
    thetaW_true.pi = W_pi;
    
    %parameters for E
    for q = 1:nZ % for each group
        if(q==1)
            E_mu{q} = [-0.3*ones(1,m);0.3*ones(1,m)];
        end
        if(q==2)
            E_mu{q} = [-0.5*ones(1,m);0.5*ones(1,m)];
        end
        E_sigma{q} = 0.1*ones(1,m);
        E_pi{q} = ones(1,2)/2;
        gmE{q} = gmdistribution(E_mu{q},E_sigma{q},E_pi{q});
    end
    
    thetaE_true.mu = E_mu;
    thetaE_true.sigma = E_sigma;
    thetaE_true.pi = E_pi;
    
    s = 1;
    while(1)
        Z(s) = randi(nZ);
        W{s} = vec2matrix(mvnrnd(W_mu{Z(s)}',W_sigma{Z(s)}),Mask0);
        W{s} = W{s} + eye(m);
        E{s} = random(gmE{Z(s)},nT)';
        Data{s} = inv(W{s})*E{s};
        
        % check extreme values
        if(max(abs(Data{s}(:)))<30)
            s = s+1;
        end
        if(s>nS)
            break;
        end
    end
    
    % preprocessing
    W_init = preprocess(Data,Mask,nZ,W,Z);
    
    for q = 1:nZ
        thetaW0{q}.mu = matrix2vec(W_init{q},Mask)';
        thetaW0{q}.Sigma = 0.01*eye(sum(Mask(:)));
        thetaW0{q}.pi = 0.5;
        for q2 = 1:2
            thetaE0{q}{1}.mu = -0.5*ones(1,m);
            thetaE0{q}{2}.mu = 0.5*ones(1,m);
            thetaE0{q}{q2}.Sigma = 0.1*eye(m);
            thetaE0{q}{q2}.pi = 0.5;
        end
        %thetaE0(q) = thetaW0{q}.pi; Note: THEY ARE EQUAL!!
    end
    
    fileID = fopen(strcat('simu1_high2_new_',num2str(rep),'.txt'),'w');
    fprintf(fileID,'%d\n',rep);
    fclose(fileID);
    
    filename = strcat('simu1_high2_new_',num2str(trial));
    [W_hat, thetaW_hat,thetaE_hat, W_save,thetaW_save,thetaE_save,Q_save] = SAEM_ins_new2_high2(Data,thetaW0,thetaE0,np,Mask,filename,thetaE_true,thetaW_true, W, Z)
    
end


