function [thetaW_old,thetaE_old,Q_old] = EM_inner_new2_high2(Data,thetaW_old,thetaE_old,gamma,W_save,Mask)
% inner EM to update parameters in M step
% model the whole W matrix and E with mixture of Gaussian in high dimensions
nS = length(Data); % number of subjects or groups
m = size(Data{1},1); % number of variables
thre = 1e-4;
% initialize the parameters and evaluate the initial value of log likelihood
Q_old = 0;
for i = length(W_save):-1:1
    gp = gamma(i);
    for ii = i+1:length(W_save)
        gp = gp*(1-gamma(ii));
    end
    Qtmp = 0;
    for j = 1:length(W_save{i})
        f = f_estimate2_ins2_high2(Data,W_save{i}{j},thetaW_old,thetaE_old,Mask);
        Qtmp = Qtmp + f;
    end
    Q_old = Q_old + gp*Qtmp/length(W_save{i});
end


% estimate the parameters with EM
count = 0;
while(1)
    count = count+1;
    if(count>15)
        break;
    end
    % E step
    alphaW = [];
    alphaE = [];
    alphaEE = [];
    for i = 1:length(W_save)
        for j = 1:length(W_save{i})
            for s = 1:nS
                X = Data{s};
                Ts = size(X,2);
                W = squeeze(W_save{i}{j}(:,:,s));
                E = W*X;
                idx = find(Mask==1);
                W_tilde = reshape(W,size(X,1)*size(X,1),1);
                W_tilde = W_tilde(idx);
                
                tmp = 0;
                for q = 1:length(thetaW_old)
                    tmp = tmp + thetaW_old{q}.pi * mvnpdf(W_tilde',thetaW_old{q}.mu,thetaW_old{q}.Sigma);
                end
                for q = 1:length(thetaW_old)
                    alphaW(i,j,s,q) = thetaW_old{q}.pi * mvnpdf(W_tilde',thetaW_old{q}.mu,thetaW_old{q}.Sigma)/tmp;
                    if(isnan(alphaW(i,j,s,q)) | alphaW(i,j,s,q)==0)
                        alphaW(i,j,s,q) = 1e-300; % to avoid numerical problem
                    end
                end
                alphaW(i,j,s,:) = alphaW(i,j,s,:)/sum(alphaW(i,j,s,:)); % to make sure the sum is 1
                
                tmp = zeros(Ts,1);
                for q = 1:length(thetaE_old)
                    for q2 = 1:length(thetaE_old{q})
                        tmp = tmp + thetaW_old{q}.pi *thetaE_old{q}{q2}.pi * mvnpdf(E',thetaE_old{q}{q2}.mu,thetaE_old{q}{q2}.Sigma);
                    end
                end
                for q = 1:length(thetaE_old)
                    tmp2 = zeros(Ts,1);
                    for q2 = 1:length(thetaE_old{q})
                        alphaEE(i,j,s,:,q,q2) = thetaW_old{q}.pi *thetaE_old{q}{q2}.pi * mvnpdf(E',thetaE_old{q}{q2}.mu,thetaE_old{q}{q2}.Sigma)./tmp;
                        % i,j,s,t,q,q2
                        tmp2 = tmp2 + thetaW_old{q}.pi *thetaE_old{q}{q2}.pi * mvnpdf(E',thetaE_old{q}{q2}.mu,thetaE_old{q}{q2}.Sigma);                        
                    end
                    alphaE(i,j,s,:,q) = tmp2./tmp;
                    if(~isempty(find(isnan(alphaE(i,j,s,:,q))==1)))
                        idx = find(isnan(alphaE(i,j,s,:,q))==1);
                        alphaE(i,j,s,idx,q) = 1e-300;
                    end
                    if(~isempty(find(alphaE(i,j,s,:,q)==0)))
                        idx = find(alphaE(i,j,s,:,q)==0);
                        alphaE(i,j,s,idx,q) = 1e-300;
                    end
                    % i,j,s,t,q
                end                
            end
        end
    end
    
    
    % M step: update the parameters
    
    % update thetaE
    for q = 1:length(thetaE_old)
        for q2 = 1:length(thetaE_old{q})
            for s = 1:nS
                X = Data{s};
                Ts = size(X,2);
                tmpnq_mu1 = zeros(Ts,m);
                tmpnq_mu2 = zeros(Ts,1);
                tmpnq_pi1 = zeros(Ts,1);
                tmpnq_pi2 = zeros(Ts,1);
                for i = length(W_save):-1:1
                    gp = gamma(i);
                    for ii = i+1:length(W_save)
                        gp = gp*(1-gamma(ii));
                    end
                    for j = 1:length(W_save{i})
                        W = squeeze(W_save{i}{j}(:,:,s));
                        E = W*X;
                        tmpnq_mu1 = tmpnq_mu1 + gp*repmat(squeeze(alphaEE(i,j,s,:,q,q2)),1,m).*E';
                        tmpnq_mu2 = tmpnq_mu2 + gp*squeeze(alphaEE(i,j,s,:,q,q2));
                        tmpnq_pi1 = tmpnq_pi1 + gp*squeeze(alphaEE(i,j,s,:,q,q2));
                        tmpnq_pi2 = tmpnq_pi2 + gp*squeeze(alphaE(i,j,s,:,q));
                    end
                end
                tmpnq_mu1s(s,:) = sum(tmpnq_mu1);
                tmpnq_mu2s(s) = sum(tmpnq_mu2);
                tmpnq_pi1s(s) = sum(tmpnq_pi1);
                tmpnq_pi2s(s) = sum(tmpnq_pi2);
            end
            if(sum(tmpnq_mu2s)==0)
                thetaE_new{q}{q2}.mu = zeros(1,m);
            else
                thetaE_new{q}{q2}.mu = sum(tmpnq_mu1s)/sum(tmpnq_mu2s);
            end
            thetaE_new{q}{q2}.pi = sum(tmpnq_pi1s)/sum(tmpnq_pi2s);
            
            for s = 1:nS
                X = Data{s};
                Ts = size(X,2);
                tmpnq_sg1 = zeros(m,m,Ts);
                tmpnq_sg2 = zeros(Ts,1);
                for i = length(W_save):-1:1
                    gp = gamma(i);
                    for ii = i+1:length(W_save)
                        gp = gp*(1-gamma(ii));
                    end
                    for j = 1:length(W_save{i})
                        W = squeeze(W_save{i}{j}(:,:,s));
                        E = W*X;
                        tmpnq_sg1 = tmpnq_sg1 + gp*mproduct3(squeeze(alphaEE(i,j,s,:,q,q2)),(E-repmat(thetaE_new{q}{q2}.mu',1,Ts)));
                        tmpnq_sg2 = tmpnq_sg2 + gp*squeeze(alphaEE(i,j,s,:,q,q2));
                    end
                end
                tmpnq_sg1s(:,:,s) = sum(tmpnq_sg1,3);
                tmpnq_sg2s(s) = sum(tmpnq_sg2);
            end
%             thetaE_new{q}{q2}.Sigma = diag(diag(sum(tmpnq_sg1s,3)))/sum(tmpnq_sg2s);
            thetaE_new{q}{q2}.Sigma = diag(diag(sum(tmpnq_sg1s,3)))/sum(tmpnq_sg2s) + 1e-5*eye(m); % to guarantee positivess
            tmp = thetaE_new{q}{q2}.Sigma;
            if(sum(sum(isnan(tmp)))>1)
               fprintf('warning NaN');
               save error1 tmpnq_sg1s tmpnq_sg2s tmpnq_sg1 tmpnq_sg2 alphaEE thetaE_new
               save error1_whole
            end
            if(sum(sum(isinf(tmp)))>1)
               fprintf('warning inf');
               save error1 tmpnq_sg1s tmpnq_sg2s tmpnq_sg1 tmpnq_sg2 alphaEE thetaE_new
               save error1_whole
            end
        end
    end
    
    
    % update thetaW
    for q = 1:length(thetaW_old)
        tmpnq_mu1 = zeros(1,length(W_tilde));
        tmpnq_mu2 = 0;
        tmpnq_sg1 = zeros(length(W_tilde),length(W_tilde));
        tmpnq_sg2 = 0;
        tmpnq_pi1 = 0;
        tmpnq_pi2 = 0;
        for i = length(W_save):-1:1
            gp = gamma(i);
            for ii = i+1:length(W_save)
                gp = gp*(1-gamma(ii));
            end
            for j = 1:length(W_save{i})
                for s = 1:nS
                    X = Data{s};
                    Ts = size(X,2);
                    W = squeeze(W_save{i}{j}(:,:,s));
                    idx = find(Mask==1);
                    W_tilde = reshape(W,size(X,1)*size(X,1),1);
                    W_tilde = W_tilde(idx);
                    
                    tmpnq_mu1 = tmpnq_mu1 + gp*alphaW(i,j,s,q)*W_tilde';
                    tmpnq_mu2 = tmpnq_mu2 + gp*alphaW(i,j,s,q);
                    tmpnq_pi1 = tmpnq_pi1 + gp*(sum(sum(alphaEE(i,j,s,:,q,:))) + alphaW(i,j,s,q));
                    tmpnq_pi2 = tmpnq_pi2 + gp*(Ts+1);
                end
            end
        end
        if(tmpnq_mu2==0)
            thetaW_new{q}.mu = zeros(1,sum(Mask(:)));
        else
            thetaW_new{q}.mu = tmpnq_mu1/tmpnq_mu2;
        end
        thetaW_new{q}.pi = tmpnq_pi1/tmpnq_pi2;
        
        for i = length(W_save):-1:1
            gp = gamma(i);
            for ii = i+1:length(W_save)
                gp = gp*(1-gamma(ii));
            end
            for j = 1:length(W_save{i})
                for s = 1:nS
                    X = Data{s};
                    W = squeeze(W_save{i}{j}(:,:,s));
                    idx = find(Mask==1);
                    W_tilde = reshape(W,size(X,1)*size(X,1),1);
                    W_tilde = W_tilde(idx);
                    tmpnq_sg1 = tmpnq_sg1 + gp*alphaW(i,j,s,q)*(W_tilde-thetaW_new{q}.mu')*(W_tilde-thetaW_new{q}.mu')';
                    tmpnq_sg2 = tmpnq_sg2 + gp*alphaW(i,j,s,q);
                end
            end
        end
%         thetaW_new{q}.Sigma = diag(diag(tmpnq_sg1/tmpnq_sg2));
        thetaW_new{q}.Sigma = diag(diag(tmpnq_sg1/tmpnq_sg2)) + 1e-5*eye(size(tmpnq_sg1,1)); % to guarantee positivess
        tmp = thetaW_new{q}.Sigma;
        if(sum(sum(isnan(tmp)))>1)
            fprintf('warning NaN');
            save error2 tmpnq_sg1 tmpnq_sg2 alphaW thetaW_new
            save error2_whole
        end
        if(sum(sum(isinf(tmp)))>1)
            fprintf('warning inf');
            save error2 tmpnq_sg1 tmpnq_sg2 alphaW thetaW_new
            save error2_whole
        end
    end
    
    
        
    
    % evaluate the log likelihoood
    Q_new = 0;
    for i = length(W_save):-1:1
        gp = gamma(i);
        for ii = i+1:length(W_save)
            gp = gp*(1-gamma(ii));
        end
        Qtmp = 0;
        for j = 1:length(W_save{i})
            f = f_estimate2_ins2_high2(Data,W_save{i}{j},thetaW_new,thetaE_new,Mask);
            Qtmp = Qtmp + f;
        end
        Q_new = Q_new + gp*Qtmp/length(W_save{i});
    end
    
    
    if(isnan(Q_new))
        'warning'
        fprintf('%d\n',count);
    end
    
    if((Q_new-Q_old)/abs(Q_old)<=thre)
        break;
    else
        Q_old = Q_new;
        thetaW_old = thetaW_new;
        thetaE_old = thetaE_new;
    end
    
end

