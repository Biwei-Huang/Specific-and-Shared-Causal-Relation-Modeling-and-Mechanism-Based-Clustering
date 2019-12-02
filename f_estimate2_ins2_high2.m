function f = f_estimate2_ins2_high2(Data,W,thetaW,thetaE,Mask)
% estimate the joint probability give Data and W
nS = length(Data); % number of subjects or groups
m = size(Data{1},1); % number of variables

f = 0;
for s = 1:nS
    Wtmp = squeeze(W(:,:,s));
    X = Data{s};
    T0 = size(X,2);
    E = Wtmp*X;
    Wtmp2 = reshape(Wtmp,m*m,1);
    idx = find(Mask==1);
    W_tilde = Wtmp2(idx);
    
    f = T0*log(abs(det(Wtmp)));
    
    Qwn = 0;
    for q = 1:length(thetaW)
        Qwn = Qwn + thetaW{q}.pi * mvnpdf(W_tilde',thetaW{q}.mu,thetaW{q}.Sigma);
    end
    if(Qwn==0)
        Qwn = 1e-300;
    end
    f = f + log(Qwn);
    
    Qen = zeros(T0,1);
    for q = 1:length(thetaE)
        for q2 = 1:length(thetaE{q})
            %             Qen = Qen + thetaE{q}.pi *thetaE{q}{q2}.pi * mvnpdf(E',thetaE{q}{q2}.mu,thetaE{q}{q2}.Sigma);
            Qen = Qen + thetaW{q}.pi *thetaE{q}{q2}.pi * mvnpdf(E',thetaE{q}{q2}.mu,thetaE{q}{q2}.Sigma);
            %             Note: thetaE{q}.pi and thetaW{q}.pi are equal
        end
    end
    if(~isempty(find(Qen==0)))
        idx = find(Qen==0);
        Qen(idx) = 1e-300;
    end
    f = f + sum(log(Qen));
    
    
end

