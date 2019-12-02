function [f] = f_estimate_ins2_high2(X,W,thetaW,thetaE,Mask,idw,S)
% estimate the function at point S
iw = idw(1);
jw = idw(2);
m = size(X,1);
T = size(X,2);
f = zeros(1,length(S));
for l = 1:length(S)
    Wtmp = W;
    Wtmp(iw,jw) = S(l); % estimate the function at Wtmp
    E = Wtmp*X;
    
    id = find(Mask==1);
    W_tilde = reshape(Wtmp,m*m,1);
    W_tilde = W_tilde(id);
    
    f(l) = f(l) + T*log(abs(det(Wtmp)));
    
    Qwn = 0;
    for q = 1:length(thetaW)
        Qwn = Qwn + thetaW{q}.pi * mvnpdf(W_tilde',thetaW{q}.mu,thetaW{q}.Sigma);
    end
    if(Qwn==0)
        Qwn = 1e-300;
    end
    f(l) = f(l) + log(Qwn);
    
    Qen = zeros(T,1);
    for q = 1:length(thetaE)
        for q2 = 1:length(thetaE{q})
            Qen = Qen + thetaW{q}.pi *thetaE{q}{q2}.pi * mvnpdf(E',thetaE{q}{q2}.mu,thetaE{q}{q2}.Sigma);
            %              Note: thetaE{q}.pi and thetaW{q}.pi are equal
        end
    end
    if(~isempty(find(Qen==0)))
        idx = find(Qen==0);
        Qen(idx) = 1e-300;
    end
    f(l) = f(l) + sum(log(Qen));
    
    if(isnan(f(l)))
        'warning'
    end
end





