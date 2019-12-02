
function W_init = preprocess2(Data0,Mask,nZ,W_true,label)
addpath(genpath(pwd))
nS = length(Data0); % number of subjects
m = size(Data0{1},1); % number of variables
nT = size(Data0{1},2); % number of data points per subject or group

cv = [];
for i = 1:nS
    c{i} = corr(Data0{i}');
    cv = [cv, matrix2vec(c{i},Mask)];
end
idx = kmeans(cv',3); % three groups

Data1 = [];
Data2 = [];
Data3 = [];
for i = 1:nS
    if(idx(i)==1)
        Data1 = [Data1 Data0{i}];
    end
    if(idx(i)==2)
        Data2 = [Data2 Data0{i}];
    end
    if(idx(i)==3)
        Data3 = [Data3 Data0{i}];
    end
end

C1 = corr(Data1');
C2 = corr(Data2');
C3 = corr(Data3');
Mask1 = ones(m,m) - eye(m);
Mask2 = ones(m,m) - eye(m);
Mask3 = ones(m,m) - eye(m);
% Mask1 = (abs(C1)>=0.1);
% Mask1 = Mask1-eye(m);
% Mask2 = (abs(C2)>=0.1);
% Mask2 = Mask2-eye(m);
Lasso1 = zeros(m,m);
Lasso2 = zeros(m,m);
Lasso3 = zeros(m,m);
for i = 1:m
    id1 = find(Mask1(i,:)~=0);
    a1 = lasso(Data1(id1,:)',Data1(i,:)');
    Lasso1(i,id1) = a1(:,50)';
    
    id2 = find(Mask2(i,:)~=0);
    a2 = lasso(Data2(id2,:)',Data2(i,:)');
    Lasso2(i,id2) = a2(:,50)';
    
    id3 = find(Mask3(i,:)~=0);
    a3 = lasso(Data3(id3,:)',Data3(i,:)');
    Lasso3(i,id3) = a3(:,50)';
end
W_init{1} = eye(m)-Lasso1;
W_init{2} = eye(m)-Lasso2;
W_init{3} = eye(m)-Lasso3;





