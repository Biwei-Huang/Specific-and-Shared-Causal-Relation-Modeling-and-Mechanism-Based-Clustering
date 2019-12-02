% estimate the initial value of W

function W_init = preprocess(Data0,Mask,nZ)

nS = length(Data0); % number of subjects
m = size(Data0{1},1); % number of variables
nT = size(Data0{1},2); % number of data points per subject or group

cv = [];
for i = 1:nS
    c{i} = corr(Data0{i}');
    cv = [cv, matrix2vec(c{i},Mask)];
end
idx = kmeans(cv',nZ); % two groups

for k = 1:nZ
    Data{k} = [];
    for i = 1:nS
        if(idx(i)==k)
            Data{k} = [Data{k} Data0{i}];
        end
    end
    Lasso{k} = zeros(m,m);
    for i = 1:m
        id = find(Mask(i,:)~=0);
        a = lasso(Data{k}(id,:)',Data{k}(i,:)');
        Lasso{k}(i,id) = a(:,50)';
    end
    W_init{k} = eye(m)-Lasso{k};
end

