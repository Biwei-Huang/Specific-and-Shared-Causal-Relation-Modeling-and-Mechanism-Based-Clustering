function Y = mproduct3(Z,X)

Y = zeros(size(X,1),size(X,1),size(X,2));
for t = 1:size(X,2)
   Y(:,:,t) = X(:,t)*X(:,t)'*Z(t);
end