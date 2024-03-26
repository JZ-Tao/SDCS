
function [G1, G2] = Ltrans(X, m, n)
  G1 = X(1:m-1,:,:)-X(2:m,:,:);
  G2 = X(:,1:n-1,:)-X(:,2:n,:);
end