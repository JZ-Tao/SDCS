function X = Lforward(P1, P2, m, n)
T = size(P1,3);
X = zeros(m,n,T);
X(1:m-1,:,:) = P1;
X(:,1:n-1,:) = X(:,1:n-1,:)+P2;
X(2:m,:,:) = X(2:m,:,:)-P1;
X(:,2:n,:) = X(:,2:n,:)-P2;


