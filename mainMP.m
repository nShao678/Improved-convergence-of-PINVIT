ng(1)
lenN = 4;
EPSet = zeros(4,lenN);
for iterMat = 1:2
for iterN = 1:lenN
n = 256*(2^iterN);

X = randn(n,n);
Y = randn(n,n);
K = zeros(n,n);
if iterMat==1
for ii = 1:n
    for jj = 1:n
        K(ii,jj) = exp(-norm(X(:,ii)-X(:,jj))/2);
    end
end
else
for ii = 1:n
    for jj = 1:n
        K(ii,jj) = (X(:,ii)'*X(:,jj)+1)^3;
    end
end
K = (X'*X+1).^3+(Y'*Y+1).^3+1i*((X'*Y+1).^3-(Y'*X+1).^3);
end

Ks = single(K);
Rs = chol(Ks);


[U,lambda] = eigs(K,2,'smallestabs');
u1 = U(:,1);
x1 = (Rs\(Rs'\u1)); % B\u1
x2 = (Rs'*(Rs*u1));
x1 = double(x1);
x2 = double(x2);

cphi2 = 1-(1/(u1'*x1)/(u1'*x2))^2;

m = 1000;
Rq = zeros(m,1);
Theta = zeros(m,1);
for ii = 1:m
    xx = randn(n,1);
    xx = xx/norm(xx);
    Rq(ii) = xx'*(K*xx);
    xb = Rs'*(Rs*xx);
    xb = double(xb);
    Theta(ii) = ((xx'*x2).^2)/(u1'*x2)./(xb'*xx);
end
EPSet(2*iterMat-1,iterN) = sum(Rq(1:ii)<lambda(2))/m;
EPSet(2*iterMat,iterN) = sum(Theta(1:ii)>cphi2)/m;

end
end

latex(sym(EPSet))