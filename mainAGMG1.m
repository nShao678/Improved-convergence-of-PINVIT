% warning off
rng(1)
kMax = 10;
cphi2 = zeros(1,kMax);
boundkappa = zeros(1,kMax);
for k = 6:kMax
    n = 2^k-2;
    e = ones(n,1);
    K = spdiags([-e,2*e,-e],-1:1,n,n);
    I = speye(n);
    K = kron(K,I)+kron(I,K);
    
    agmg(K,[],1,[],[],[],[],1);
    B = @(x) agmg([],x,[],[],[],[],[],3); % set up for Bx=b, where B ~ A^{-1}
    [u1,~] = eigs(K,1,'smallestabs');
    x1 = B(u1); % B\u1
    [x2,~] = pcg(@(x) B(x),u1,[],[],@(x) K*x); % B*u1
    cphi2(k) = 1-(1/(u1'*x1)/(u1'*x2))^2;

    n = size(K,1);
    nu_max = eigs(@(x) K*(B(K*x)),size(K,1),K,1,'largestreal','IsSymmetricDefinite',1);
    nu_min = eigs(@(x) K*(B(K*x)),size(K,1),K,1,'smallestreal','IsSymmetricDefinite',1);
    kappa_nu = nu_max/nu_min;
    boundkappa(k) = 1-1/kappa_nu;

end
latex(sym(cphi2(6:10)))
latex(sym(boundkappa(6:10)))
latex(sym(cphi2(6:10)./boundkappa(6:10)))
save('dataAGMG1','cphi2','boundkappa')