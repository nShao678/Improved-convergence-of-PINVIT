rng(1)
warning('off');
addpath(genpath(pwd))
NMax = 6;
cphi2 = zeros(1,NMax);
boundkappa = cphi2;
for N = NMax
    tic
    n = 8-N;
    [K,M,B] = genDD(n,N);
    [u1,lambda] = eigs(K,M,1,'smallestabs');
    x1 = B(M*u1); % B\u1
    [x2,~] = pcg(@(x) B(x),u1,[],[],@(x) K*x); % B*u1
    cphi2(N) = 1-1/((u1'*M*x1)*(u1'*x2));
    nu_max = eigs(@(x) K*(B(K*x)),size(K,1),K,1,'largestreal','IsSymmetricDefinite',1,'IsFunctionSymmetric',1,'Tolerance',1e-3);
    nu_min = eigs(@(x) K*(B(K*x)),size(K,1),K,1,'smallestreal','IsSymmetricDefinite',1,'IsFunctionSymmetric',1,'Tolerance',1e-3);
    kappa_nu = nu_max/nu_min;
    boundkappa(N) = 1-1/kappa_nu;
    toc
end
latex(sym(cphi2(2:6)))
latex(sym(boundkappa(2:6)))
latex(sym(cphi2(2:6)./boundkappa(2:6)))
save('data1','cphi2','boundkappa')



nMax = 6;
cphi2 = zeros(1,nMax);
boundkappa = cphi2;
bd = cphi2;
err = zeros(2,nMax);
type = 0;
for n = 2:nMax
    tic
    N = 2;
%     n = 8-N;
    [K,M,B] = genDD(n,N);
    [u1,lambda] = eigs(K,M,1,'smallestabs');
    x1 = B(M*u1); % B\u1
    [x2,~] = pcg(@(x) B(x),u1,[],[],@(x) K*x); % B*u1
    cphi2(n) = 1-1/((u1'*M*x1)*(u1'*x2));

    nu_max = eigs(@(x) K*(B(K*x)),size(K,1),K,1,'largestreal','IsSymmetricDefinite',1,'IsFunctionSymmetric',1);
    nu_min = eigs(@(x) K*(B(K*x)),size(K,1),K,1,'smallestreal','IsSymmetricDefinite',1,'IsFunctionSymmetric',1);
    kappa_nu = nu_max/nu_min;
    boundkappa(n) = 1-1/kappa_nu;
    toc
end
latex(sym(cphi2(2:6)))
latex(sym(boundkappa(2:6)))
latex(sym(cphi2(2:6)./boundkappa(2:6)))
save('data2','cphi2','boundkappa')
