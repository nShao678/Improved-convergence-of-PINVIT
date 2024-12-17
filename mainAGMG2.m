% warning off
rng(1)
kMax = 10;
cphi2 = zeros(1,kMax);
boundkappa = zeros(1,kMax);
load('dataAGMG1.mat')
iterMax = 1000;
kMax = 10;
EP = zeros(2,kMax);
for k = 6:kMax
    tic
    n = 2^k-2;
    e = ones(n,1);
    K = spdiags([-e,2*e,-e],-1:1,n,n);
    I = speye(n);
    K = kron(K,I)+kron(I,K);
    
    agmg(K,[],1,[],[],[],[],1);
    B = @(x) agmg([],x,[],[],[],[],[],3); % set up for Bx=b, where B ~ A^{-1}
    [U,Lambda] = eigs(K,2,'smallestabs');
    lambda2 = Lambda(2,2);
    len = size(K,1);
    u1 = U(:,1);
%     x1 = B(u1); % B\u1
    [x2,~] = pcg(@(x) B(x),u1,[],[],@(x) K*x); % B*u1
    u1x2 = u1'*x2;

    for iter = 1:iterMax

    w = randn(len,1);
    ww = B(w);
    distB = (w'*u1)^2/(w'*ww)/(u1x2);
    if distB>cphi2(k)
        EP(1,k) = EP(1,k)+1;
    end
    rho = ww'*K*ww/(ww'*ww);
    if rho<lambda2
        EP(2,k) = EP(2,k)+1;
    end
    end
    toc

end
latex(sym(EP(:,6:10)))
save('dataEPAGMG','EP')