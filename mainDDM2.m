% clear all
%close all
rng(1)
warning('off');
load('data1.mat')
addpath(genpath(pwd))
sympref('FloatingPointOutput',true);
iterMax = 1000;
NMax = 6;
EP = zeros(2,NMax);



for N = 2:NMax
    tic
    n = 8-N;
    [K,M,B] = genDD(n,N);
    [U,Lambda] = eigs(K,M,2,'smallestabs');
    lambda2 = Lambda(2,2);
    len = size(K,1);
    u1 = U(:,1);
%     x1 = B(M*u1); % B\u1
    [x2,~] = pcg(@(x) B(x),u1,[],[],@(x) K*x); % B*u1
    u1x2 = u1'*x2;
    for iter = 1:iterMax

    w = randn(len,1);
    ww = B(M*w);
    distB = (w'*M*u1)^2/(w'*M*ww)/(u1x2);
    if distB>cphi2(N)
        EP(1,N) = EP(1,N)+1;
    end
    rho = ww'*K*ww/(ww'*M*ww);
    if rho<lambda2
        EP(2,N) = EP(2,N)+1;
    end
    end
    toc

end
save('dataEP1','EP')


load('data2.mat')
iterMax = 1000;
nMax = 6;
EP = zeros(2,nMax);

type = 0;

N = 2;
for n = 2:nMax
    tic
%     N = 2;
    
    [K,M,B] = genDD(n,N);
    [U,Lambda] = eigs(K,M,2,'smallestabs');
    lambda2 = Lambda(2,2);
    len = size(K,1);
    u1 = U(:,1);
%     x1 = B(M*u1); % B\u1
    [x2,~] = pcg(@(x) B(x),u1,[],[],@(x) K*x); % B*u1
    u1x2 = u1'*x2;
    for iter = 1:iterMax

    w = randn(len,1);
    ww = B(M*w);
    distB = (w'*M*u1)^2/(w'*M*ww)/(u1x2);
    if distB>cphi2(n)
        EP(1,n) = EP(1,n)+1;
    end
    rho = ww'*K*ww/(ww'*M*ww);
    if rho<lambda2
        EP(2,n) = EP(2,n)+1;
    end
    end
    toc

end
save('dataEP2','EP')
