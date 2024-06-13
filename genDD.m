function [K,M,T] = genDD(n,N)
coordinates=[0 0; 1 0; 0 1; 1 1];
elements=[1 2 3; 2 4 3];   
dirichlet=[1 2; 2 4; 4 3; 3 1];
for iter = 1:N
    [coordinates,elements,dirichlet]=refinement_uniform(coordinates,elements,dirichlet);
end  
coordinatesN = coordinates;

for iter = 1:n
    [coordinates,elements,dirichlet]=refinement_uniform(coordinates,elements,dirichlet);
end  
[K,areas] = stiffness_matrixP1_2D(elements,coordinates, ones(size(coordinates,1),1)); 
M = mass_matrixP1_2D(elements,areas,ones(size(coordinates,1),1));
NN = 2^N;
idxI = cell(NN^2,1);
for iter1 = 1:NN
    for iter2 = 1:NN
        idxI{iter1+(iter2-1)*NN} = find( ...
            coordinates(:,1)<(iter1+0.5)/NN&coordinates(:,2)<(iter2+0.5)/NN ...
            &coordinates(:,1)>(iter1-1.5)/NN&coordinates(:,2)>(iter2-1.5)/NN ...
            &coordinates(:,1)<1&coordinates(:,2)<1 ...
            &coordinates(:,1)>0&coordinates(:,2)>0 ...
            );
    end
end
idx = unique(cell2mat(idxI));
idxN = find(coordinatesN(:,1)==0|coordinatesN(:,1)==1|coordinatesN(:,2)==0|coordinatesN(:,2)==1);
idxN = setdiff(1:size(coordinatesN,1),idxN);
vH0 = eye(size(coordinatesN,1));
vH = zeros(size(K,1),length(idxN));
for iter = 1:length(idxN)
    F = scatteredInterpolant(coordinatesN(:,1),coordinatesN(:,2),vH0(:,idxN(iter)));
    vH(:,iter) = F(coordinates(:,1),coordinates(:,2));
end
KH = vH'*K*vH;
T = @(x) ddm(x,K,KH,vH,idxI,idx);
K = K(idx,idx);
M = M(idx,idx);






end
