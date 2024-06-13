function [y,yf,yH] = ddm(r,K,KH,vH,idxI,idx0)
x = r;
r = zeros(size(K,1),1);
yf = zeros(size(K,1),1);

r(idx0) = x;
rH = vH'*r;
uH = KH\rH;
yH = vH*uH;

yy = cell(length(idxI),1);
for iter = 1:length(idxI)
    yy{iter} = zeros(length(yH),1);
    if ~isempty(idxI{iter})
        idx = idxI{iter};
        yy{iter}(idx) = K(idx,idx)\r(idx);
    end
end
for iter = 1:length(idxI)
    yf = yf+yy{iter};
end
yf = yf(idx0);
yH = yH(idx0);
y = yf+yH;
end

