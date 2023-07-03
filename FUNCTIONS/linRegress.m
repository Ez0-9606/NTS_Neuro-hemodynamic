function [b, rg, coefRg, coefSh] = linRegress(x,y)
    b = regress(y, x);
    
    rg = zeros(size(x,1),1);
    for i = 1:numel(b)
        rg = rg + b(i)*x(:,i);
    end
    shuf = rg(randperm(numel(y)));

    corRg = corrcoef(rg, y);
    corSh = corrcoef(shuf, y);

    coefRg = corRg(1,2);
    coefSh = corSh(1,2);
end