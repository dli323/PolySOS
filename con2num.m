function [c, c0, pw] = con2num(conpol, X)
Xc = sort(symvar(conpol));
nc = length(Xc);
degc = double(feval(symengine, 'degree', conpol));
tc = degc/2;

item = feval(symengine, 'poly2list', conpol, X);
nmon = length(item);

c0 = 0;

for i = 1:nmon
    item_i = item(i); 
    if ~all(double(item_i(2)) == 0)
        c(i) = double(item_i(1));
        pw(i,:) = double(item_i(2)); 
    else
        c0 = double(item_i(1));
    end
end

end