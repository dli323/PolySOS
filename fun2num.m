function [c, c0, pw] = fun2num(pol, X)
Xf = sort(symvar(pol));
n = length(Xf);
degree = double(feval(symengine, 'degree', pol));
t = degree/2;

item = feval(symengine,'poly2list', pol, X);
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