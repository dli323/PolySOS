function [y,objv,sieve_info] = PolySOS(obj, cstr, order)

% Make the vector of objective functions and constraints as a column vector
if size(obj, 1) == 1
    obj = obj.';
end

if size(cstr, 1) == 1
    cstr = cstr.';
end

Nobj = size(obj,1);
Ncstr = size(cstr,1);

Xo = sort(symvar(obj));
Xc = sort(symvar(cstr));
X = union(Xo, Xc);
n = length(X);

Xoe = cell(Nobj, 1);
varset = cell(Nobj, 1);
for i = 1:Nobj
    Xoe{i} = sort(symvar(obj(i)));
    [~,~,varset{i}] = intersect(Xoe{i}, X);
    varset{i} = varset{i}' - 1;
end

d = zeros(Nobj + Ncstr, 1);

for i = 1:Nobj
    d(i) = double(feval(symengine, 'degree', obj(i)));
end

for i = 1:Ncstr
    d(Nobj+i) = double(feval(symengine, 'degree', cstr(i)));
end

% dmax = max(d);
dmax = order;
t = dmax/2;


% Generate selection matrix/matrices corresponding to objective function
Nvec = nchoosek(n + dmax, dmax);

if Nobj == 1
    Nmat = nchoosek(n + t, t);
    [row, col] = matGen(n, t);
    Aall = sparse(row,col,ones(length(row),1));
    dimsA = Nmat;
    objAll = obj;
else
    Aall = [];
    dimsA = [];
    objAll = sym(0);
    for i = 1:Nobj
        Nmat = nchoosek(length(varset{i}) + t, t);
        [row, col] = matGen(n, t, varset{i});
        Aall = [Aall, sparse(row,col,ones(length(row),1),Nvec,Nmat^2)];
        dimsA = [dimsA, Nmat];
        objAll = objAll + obj(i);
    end
end

Matall = Aall;
dims = dimsA;

clear Aall;

% Generate coefficient vector
symLIST = @(varargin)feval(symengine,'DOM_LIST',varargin{:});
varcell = num2cell(X);

[coeff, const, pw] = fun2num(objAll, symLIST(varcell{:}));
b = vecGen(dmax, pw, coeff);
b = -sparse(b);

% Generate selection matrix/matrices corresponding to constraints
for i = 1 : Ncstr    
    [cstrcf, cstrcs, cstrpw] = con2num(cstr(i), symLIST(varcell{:}));
    [crow, ccol, cval] = conGen(n, t, cstrpw, cstrcf, cstrcs);
    nB = max(ccol);
    Ball = sparse(crow,ccol,cval,Nvec,nB);
    dimsB = sqrt(nB);
    Matall = [Matall, Ball];
    dims = [dims, dimsB];
    clear Ball;
end

Mat = -Matall(2:end, :);
Matc = Matall(1, :);

clear Matall;

pars.fid = 1;
pars.eps = eps;
pars.bigeps = eps;
pars.sdp = 1;
pars.vplot = 0;

K.s = dims;
% [x,y,info] = sedumi(Mat, b, Matc, K, pars);

prob = convert_sedumi2mosek(Mat, b, Matc, K);

sieve_option.maxiter = intmax;
sieve_option.epsilon = eps;
sieve_option.cholEPS = 0;

[probr, sieve_info] = SieveSDP(prob, sieve_option);
solve_r = mosekCall(probr);

if isfield(solve_r,'sol')
    y = solve_r.sol.itr.y;
    objv = -solve_r.obj1;
    info = solve_r.info;
else
    % solve the original SDP
    solve_original = mosekCall(prob);
    y = solve_original.sol.itr.y;
    objv = -solve_original.obj1;
    info = solve_original.info;
end

end