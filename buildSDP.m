function [prob] = buildSDP(obj, cstr, userOpts)
% Build SDP in sedumi format

%---------------------------------------------------------------------
% Error checking
%---------------------------------------------------------------------
narginchk(1,3);

%---------------------------------------------------------------------
% Default options
%---------------------------------------------------------------------
opts.FR = 0;        % Facial reduction
                    % 0: Don't apply facial reduction
                    % 1: Apply facial reduction
opts.order = 0;     % Relaxation order
                    % Comments: for some problems, facial reduction seems
                    % cannot work toghter with higher order relaxation.

%---------------------------------------------------------------------
% Process Inputs
%---------------------------------------------------------------------
if nargin < 2
    % unconstrained optimization problem
    cstr = sym([]);
elseif nargin > 2
    opts = setUserOpts(opts, userOpts);
end

if ~isa(obj, 'sym')
    obj = sym(obj);
end

if ~isa(cstr, 'sym')
    cstr = sym(cstr);
end

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

dmax = max(d);
t = ceil(dmax/2) + opts.order;

% Generate selection matrices corresponding to objective function
Nvec = nchoosek(n + 2*t, 2*t);

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

Matall = Aall;
dims = dimsA;

clear Aall;

% Generate selection matrices corresponding to constraints
symLIST = @(varargin)feval(symengine,'DOM_LIST',varargin{:});
varcell = num2cell(X);

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

% Generate coefficient vector
[coeff, const, pw] = fun2num(objAll, symLIST(varcell{:}));
b = vecGen(2*t, pw, coeff);
b = -sparse(b);

Mat = -Matall(2:end, :);
Matc = Matall(1, :);

K.s = dims;

%---------------------------------------------------------------------
% Facial Reduction
%---------------------------------------------------------------------
FRdone = 0;
if opts.FR
    pMosek = convert_sedumi2mosek(Mat, b, Matc, K);
    [pMosekRd, sieveInfo] = SieveSDP(pMosek);
    if ((sieveInfo.n_pre > sieveInfo.n_post) && (sieveInfo.m_pre > sieveInfo.m_post))
        FRdone = 1;
    end
end

if FRdone
    [prob.A, prob.b, prob.c, prob.K] = convert_mosek2sedumi(pMosekRd);
    prob.A = prob.A';
    prob.c = prob.c';
    prob.info.mon = [1; sieveInfo.undeleted + 1]; % add the constant term
else
    prob.A = Mat;
    prob.b = b;
    prob.c = Matc;
    prob.K = K;
    prob.info.mon = [1:Nvec]';
end

prob.info.nvar = n;                 % number of variables
prob.info.deg = dmax;               % degree
prob.info.order = opts.order;       % relaxation order
prob.info.const = const;            % the constant term

clear Matall;

end