function  sol = solGet(probSDP, solSDP, nSol)

% Get the solution from the moment matrix calculated from the dual problem
droptol = 1e-3;
choltol = 1e-6;

nd = probSDP.K.s(1)^2;
matD = momGen(probSDP.A(:,1:nd), solSDP.dual);
matD = (matD + matD')/2;
nvar = probSDP.info.nvar;
deg = probSDP.info.deg;
order = probSDP.info.order;

[Phi, Lambda] = eig(matD);
[Lambda, eigInd] = sort((diag(Lambda)),'descend') ;
Phi = Phi(:,eigInd);
nlam = length(Lambda);

% Set rank of M
if(nargin > 2)
    rankM = nSol;
else
    dropInd = find(Lambda(2:nlam)./Lambda(1:nlam-1) < droptol);
    if ~isempty(dropInd)
        rankM = dropInd(1);
    else
        rankM = nlam;
    end
end

flag = 0;
while ~flag
    ind = 1:rankM;
    
    V = Phi(:,ind)*diag(sqrt(Lambda(ind)));
    
    [U, monb] = rcef(V, choltol);
    
    t = ceil(deg/2) + order;
    monIndDense = multiMatGen(nvar, t, monb);
    
    if max(max(monIndDense)) <= size(U,1)
        flag = 1;
        break;
    else
        rankM = rankM - 1;
    end
    
end

N = cell(nvar, 1);
for i = 1:nvar
    [~, monInd] = intersect(probSDP.info.mon, monIndDense(i, :));
    N{i} = U(monInd,:);
end

% compute common roots with the algorithm described in
% [Corless, Gianni, Trager. A reordered Schur factorization method
% for zero-dimensional polynomial systems with multiple roots.
% Proc. ISSAC (Maui), W. Kuechlin (Ed) pp. 133-140, 1997]
% random combinations of multiplication matrices
coef = rand(nvar, 1);
coef = coef / sum(coef);
MM = zeros(rankM);
for i = 1:nvar,
    MM = MM + coef(i)*N{i};
end

% ordered Schur decomposition of
% random combinations of multiplication matrices
[Q, T] = orderschur(MM);

% retrieve optimal vectors
% it is assumed than there is no multiple root
sol = cell(rankM, 1);

for i = 1:rankM
    sol{i} = zeros(nvar,1);
    for j = 1:nvar
        sol{i}(j) = Q(:,i)'*N{j}*Q(:,i);
    end;
end

end