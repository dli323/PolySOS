function [sol] = solveSDP(prob, solver)
% Solve SDP formulated from polynomial optimization problem

A = prob.A;
b = prob.b;
c = prob.c;
K = prob.K;

% Checks on specified solver type and method
solver = lower(solver);
if ~any(strcmpi(solver,{'sedumi','mosek','sdpa','sdpnal+','cdcs','afom'}))
    error('Unknown opts.solver. Please use "sedumi", "mosek", "sdpa", "sdpnal+" or "cdcs".')
end

switch solver
    case {'sedumi'}
        pars.fid = 1;
        pars.eps = eps;
        pars.bigeps = eps;
        pars.sdp = 1;
        pars.vplot = 0;
        [x,y,info] = sedumi(A, b, c, K, pars);
        sol.primal = x;
        sol.dual = y;
        sol.info = info;
        
    case {'mosek'}
        probMosek = convert_sedumi2mosek(A, b, c, K);
        res = mosekCall(probMosek);
        
        sol.primal = res.sol.itr.barx;
        sol.dual = res.sol.itr.y;
        sol.info = res.info;
        
    case {'sdpa'}
        SedumiToSDPA('SDP-Data.dat-s', A, b, c, K, '%16.24e');
        [mDIM,nBLOCK,bLOCKsTRUCT,c,F] = read_data('SDP-Data.dat-s');
        [objVal,y,x,z,info] = sdpam(mDIM,nBLOCK,bLOCKsTRUCT,c,F);
        delete('SDP-Data.dat-s');
        
        sol.primal = x;
        sol.dual = y;
        sol.info = info;
        
        
    case {'sdpnal+'}
        opts.maxiter = 20000;
        opts.tol = 1e-6;
        
        [blk,At,cc,bb] = read_sedumi(A, b, c, K);
        [objv,Xv,sv,yv,Z1v,Z2v,y2v,vv,info,runhist] = ...
            sdpnalplus(blk,At,cc,bb,[],[],[],[],[],opts);
        
        sol.primal = Xv;
        sol.dual = yv;
        sol.info = info;
        
    case {'cdcs'}
        opts.solver = 'primal';
        opts.maxIter = 20000;
        opts.relTol = 1e-6;
        [x,y,z,info] = cdcs(A',b,c,K,opts);
        
        sol.primal = x;
        sol.dual = y;
        sol.info = info;
        
    case {'afom'}
        opts.maxIter = 100000;
        opts.relTol = 1e-6;
        opts.rescale = 1;
        [x,y,z,info] = afom(A,b,c,K,opts);
        
        sol.primal = x;
        sol.dual = y;
        sol.info = info;
        
    otherwise
        
        error('Unknown solver.')
        
end


end