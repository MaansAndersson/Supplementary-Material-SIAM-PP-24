function [RESVEC, INNER] = complexsolvetimer(solver_type,max_it,tol_res,C,c,xref,M0,M1,x0,inner_solver,max_inner_iter, repeats, log_file)
    A = real(C);
    B = imag(C);
    n = size(A);
    
    clear global inner_iter;
    
    timings = true;
    times = zeros(repeats, 1);
    
    for i=1:repeats
        %We want to keep track of the number of inner iterations used by the
        %different preconditioners. Since these are passed as a function handle
        %to funcitons like GMRES, we cannot control the return value, so use a
        %global variable instead.
        global inner_iter
        inner_iter = []; 
    
        tic;
        switch solver_type
            case {'GMRES', 'gmres'}
                [X,FLAG,RELRES,ITER,RESVEC] = gmres(C,c,[],tol_res,max_it,[],[],x0);
                ITER = ITER(end);
                INNER = NaN;
            case {'PMHSS-GMRES', 'phmss-gmres'}
                switch inner_solver
                    case {'\', 'backslash', 'BACKSLSASH'}
                        inner_iter = NaN;
                        [X,FLAG,RELRES,ITER,RESVEC] = gmres(C,c,[],tol_res,max_it,@pmhss_prec,[],[],A,B);
                    case {'CG', 'cg', 'PCG', 'pcg'}
                        [X,FLAG,RELRES,ITER,RESVEC] = gmres(C,c,[],tol_res,max_it,@pmhss_precCG__,[],[],A,B,max_inner_iter);
                    otherwise
                        inner_iter = NaN;
                        [X,FLAG,RELRES,ITER,RESVEC] = gmres(C,c,[],tol_res,max_it,@pmhss_prec,[],[],A,B);
                end
                ITER = length(RESVEC);
                INNER = 1;
            case {'PRESB-GMRES', 'presb-gmres'}
                Ablock = [A -B; B A];
                bblock = [real(c); imag(c)];
                switch inner_solver
                    case {'\', 'backslash', 'BACKSLASH'}
                        inner_iter = NaN;
                        [X,FLAG,RELRES,ITER,RESVEC] = gmres(Ablock,bblock,[],tol_res,max_it,@blkprec,[],[],A,B); %,x0);
                    case {'CG', 'cg', 'PCG', 'pcg'}
                        [X,FLAG,RELRES,ITER,RESVEC] = gmres(Ablock,bblock,[],tol_res,max_it,@blkprecCG__,[],[],A,B,max_inner_iter); %,x0);
                    otherwise
                        inner_iter = NaN;
                        [X,FLAG,RELRES,ITER,RESVEC] = gmres(Ablock,bblock,[],tol_res,max_it,@blkprec,[],[],A,B); %,x0);
                end
                INNER = 1;
                X = X(1:(length(X)/2))+1i*X((length(X)/2+1):end);
                ITER = length(RESVEC);
            case {'QMR', 'qmr'}
                [X,FLAG,RELRES,ITER,RESVEC] = qmr(C,c,tol_res,max_it,M0,M1);
                ITER = ITER(end);
                INNER = 1;
            case {'AA-PMHSS'}
                [X,ITER,INNER,RESVEC,RELRES] = AA_PMHSS(C,x0,c,tol_res,max_it,M0,M1,inner_solver,max_inner_iter,timings);
                ITER = ITER(end);
            case {'PMHSS-it', 'pmhss-it', 'PMHSS-ITERATION'}
                [X,ITER,INNER,RESVEC,RELRES]= PMHSS_it(C,x0,c,tol_res,max_it,M0,M1,inner_solver,max_inner_iter,timings);
            otherwise
                disp(['This solver type is not implemented! ', solver_type ])
        end
        times(i) = toc;
    end
    RES = norm(c-C*X);
    RELRES = 1/norm(c)*RES;
    ERROR = norm(xref-X);
    warmstart_iters = repeats / 3;
    times = times(warmstart_iters:end);
    
    TIME = mean(times);
    TIME_STD = std(times);
    STD_PCT = (TIME_STD / TIME) * 100;
    
    if ~isempty(inner_iter)
        INNER = inner_iter;
    end
    
    switch solver_type
        
        case {'AA-PMHSS','PMHSS-it', 'pmhss-it', 'PMHSS-ITERATION'}
            PRECRES = norm((A+B)\(c-C*X));
            format = ('N = %d & %d & %d & %.4g ± %.4g%% & %g & %g & %g \\\\ \n');
            output = [n(1) ITER sum(INNER) TIME STD_PCT RES ERROR PRECRES];
        case {'PMHSS-GMRES', 'phmss-gmres', ...
                'PRESB-GMRES','presb-gmres'}
            PRECRES = RESVEC(end);
            format = ('N = %d & %d & %d & %.4g ± %.4g%% & %g & %g & %g \\\\ \n');
            output = [n(1) ITER sum(INNER) TIME TIME_STD RES ERROR PRECRES];
        otherwise
            format = ('N = %d & %d & %.4g ± %.4g%% & %g & %g & - \\\\ \n');
            output = [n(1) ITER TIME STD_PCT RES ERROR];
    end
    fprintf(log_file, format,output);
    fprintf(format,output);
end

%These are the actual preconditioners, beginning with PRESB
function w = blkprecCG__(v,A,B,max_iter)
    global inner_iter
    H=A+B;
    
    s = size(A(:,1));
    p = v(1:s);
    q = v((s(1)+1):end,1);
    
    [h, ~, ~, iter_cg_1, ~] = pcg(H,(p+q),1e-12,max_iter,[],[]); 
    [y, ~, ~, iter_cg_2, ~] = pcg(H,(q-B*h),1e-12,max_iter,[],[]);
    x = h-y;
    inner_iter = [inner_iter, iter_cg_1 + iter_cg_2];
    w = [x;
         y];
end
%And PMHSS iteration based preconditioner.
function w = pmhss_precCG__(v,A,B,max_iter)
    global inner_iter
    [z, ~, ~, iter_cg, ~] = pcg(A+B,v,1e-12,max_iter,[],[]); 
    % w = 0.5*(1-1i)*z; 
    w = z; 
    inner_iter = [inner_iter, iter_cg];
end
