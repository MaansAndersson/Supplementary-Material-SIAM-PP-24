% clc
%%
log_names = ["log_padebai_max_inner.txt", "log_helmholtz_max_inner.txt", "log_dfda_max_inner.txt"];
%log_names = ["slask_1.txt", "slask_2.txt", "slask_3.txt"];
matrix_sizes = [100, 200, 300];
disp('-----')
repeats = 15;
for i = [1 2 3]
    log_file = fopen(log_names(i), 'w');
    for m = matrix_sizes
        max_iter = round(m);
        rel_tol = 1e-8;
        
        switch i
            case 1
                [C, c] = PadeBai(m);
            case 2
                [C, c] = IHelmholtz(m,0,0.02);
            case 3
                [C, c] = DFDA(m,pi,0.02);
        end
        
        
        % tic
        xref=C\c;
        % toc
        xin = 0*c;
        
        fprintf(log_file, "------ AA-PMHSS --------\n");
        fprintf("------ AA-PMHSS --------\n");
        [RESVEC0, INNER_AA] = complexsolvetimer('AA-PMHSS',max_iter,rel_tol,C,c,xref,[],[],xin,'cg',m*m, repeats, log_file);
        fprintf(log_file, "------ AA-PMHSS (m) --------\n");
        fprintf("------ AA-PMHSS (m)--------\n");
        [RESVEC0, INNER_AA] = complexsolvetimer('AA-PMHSS',max_iter,rel_tol,C,c,xref,[],[],xin,'cg',m, repeats, log_file);

        fprintf(log_file, "------ PMHSS --------\n");
        fprintf("------ PMHSS --------\n");
        [RESVEC0, INNER_PMHSS] = complexsolvetimer('PMHSS-ITERATION',max_iter,rel_tol,C,c,xref,[],[],xin,'cg',m*m, repeats, log_file);
        fprintf(log_file, "------ PMHSS (m)--------\n");
        fprintf("------ PMHSS (m)--------\n");
        [RESVEC0, INNER_PMHSS] = complexsolvetimer('PMHSS-ITERATION',max_iter,rel_tol,C,c,xref,[],[],xin,'cg',m, repeats, log_file);

        fprintf(log_file, "------ PMHSS-GMRES --------\n");
        fprintf("------ PMHSS-GMRES --------\n");
        [RESVEC_PMHSSG, INNER_PMHSS_GMRES] = complexsolvetimer('PMHSS-GMRES',max_iter,rel_tol,C,c,xref,[],[],xin,'cg',m*m, repeats, log_file);
        fprintf(log_file, "------ PMHSS-GMRES (m)--------\n");
        fprintf("------ PMHSS-GMRES (m)--------\n");
        [RESVEC_PMHSSG, INNER_PMHSS_GMRES] = complexsolvetimer('PMHSS-GMRES',max_iter,rel_tol,C,c,xref,[],[],xin,'cg',m, repeats, log_file);
        
        fprintf(log_file, "------ PRESB-GMRES --------\n");
        fprintf("------ PRESB-GMRES --------\n");
        [RESVEC_PRESB, INNER_PRESB_GMRES] = complexsolvetimer('PRESB-GMRES',max_iter,rel_tol,C,c,xref,[],[],xin,'cg',m*m, repeats, log_file);
        fprintf(log_file, "------ PRESB-GMRES (m) --------\n");
        fprintf("------ PRESB-GMRES (m)--------\n");
        [RESVEC_PRESB, INNER_PRESB_GMRES] = complexsolvetimer('PRESB-GMRES',max_iter,rel_tol,C,c,xref,[],[],xin,'cg',m, repeats, log_file);
        
        fprintf(log_file, "------ GMRES --------\n");
        fprintf("------ GMRES --------\n");
        [RESVEC_GMRES, INNER0] = complexsolvetimer('gmres',max_iter*10,rel_tol,C,c,xref,[],[],xin,'',0, repeats, log_file);
    end
    
    fclose(log_file);
end

