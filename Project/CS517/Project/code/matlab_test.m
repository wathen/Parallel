clear
close all
load('matrix/Wdim.mat');
addpath('IFIS_AMG/')

for i = 1:5
    fprintf('   %1.0f  \n',i)
    Amat = sprintf('matrix/A%i.mat',Wdim(i));
    Pmat = sprintf('matrix/P%i.mat',Wdim(i));
    Bmat = sprintf('matrix/B%i.mat',Wdim(i));
    Umat = sprintf('matrix/U%i.mat',Wdim(i));
    A = load( Amat );
    P = load (Pmat);
    b = load(Bmat);
    u = load(Umat);
        
    A = struct2cell(A);
    A = A{1};        
    P = struct2cell(P);
    P = P{1};
    b = struct2cell(b);
    b = b{1}';
    u = struct2cell(u);
    u = u{1};

    r = ichol(P);
    
    AA = (r'\A)/r;
    b = r'\b;
%     AA = P\A;
%     b = P\b;
    levels = 10;
    
    [grid_data] = amg_grids_setup(AA,1,levels);
    [smoother_params] = amg_smoother_params(grid_data,'ILU',1);
    [smoother_data] = amg_smoother_setup(grid_data, smoother_params);
    x = sparse(length(b),1);
    res = norm(AA*x-b)/norm(b);
    jj = 0;
    M = @(x) amg_v_cycle(b, grid_data, smoother_data,x);
%     [X,FLAG,RELRES,ITER = pcg(AA,b,1e-10,1000,M);
    
    while res > 1e-6
        [x] = amg_v_cycle(b, grid_data, smoother_data,levels, x);
%         [x,res] = MGV3D(A,b,x,levels,jj);
        res = norm(AA*x-b)/norm(b)
        jj = jj+1;
    end
    
    fprintf('\n  N = %2i: residual = %1.4e after %1.0i iterations\n\n',Wdim(i),res,jj);
    norm(u-r\x,2)
    
end
