function PoissonTest();

outer = 0;

for NN = 4:8
    outer = outer+1;
    
    % setting up mesh size
    nn = 2^NN-1;
    
    % storing information for results table
    MatrixSize(outer) = nn^2;
    GridSize(outer) = nn;
    
    % Setting up cell array for fine and coarse operators
    A = cell(1,NN-1);
    
    levels = NN-1;
    levelnum = 1;
    iii = 0;
    tic
    % Creating sequence of coarse grid operatos
    for ii = NN:-1:2
        nN = 2^ii-1;
        iii = iii+1;
        A{iii} = laplace(nN);
    end
    
    % Setting starting value and rhs vector
    x0 = zeros(nn^2,1);
    x = x0;
    b = A{1}*2*ones(nn^2,1);
    
    % Running CG with V-cycle preconditioner
    [~,~,~,ITER] = pcg(A{1},b,1e-8,2000,@(b) MGV2D(A,b,x0,levels,levelnum),[],x0);
    mgtime = toc;
    
    
    tic
    % Computing incomplete cholesky factorisation
    R= ichol(A{1});
    
    % Running CG with and incomplete cholesky preconditioner
    [~,~,~,ITER1] = pcg(A{1},b,1e-8,20000,R',R,x0);
    ilutime = toc;
    
    tic
    A{1}\b;
    DIRECTtime(outer) =toc;
    
    
    MGiters(outer) = ITER;
    ILUiters(outer) = ITER1;
    MGtime(outer) = mgtime;
    ILUtime(outer) = ilutime;
end
data = [GridSize(:),MatrixSize(:),MGiters(:),MGtime(:),ILUiters(:),ILUtime(:),DIRECTtime(:)];

% Set up some options
tblOpts = {'header',{'Grid Size','Matrix Size',...
    'iters','time','iters','time','time'},'format',{'%1.0i','%1.0i','%1.0i'...
    ,'%2.4f','%1.0i','%2.4f','%2.4f'},'align','center','delim','|',...
    'printRow',true};

for ii = 1:size(data,1);
    table('',data(1:ii,:),tblOpts{:}...
        ,'finalRow',ii == size(data,1));
end


    function [x,res] = MGV2D(A,b,x,level,levelnum)
        
        %defining damping parameter
        theta = 3/4;
        n = length(b);
        
        coarsest = 1;
        e1 = ones(n,1);
        if level == coarsest
            % Preform exact solve on coarsest grid
            x = A{levelnum}\b;               
            r = b - A{levelnum}*x;
            
        else 
            
            % relax using damped Jacobi
            Dv = (diag(A{levelnum}));
            r = b - A{levelnum}*x;
            for i=1:3
                x = x + r./Dv;
                r = b - A{levelnum}*x;
            end
            % restrict residual from r to rc on coarser grid
            N = sqrt(length(b));
            r = reshape(r,N,N);
            N1 = (N+1)/2 - 1; n = N1^2;   
            
            % rectriction in a matrix-free way (same as rc = P'*r)
            rc = r(2:2:N-1,2:2:N-1) + .5*(r(3:2:N,2:2:N-1)+r(1:2:N-2,2:2:N-1) +...
                r(2:2:N-1,3:2:N)+r(2:2:N-1,1:2:N-2)) + .25*(r(3:2:N,3:2:N)+...
                r(3:2:N,1:2:N-2)+r(1:2:N-2,3:2:N)+r(1:2:N-2,1:2:N-2));
            rc = reshape(rc,n,1);
            vc=sparse( n,1);
            levelnum = levelnum +1;
            
            % recursively applying two-grid cycle
            [vc,r] = MGV2D(A,rc,vc,level-1,levelnum);
            levelnum = levelnum-1;
            % prolongate/interpolate correction from vc to v on finer grid
            v = sparse(N,N);
            vc = reshape(vc,N1,N1);
            v(2:2:N-1,2:2:N-1) = vc;
            vz = [sparse(1,N);v;sparse(1,N)];
            vz = [sparse(N+2,1),vz,sparse(N+2,1)];
            % prolongate in a matrix-free way 
            v(1:2:N,2:2:N-1) = .5*(vz(1:2:N,3:2:N)+vz(3:2:N+2,3:2:N));
            v(2:2:N-1,1:2:N) = .5*(vz(3:2:N,1:2:N)+vz(3:2:N,3:2:N+2));
            v(1:2:N,1:2:N) = .25*(vz(1:2:N,1:2:N)+vz(1:2:N,3:2:N+2)+...
                vz(3:2:N+2,3:2:N+2)+vz(3:2:N+2,1:2:N));
            n = N^2;
            x = x + reshape(v,n,1);
            
            % relax using damped Jacobi
            Dv = diag(A{levelnum});
            for i=1:3
                r = b - A{levelnum}*x;
                x = x + theta*r./Dv;
            end
            
        end
        res = norm(b - A{levelnum}*x);
        
    end


    function [A] = laplace(n)
        % function to produce 2D laplacian matrix
        
        e1 = ones(n,1);
        A = spdiags([-e1 2*e1 -e1],[-1:1],n,n);
        
        A = kron(A,speye(n)) + kron(speye(n),A);
    end

end