function A2Q2()
clc

NN = [100,128,256,512,1024,2048,4096];
kk = [NN/2;NN];

for JJ = 1:size(kk,1)
    for II = 1:length(NN)
        
        % Defining descrete Laplacian
        N = NN(II);
        A = Laplacian(N);
        
        % Defining eigenvalues of descrete Laplacian
        Eigen = @(ii,jj,N) 4-2*(cos(ii*pi/(N+1))+cos(jj*pi/(N+1)));
        SmallestEigA = [Eigen(1,1,N);Eigen(2,1,N);Eigen(1,2,N)];
        LargestEigA = [Eigen(N,N,N);Eigen(N-1,N,N);Eigen(N,N-1,N)];
        
        b = randn(N^2,1);
        
        % Calling Lanczos code
        k = kk(JJ,II);
        [T] = lancz(A, b, k);
        
        % Calculating eigenvalues of T
        OPTS.maxit = 1e6;
        OPTS.tol = 1e-10;
        
        EIG = eig(full(T));
        EIG = sort(EIG);
        SmallestEigT = EIG(1:3);
        LargestEigT =  EIG(end-2:end);
        
        % Sort eigenvalues
        SSEA = sort(SmallestEigA);
        SLEA = sort(LargestEigA);
        SSET = sort(SmallestEigT);
        SLET = sort(LargestEigT);
        
        % Defining table data
        data = [SSEA,SSET,abs(SSEA-SSET)./abs(SSEA),SLEA,SLET,abs(SLEA-SLET)./abs(SLEA)];
        
        % Set up some options
        tblOpts = {'header',{'Smallest Eig A','Smallest Eig T',...
            'inf-norm rel','Largest Eig A','Largest Eig T'...
            ,'inf-norm rel'},'format',{'%1.4e','%1.4e','%1.4f'...
            ,'%1.6f','%1.6f','%1.4e'},'align','center','delim','|',...
            'printRow',true};
        
        for ii = 1:size(data,1);
            table(['Table of Eigenvalues for n = ',num2str(NN(II)^2),...
                ' and k = ',num2str(k)],data(1:ii,:),tblOpts{:}...
                ,'finalRow',ii == size(data,1));
        end
        
    end
    fprintf('\n\n\n\n\n')
end

    function [A] = Laplacian(n)
        % Creating discretised Laplacian
        
        e = ones(n,1);
        
        % Creating sparse diagonal matrices
        I = spdiags(e,0,n,n);
        I1 =spdiags(e,1,n,n);
        I2 = spdiags(e,-1,n,n);
        
        % Creating 1D Convection-Diffusion matricies
        A1D = 2*I -1*I1 - 1*I2;
        
        % Creating 2D Convection-Diffusion matrix
        A = kron(I,A1D)+kron(A1D,I);
        
    end



    function [T,Q] = lancz(A, b, k)
        %function [T,Q] = lancz(A, b, k)
        %
        % Function the performs the Lanczos process
        %
        % Input:
        %        A - Symmetic matrix
        %        b - initial guess
        %        A - number of steps in the Lanczos algorithm
        %
        % Output:
        %        T - Symmetic Hessenberg matrix (Tridiagonal)
        %        Q - (OPTIONAL) orthogonal basis
        
        n = length(b);
        qprev = sparse(n,1);
        q = b / norm(b);
        beta = [];
        alpha = [];
        
        if nargout == 2
            Q = [];
        end
        
        for i = 1:k
            v = A*q;
            alpha(i) = q' * v;
            
            if i == 1
                v = v - alpha(i)*q;
            else
                v = v - beta(i-1)*qprev - alpha(i)*q;
            end
            beta(i) = norm(v);
            qprev = q;
            
            if nargout == 2
                Q = [Q,q];
            end
            
            if (abs(beta(i)) < 1e-10)
                break
            end
            q = v / beta(i);
        end
        beta = beta(:);
        T = spdiags([beta alpha(:) [0;beta(1:end-1)]],[-1:1],i,i);
    end


end
