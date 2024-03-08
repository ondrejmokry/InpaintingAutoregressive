function [Basis,sparsity_init,sparsity_final] = basis_opt_new(X_tr,level_init,epsilon)
 
 %%	Optimizes Unitary Basis 
 %%
 %%  Copyright: (c) Georg Taub�ck, 2006-2010; adapted 2020
 %
 % Basis              ... orthonormal basis that sparsifies all
 %                    training signals X_tr; represented by unitary NxN matrix
 % sparsity_init  ... l1-norm of (all) coefficients using the identity
 %                     matrix as sparsifying matrix
 % sparsity_final ... l1-norm of (all) coefficients when all training signals are
 %                    multiplied by orthonormal basis Basis generated by
 %                    algorithm
 % X_tr           ... NxM matrix containing M training signals as vectors
 % 
 %                    Goal: Y_tr := U*X_tr should be as sparse as possible: minimize
 %                    norm(vec(Y_tr),1) over set of all unitary matrices U
 %
 % epsilon        ... threshold after which algorithm is terminated (if level falls below epsilon)
 % level_init     ... start value of level (see epsilon)
 %

 % initalization

[N,M]=size(X_tr);

I=eye(N);

Aopt=I;
sparsity_old=Inf;
sparsity=norm(vec(X_tr),1); 
sparsity_init=sparsity;
level=level_init;

cnt=1;
cnt_max=20; %maximum number of iterations
off_diags=1; % determines number of off-diagonals of Hermitian matrix (off_diags=1 corresponds to main diagonal plus 1 upper and 1 lower off-diagonal)

 % cvx_solver sdpt3 % only this solver seems to work

 % fprintf('Starting basis optimization ...\n');

while level>epsilon
    level_power=log2(level);
    while sparsity < sparsity_old     
        Y=Aopt*X_tr;
     
        cvx_begin quiet
            variable A(N,N) hermitian banded(off_diags,off_diags); 
            minimize(norm(vec((I+j*2*pi*A)*Y),1));  % this correspods to an averaged l1-norm over all training signals; 
            %minimize(norm((I+j*2*pi*A)*Y,1)); %the maximum of the l1-norms of the training signals could be also possible (but is extremely slow) 
            subject to
                max(max(abs(A))) <= level % this could be replaced by another norm
        cvx_end     
        
     
        sparsity_old_save=sparsity_old;
        sparsity_old=sparsity;
        Aopt_old=Aopt;
        Aopt=expm(j*2*pi*A)*Aopt;
        sparsity=norm(vec(Aopt*X_tr),1);
        cnt=cnt+1;
        if cnt>cnt_max
            break;
        end    
    end
    if cnt>cnt_max
            break;
        end    
    level=level/2;
    Aopt=Aopt_old;
    sparsity=sparsity_old;
    sparsity_old=sparsity_old_save;
end
sparsity_final=norm(vec(Aopt*X_tr),1);

Basis=Aopt;

%end of function
