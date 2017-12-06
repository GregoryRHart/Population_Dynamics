function FIM_Ising

N = 10000;      % # of sequences
m = 5;          % # of positions
tol = eps;      % sparsity tolerance in FIM elements
neigs = 4;      % # largest magnitude eigenvectors of FIM to compute

lambda_h = 0;   % regularization parameter on {h}
lambda_J = 0;   % regularization parameter on {J}
K = N;          % effective number of sequences



% errorchecking inputs
if neigs > m+nchoosek(m,2)
    neigs = m+nchoosek(m,2);    % setting # evecs to return to size of FIM
end


% generating a random multiple sequence alignment (MSA)
%   MSA(i,j) = 0 => wild-type amino acid in position j of sequence i 
%   MSA(i,j) = 1 => mutant amino acid in position j of sequence i 
MSA = dlmread('msa.txt') 


% constructing sparse Fisher information matrix (FIM) 
FIM = sparse(m+nchoosek(m,2),m+nchoosek(m,2))


% computing FIM matrix elements
for i=1:m
    for j=i:m
        for k=1:m
            for l=k:m

                % error checking row indexing on (i,j)
                if (i<1 || i>m)
                    error('i not in range [1,m]');
                end
                if (j<1 || j>m)
                    error('j not in range [1,m]');
                end
                if (j < i)
                    error('j not in range [i,m]');
                end

                % computing row associated with (i,j)
                row = (i-1)*m - (i-1)*(i-2)/2 + (j-i) + 1;


                % error checking col indexing on (k,l)
                if (k<1 || k>m)
                    error('k not in range [1,m]');
                end
                if (l<1 || l>m)
                    error('l not in range [1,m]');
                end
                if (l < k)
                    error('l not in range [k,m]');
                end

                % computing col associated with (k,l)
                col = (k-1)*m - (k-1)*(k-2)/2 + (l-k) + 1;


                % error checking on row and col
                if (col < row)
                    %error('col < row - no need to compute lower triangle of FIM since symmetric');
                    continue;
                end
                
                % printing (i,j,k,l) and (row,col) to screen
                %fprintf('(i,j,k,l) = (%d,%d,%d,%d)\n',i,j,k,l);
                %fprintf('(row,col) = (%d,%d)\n',row,col);
                
                % computing FIM(row,col)
                P_ijkl = mean( MSA(:,i) .* MSA(:,j) .* MSA(:,k) .* MSA(:,l) );
                P_ij = mean( MSA(:,i) .* MSA(:,j) );
                P_kl = mean( MSA(:,k) .* MSA(:,l) );
                
                if (i == k) && (j == l)
                    if (i == j)
                        FIM_row_col = -(P_ij*P_kl - P_ijkl - 2*lambda_h/K);
                    else
                        FIM_row_col = -(P_ij*P_kl - P_ijkl - 2*lambda_J/K);
                    end
                else
                    FIM_row_col = -(P_ij*P_kl - P_ijkl);
                end
                
                if FIM_row_col > tol
                    FIM(row,col) = FIM_row_col;
                end
                
            end
        end
    end
end

% completing bottom triangle by symmetry
FIM = FIM + tril(transpose(FIM),-1);

FIM

% diagonalizing
[V,D,flag] = eigs(FIM,neigs,'LM')

