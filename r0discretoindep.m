function r0discretoindep=r0discretoindep(n,arrayg,arrayd,arrayb)

    % Define symbolic variables for b_ij parameters, k_i, N_i, and T_i coefficients
    b = sym('b', [n, n]);
    k = sym('k', [1, n]);
    N = sym('N', [1, n]);
    G = sym('g', [1, n]);
    
    % Construct the V matrix
    V = sym(zeros(n, n));
    for i = 1:n
        V(i, i) = G(i);  % Diagonal elements are T_i
    end
    % Construct the F matrix
    F = sym(zeros(n, n));
    for i = 1:n
        for j = 1:n
            F(i, j) = b(j, i);  % F_ij = b_ji * N_j
        end
    end
    
    % Display the results
    disp(V);
    disp(F);
    K=F*V^-1;
    size(G)
    size(arrayg)
    % Substitute numerical values into the symbolic matrices
    V = subs(V, G, arrayg);
    bmatrix=repmat(arrayb, n, 1);
    F = subs(F, b, bmatrix);
    
    % Convert symbolic matrices to numeric matrices
    V_numeric = double(V);
    F_numeric = double(F);
    
    % Compute the matrix K
    K = F_numeric / V_numeric;
    
    % Calculate the eigenvalues of K
    eigenvalues = eig(K);
    
    % Find the maximum eigenvalue
    r0discretoindep = max(eigenvalues);
    
    % Display the results

    disp(r0discretoindep);

