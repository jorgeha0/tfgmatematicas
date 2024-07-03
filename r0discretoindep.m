function r0discretoindep=r0discretoindep(n,arrayg,arrayd,arrayb)

    b = sym('b', [n, n]);
    k = sym('k', [1, n]);
    N = sym('N', [1, n]);
    G = sym('g', [1, n]);
    
    V = sym(zeros(n, n));
    for i = 1:n
        V(i, i) = G(i); 
    end
    F = sym(zeros(n, n));
    for i = 1:n
        for j = 1:n
            F(i, j) = b(j, i);  
        end
    end

    K=F*V^-1;
    V = subs(V, G, arrayg);
    bmatrix=repmat(arrayb, n, 1);
    F = subs(F, b, bmatrix);
    
    V_numeric = double(V);
    F_numeric = double(F);
    
    K = F_numeric / V_numeric;
    
    eigenvalues = eig(K);
    
    r0discretoindep = max(eigenvalues);
    

    disp(r0discretoindep);

