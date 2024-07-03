function r0discretoedades=r0discretoedades(n,arrayg,arraym,arrayb,arrayk)

    b = sym('b', [n, n]);
    k = sym('k', [1, n]);
    N = sym('N', [1, n]);
    T = sym('T', [1, n]);
    
    V = sym(zeros(n, n));
    for i = 1:n
        V(i, i) = T(i);
        if i < n
            V(i+1, i) =-k(i); 
        end
    end
    F = sym(zeros(n, n));
    for i = 1:n
        for j = 1:n
            F(i, j) = b(j, i) * N(j);  % F_ij = b_ji * N_j
        end
    end

    K=F*V^-1;
    Tarray=arrayg+arraym+arrayk;
    Barray=arraym+arrayk;
    Tarray(n)=arrayg(n)+arraym(n);
    V = subs(V, T, Tarray);
    V = subs(V, k, arrayk);
    bmatrix=repmat(arrayb, n, 1);
    a=1/n;
    F = subs(F, b, bmatrix);
    arrayn=ones(1,n);
    arrayn(1)=a/Barray(1);
    for i=2:n
        arrayn(i)=arrayk(i-1)*arrayn(i-1)/Barray(i);
    end
    F = subs(F, N, arrayn);
    
    V_numeric = double(V);
    F_numeric = double(F);
    
    K = F_numeric / V_numeric;
    
    eigenvalues = eig(K);
    
    r0discreto = max(eigenvalues);
    

    disp(r0discreto);

    r0discretoedades=max(eig(K));
