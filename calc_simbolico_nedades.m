%calcular de forma general la K y los autovalores en el caso de
%compartimentos de edad (con factor envejecimiento)

n=3; %numero de edades
%construccion de variables
syms N real
g = sym('g', [1, n], 'real'); %gamma
A = sym('A', [1, n], 'real'); %valores asintóticos de poblaciones
B=sym('B',[n n], 'real'); %beta
F=sym('F',[n n], 'real'); %matriz F
T=sym('T',[1 n], 'real'); %factores T
k=sym('k',[1 n], 'real');
assume(A >= 0);
assume(B >= 0);
assume(F >= 0);
assume(T >= 0);
assume(k >= 0);


for i=1:n
    for j=1:n
        F(i,j)=A(j)*B(j,i);
    end
end
i=0; j=0;
V=diag(T);
for i=1:n
    for j=1:n
        if j==i-1
            V(i,j)=-k(i-1);
        end
    end
end

K=F*V^-1 %matriz NGM
L=eig(K) %autovalores de K