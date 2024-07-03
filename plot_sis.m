gamma = 1;
beta = 1/4;
K = 1 - gamma / beta;
r = beta - gamma;
I0 = 0.2;

f = @(x) K * I0 ./ (I0 + (K - I0) .* exp(-r .* x));

I = linspace(0, 1, 20);
S = 1 - I;

dI = beta * I .* S - gamma * I;
dS = -dI;

figure;
quiver(S, I, dS, dI, 'AutoScale', 'on', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
hold on;

h_point = plot(1,0, 'bo', 'MarkerFaceColor', 'r', 'MarkerSize', 6);

xlabel('Susceptibles (S)');
ylabel('Infectados (I)');
title('Diagrama de fases del modelo SIS, \beta=1/4, \gamma=1');
legend(h_point, 'Punto (1, 0)');
axis([0 1 0 1]);
grid on; 
hold off;
