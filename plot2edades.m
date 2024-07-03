gamma1 = 1; gamma2 = 2;
k1 = 1/0.3; k2 = 1/0.7; %asumiendo que la edad máxima es de 1 unidad temporal
beta = [2.6 1.4; 1.3 0.5];
eta1 = 0.2; eta2 = 0.1;
mu1 = 0.3; mu2 = 0.7;
a = 4;
n1 = a/(mu1 + gamma1);
n2 = n1 * k1 / (gamma2 + mu2);
F = [beta(1,1) * n1 beta(2,1) * n2; beta(1,2) * n1 beta(2,2) * n2];
V = [gamma1 + k1 + mu1 0; -k1 gamma2 + mu2];
M = F / V;
eig(M)

y0 = [0.7, 0.8, 0.3, 0.2, 0, 0];

tspan = [0 200];

odefun = @(t, y) sirs2edadessystem(t, y, gamma1, gamma2, eta1, eta2, beta, n1, n2, a, mu1, mu2, k1);
[t, y] = ode45(odefun, tspan, y0);
s1 = y(:, 1);
s2 = y(:, 2);
i1 = y(:, 3);
i2 = y(:, 4);
ds1 = diff(s1);
di1 = diff(i1);
ds2 = diff(s2);
di2 = diff(i2);

figure;
plot(s1, i1, 'r');
hold on;
quiver(s1(1:end-1), i1(1:end-1), ds1, di1, 'AutoScale', 'on', 'AutoScaleFactor', 1, 'Color', 'r'); % Arrows in red
xlabel('S_1');
ylabel('I_1');
title('Modelo SIRS con dos grupos de edad, R_0>1');
hold off;

figure;
plot(s2, i2, 'r'); 
hold on;
quiver(s2(1:end-1), i2(1:end-1), ds2, di2, 'AutoScale', 'on', 'AutoScaleFactor', 1, 'Color', 'r'); % Arrows in red
xlabel('S_2');
ylabel('I_2');
title('Modelo SIRS con dos grupos de edad, R_0>1');
hold off;


function dydt = sirs2edadessystem(t, y, gamma1, gamma2, eta1, eta2, beta, n1, n2, a, mu1, mu2, k1)
    s1 = y(1); s2 = y(2); i1 = y(3); i2 = y(4); r1 = y(5); r2 = y(6);
    dydt = [
        -beta(1,1)*s1*i1 - beta(2,1)*s1*i2*n2/n1 + eta1*r1 - k1*s1 + a - mu1*s1;
        -beta(1,2)*s2*i1*n1/n2 - beta(2,2)*s2*i2 + eta2*r2 + k1*s1 - mu2*s2;
        beta(1,1)*s1*i1 + beta(2,1)*s1*i2*n2/n1 - gamma1*i1 - k1*i1 - mu1*i1;
        beta(1,2)*s2*i1*n1/n2 + beta(2,2)*s2*i2 - gamma2*i2 + k1*i1 - mu2*i2;
        gamma1*i1 - eta1*r1 - k1*r1 - mu1*r1;
        gamma2*i2 - eta2*r2 + k1*r1 - mu2*r2
    ];
end