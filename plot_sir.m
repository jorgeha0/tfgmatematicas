gamma =1;
beta = 4;

s = linspace(0, 1, 20);
i = linspace(0, 1, 20);
[X1, X2] = meshgrid(s, i);

valid_points = (X1 + X2) < 1;

X1_valid = X1(valid_points);
X2_valid = X2(valid_points);

dx = arrayfun(@(x, y) sirssystema([x, y], gamma, beta), X1_valid, X2_valid, 'UniformOutput', false);
dx1 = cellfun(@(x) x(1), dx);
dx2 = cellfun(@(x) x(2), dx);

h = quiver(X1_valid, X2_valid, dx1, dx2, 'AutoScale', 'on', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);

h.Color = [0 0.4470 0.7410];

xlabel('S');
ylabel('I');
axis([0 1 0 1]);
title('Diagrama de fases del modelo SIR, \beta=4, \gamma=1'); hold on;
grid on; 

function sirssystem = sirssystema(y, gamma, beta)
    sirssystem = [-y(1) * beta * y(2), y(1) * beta * y(2) - gamma * y(2)];
end