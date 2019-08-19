function prac11()
close all;
clear;
clc;

% part_1
% H(Q) = a - b * Q^2
n = 15;
Q = 5000 * rand(n, 1);
Dispersion = 36;
ab_true = [32.5; 1];
Hk = ab_true(1) - (Q.^2).* ab_true(2) + randn(n, 1) * sqrt(Dispersion);

% R, Q = QQ, betta_apriornoe
W = diag(rand(n, 1));
R = W' * W;
W_1 = diag(rand(2, 1));
QQ = W_1'* W_1;
betta_apriori = [32; 2];

% MMNK - evaluation
X = [Q.^0, Q.^2];
K = ((QQ + X'* R * X)^(-1)) * X'* R;
betta_eval = betta_apriori + K * (Hk - X * betta_apriori);

% graphics
dim_a = -500:10:500;
dim_b = -10:1:10;

[Dim_a, Dim_b] = meshgrid(dim_a, dim_b);
grid_matrix = [Dim_a(:) Dim_b(:)];
for i = 1:size(grid_matrix, 1)
    J(i) = functional(Hk, X, betta_apriori, R, QQ, grid_matrix(i, :)');
end
J = reshape(J, length(dim_b), length(dim_a));

Jab = functional(Hk, X, betta_apriori, R, QQ, betta_eval); 

figure(1); 
hold on;
surf(dim_a, dim_b, J); xlabel('a'); ylabel('b');
plot3(betta_eval(1), betta_eval(2), Jab, 'r*');
hold off;
figure(2);
hold on;
contour(dim_a, dim_b, log(log(J))); xlabel('a'); ylabel('b');
plot3(betta_eval(1), betta_eval(2), log(log(Jab)), 'r*');
legend('Contour Lines', 'Solution with MMNK');
hold off;

% part_2
% to quad form
Hk1 = 2 * (X' * R * X + QQ);
fk = -2 * (Hk' * R * X + betta_apriori' * QQ)';

% use of quadprog
abq = quadprog(Hk1, fk);
% limits betta_min, betta_max
abmax = [30; 15];
abmin = [10; 4];
ek = [1 0; 0 1];
A = [ek; -ek];
ab = [abmax; -abmin];
abql = quadprog(Hk1, fk, A, ab);
blim(1:length(dim_a)) = abmin(length(abmin));

figure(3);
hold on;
contour(dim_a, dim_b, log(log(J))); xlabel('a'); ylabel('b');
plot(dim_a, blim, 'b', abql(1), abql(2), '*b');
plot3(abq(1), abq(2), Jab, '*g', betta_eval(1), betta_eval(2), log(log(Jab)), 'or');
legend('Contour Lines', 'Min lim for 2 koef', 'Solution with limits', 'Solution without limits', 'Solution with MMNK');
hold off;

end

function J = functional(H, X, betta_apriori, R, Q, betta_eval)
J = (H - X * betta_eval)' * R * (H - X * betta_eval) + (betta_eval - betta_apriori)' * Q * (betta_eval - betta_apriori);
end
