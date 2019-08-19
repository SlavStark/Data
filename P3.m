function prac3()
close all;
clear;
clc;
warning off;

K_F = 8; % êîýô óñèëåíèÿ
N = 400;
Rk = 10;
Qk = 1;
P = 50;
SP = 10;
A_1 = [-4 -15 -0.4]; % êîðíè õàðàêòåðèñòè÷åñêîãî ìíîãî÷ëåíà
p = poly(A_1); % êîýô ïîëèíîìà
dt = 0.1; % ïåðèîä äèñêðåòèçàöèè dt
% äèñêðåòíàÿ ìîäåëü
Ad = [1 - p(2)/p(1) * dt -p(3)/p(1) * dt -p(4)/p(1) * dt; dt 1 0; 0 dt 1];
Bd = [K_F/p(1) * dt; 0; 0];
Cd = [0 0 1];

% part 1
M = get_M_matrix(Ad, Bd, Cd, N);
L = get_L_matrix(Ad, Cd, N);
u = [[0] [ones(1, N-1)]]';  % åäèíè÷íîå ñòóïåí÷àòîå âîçäåéñòâèå
time = 0:dt:dt*(size(u, 1)-1);
NY = zeros(size(Cd, 2), 1);  % íà÷àëüíûå óñëîâèÿ
Y_1 = L * NY + M * u;    % ðåàêöèÿ îáúåêòà ïî ôîðìóëå èç 2 ïðàêòèêè

% ïåðåâîä â ìîäåëü ñ ïðèðàùåíèÿìè
Az = [Ad zeros(3, 1); Cd*Ad 1];
Bz = [Bd; Cd*Bd];
Cz = [0 0 0 1];
v = [[1] [zeros(1, N - 1)]]';
Mz = get_M_matrix(Az, Bz, Cz, N);
Lz = get_L_matrix(Az, Cz, N);
p = zeros(3 + 1, 1);
Y_2 = Lz * p + Mz * v;

figure(1);
hold on;
plot(time, Y_1, 'g', time, Y_2, 'r --');
grid on;
title('Ðåàêöèÿ íà ñêà÷îê');
legend('èñõîäíàÿ ìîäåëü', 'ìîäåëü â ïðèðàùåíèÿõ');
hold off;

% part 2
Q = get_Q_matrix(Qk, P);
R = get_R_matrix(Rk, P);
Mz = get_M_matrix(Az, Bz, Cz, P);
Lz = get_L_matrix(Az, Cz, P);
[Xk, Yk, u_mpc, v_mpc, v, pk] = get_zeros_matrix(N, 3);
for k = 2 : N
v_mpc(k) = mpc_calc_u(R, Q, pk(:, k - 1), Lz, Mz, SP);    % óïðàâëÿþùåå âîçäåéñòâèå MPC - ðåãóëÿòîðà
u_mpc(k) = u_mpc(k - 1) + v_mpc(k);
Xk(:,k) = Ad * Xk(:,k-1) + Bd * u_mpc(k)';   % íîâîå ñîñòîÿíèå ñèñòåìû
Yk(k) = Cd * Xk(:, k);  % âûõîä ñèñòåìû
pk(:,k ) = [Xk(:, k) - Xk(:, k - 1); Yk(k)];
end

% part 3
[Xk, Yk_qp, u_qp, v_qp, v_qp, pk] = get_zeros_matrix(N, 3);
for k = 2 : N
v_qp(k) = mpc_quadprog(pk(:,k-1), SP, R, Q, Lz, Mz);
u_qp(k) = u_qp(k - 1) + v_qp(k);
Xk(:,k) = Ad * Xk(:,k-1) + Bd * u_qp(k)';   % íîâîå ñîñòîÿíèå ñèñòåìû
Yk_qp(k) = Cd * Xk(:, k);  % âûõîä ñèñòåìû
pk(:,k ) = [Xk(:, k) - Xk(:, k - 1); Yk_qp(k)];
end

% figure;
% subplot(2, 1, 1);
% plot(time, Yk, time, ones(1, N)*SP, time, Yk_qp, '--');
% subplot(2, 1, 2);
% plot(time, u_mpc, time, u_qp, '--');

% quadprog ñ îãðàíè÷åíèÿìè
u_min = -1; u_max = 35;
v_min=-0.2; v_max=0.75;
y_min = 0; y_max = 9;
[Xk, Yk_qp_c, u_mpc_c, v_mpc_c, v, pk] = get_zeros_matrix(N, 3);
for k = 2 : N
v_mpc_c(k) = mpc_constrained(y_min, y_max, u_min, u_max, v_min, v_max, pk(:,k-1), 4, R, Q, Lz, Mz, P, u_mpc_c(k-1), SP);
u_mpc_c(k) = u_mpc_c(k - 1) + v_mpc_c(k);
Xk(:,k) = Ad * Xk(:,k-1) + Bd * u_mpc_c(k)';   % íîâîå ñîñòîÿíèå ñèñòåìû
Yk_qp_c(k) = Cd * Xk(:, k);  % âûõîä ñèñòåìû
pk(:,k ) = [Xk(:, k) - Xk(:, k - 1); Yk_qp_c(k)];
end

plotting(2, time, Yk, N, SP, Yk_qp, Yk_qp_c, u_mpc, u_qp, u_mpc_c)

end

function R = get_R_matrix(Rk, P)

v = zeros(1, size(Rk, 1)*P);
for i = 0 : P - 1
    v(1, 1+size(Rk)*i : size(Rk)*(i+1)) = diag(Rk);
end
R = diag(v);

end

function Q = get_Q_matrix(Qk, P)

v = zeros(1, size(Qk, 1)*P);
for i = 0 : P - 1
    v(1, 1+size(Qk)*i : size(Qk)*(i+1)) = diag(Qk);
end
Q = diag(v);

end

function L = get_L_matrix(A, C, P) 

L = zeros(P, size(A, 1));
for i = 1 : P
    L(i, :) = C * A^i;
end

end

function M = get_M_matrix(A, B, C, P) 

M = zeros(P, P);
for i = 1 : P
    for j = 1 : i
        M(i, j) = C * A^(i - j) * B;
    end
end

end

function [Xk, Yk, u_mpc, v_mpc, v, pk] = get_zeros_matrix(N, n)

Xk = zeros(n, N);
Yk = zeros(N, 1);
u_mpc = zeros(N, 1);
v_mpc = zeros(N, 1);
v = zeros(N, 1);
pk = zeros(n + 1, N);

end

function u_mpc = mpc_calc_u(R, Q, Xk, L, M, SP)

u_mpc = (-inv(Q + M' * R * M) * M' * R) * (L * Xk - SP);
u_mpc = u_mpc(1);

end

function v = mpc_quadprog(pk, SP, R, Q, L, M)

f = (L * pk - SP)' * R * M;    
H = M' * R * M + Q;
v = quadprog(H, f);
v = v(1);

end

function v = mpc_constrained(y_min, y_max, U_min, U_max, v_min, v_max, p_k, m, R, Q, L, M, P, Uk_1, SP)

[Au, Bu] = control_constraint(P, U_min, U_max, Uk_1);
[Av, Bv] = speed_limit(v_min, v_max, P);
[Ay, By] = regulatory_limit(M, y_min, y_max, L, p_k);

A = [Av; Au; Ay];
B=[Bv; Bu; By];
H = M' * R * M + Q;
f = (L * p_k - SP)' * R * M;
v = quadprog(H, f, A, B);
v = v(1);

end

function [Au, Bu] = control_constraint(P, U_min, U_max, Uk_1)

% îãðàíè÷åíèÿ íà óïðàâëÿþùåå âîçäåéñòâèå
Mu = get_M_matrix(1, 1, 1, P);
Au = [Mu; - Mu];
M0 = ones(P,1);
U_max = U_max * ones(P, 1);
U_min = U_min * ones(P,1);  
Bu = [U_max - M0 * Uk_1; - U_min + M0 * Uk_1];

end

function [Av, Bv] = speed_limit(v_min, v_max, P)

%îãðàíè÷åíèÿ íà ñêîðîñü èçìåíåíèÿ u
v_min = ones(P, 1) * v_min;
v_max = ones(P, 1) * v_max;
Av = [eye(P); - eye(P)];
Bv = [v_max; - v_min];

end

function [Ay, By] = regulatory_limit(M, y_min, y_max, L, p_k)

% îãðàíè÷åíèÿ íà ðåãóëèðóåìóþ âåëè÷èíó
Ay = [M; - M];
By = [y_max - L * p_k; - y_min + L * p_k];

end

function plotting(NumFig, time, Yk, N, SP, Yk_qp, Yk_qp_c, u_mpc, u_qp, u_mpc_c)

figure(NumFig);
subplot(2, 2, 1);
plot(time, Yk, 'g', time, ones(1, N)*SP, 'b', time, Yk_qp, 'r --');
grid on;
title('Ïåðåõîäíûé ïðîöåñ');
legend('MPC', 'Óñòàâêà', 'MPC quadprog áåç îãðàí');
subplot(2, 2, 2);
plot(time, Yk_qp_c);
grid on;
title('Ïåðåõîäíûé ïðîöåñ');
legend('MPC quadprog ñ îãðàí');
subplot(2, 2, 3);
plot(time, u_mpc, 'g', time, u_qp, 'r --');
grid on;
title('Óïðàâëÿþùåå âîçäåéñòâèå');
legend('MPC', 'MPC quadprog áåç îãðàí');
subplot(2, 2, 4);
plot(time, u_mpc_c);
grid on;
title('Óïðàâëÿþùåå âîçäåéñòâèå');
legend('MPC quadprog ñ îãðàí');

end