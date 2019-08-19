function prac2()
clear;
clc;
station_number = 4;
pump_number = 4;
[Q, H, N] = read_pump_data(station_number, pump_number);
SourcePicture(Q, H, 1);
SourcePicture(Q, N, 2);
PolyCount(Q, H, 3, 1);
PolyCount(Q, N, 3, 2);

%H(Q) = a - b * Q^2
n = length(H);
betta = rand(2, 1);
for i_2 = 1:n
    X(i_2, :) = [Q(i_2, 1)^0,-Q(i_2, 1)^2];    
end
e_betta = H - X * betta;
J_betta = sum((e_betta') * e_betta);
betta_optimal = (((X') * X)^(-1)) * (X') * H;
reg_opt_H = regress(H, X);
a = betta_optimal(1);
b = betta_optimal(2);
H_end = a - b * Q.^2;
e1 = H - X * betta_optimal;
Q_Range1 = (0:1.5 * max(Q));
H_end1=a - b * Q_Range1.^2;
figure(1);
hold on;
plot(Q_Range1, H_end1, 'g:', 'DisplayName', 'H(Q)2.2', 'LineWidth', 3);
hold off;
figure(3);
subplot(5, 1, 5);
plot(Q, e1, 'o', 'DisplayName', 'e for H(Q)2.2');
figure(5);
subplot(5, 1, 5);
hist(e1);
title('e for H(Q)2.2');
%N(Q) = k1 * Q - k2 * Q^2
n_N = length(N);
betta_N = rand(2, 1);
for i_3 = 1:n_N
    X_N(i_3, :) = [Q(i_3, 1)^1,-Q(i_3, 1)^2];
end
e_betta_N = N - X_N * betta_N;
J_betta_N = sum((e_betta_N') * e_betta_N);
betta_optimal_N = (((X_N') * X_N)^(-1)) * (X_N') * N;
reg_opt_N = regress(N, X_N);
k1 = betta_optimal_N(1);
k2 = betta_optimal_N(2);
N_end = k1 * Q - k2 * Q.^2;
N_end1=k1 * Q_Range1 - k2 * Q_Range1.^2;
e2 = N - X_N * betta_optimal_N;
figure(2);
hold on;
plot(Q_Range1, N_end1, 'g:', 'DisplayName', 'N(Q)2.2', 'LineWidth', 3);
hold off;
figure(4);
subplot(5, 1, 5);
plot(Q, e2, 'o', 'DisplayName', 'e for N(Q)2.2');
figure(6);
subplot(5, 1, 5);
hist(e2);
title('e for N(Q)2.2');

for i_f = 1:6
    figure(i_f);
    legend('show');
end
end

function [Q, H, N] = read_pump_data(station_number, pump_number)
filename = ['.\ÍÀ_ñòàò_csv\' get_pump_string(station_number, pump_number) '.csv'];
fid = fopen(filename);
data = textscan(fid, '%s%s%s', 'delimiter', ';');
fclose(fid);
% Convert ',' to '.'
data = cellfun( @(x) str2double(strrep(x, ',', '.')), data, 'uniformoutput', false);
data = cell2mat(data);
Q = data(:, 1);
H = data(:, 2);
N = data(:, 3);
end

function s = get_pump_string(station_number, pump_number)
s = [num2str(station_number, '%02d') '_' num2str(pump_number)];
end

function PolyCount(Q, HN, n, Numb_Of_Fig)
for i = 0:n
    p_HNQ = polyfit(Q, HN, i);
    Q_Range = (0:1.5 * max(Q));
    f = polyval(p_HNQ, Q_Range);
    figure(Numb_Of_Fig);
    hold on;
    ylim([0 max(HN) + 50]);
    plot(Q_Range, f, '--', 'DisplayName', num2str(i));
    hold off;
    X0 = Q.^0;
    X1 = [Q.^1,Q.^0];
    X2 = [Q.^2,Q.^1,Q.^0];
    X3 = [Q.^3,Q.^2,Q.^1,Q.^0];
    eval(['X_Count = X',num2str(i),';']);
    yy0 = p_HNQ.*X_Count;
    yy1 = sum(yy0,2);
    yy1e1 = HN - yy1;
    figure(Numb_Of_Fig + 2);
    hold on;
    subplot(5,1,i+1);
    title(['poly ',num2str(i)]);
    hold off;
    plot(Q, yy1e1, 'o', 'DisplayName', 'e for poly 3');
    figure(Numb_Of_Fig + 4);
    subplot(5, 1, i+1);
    hist(yy1e1);
    title(['poly ',num2str(i)]);
end
end

function SourcePicture(Q, HN, o)
figure(o);
grid on;
hold on;
plot(Q, HN, 'r-o', 'LineWidth', 2, 'DisplayName', 'Source');
hold off;
end