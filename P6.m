function prac6()
clear;
clc;

n = 100;
systematic1 = 60;
systematic2 = 40;
dispersion1 = 4;
dispersion2 = 2;
[P1, P2] = generate_sensor_data(n, systematic1, dispersion1, systematic2, dispersion2);
P_Dif = P1 - P2;

t = linspace(0, 240, n)';

subplot(2, 1, 1);
hold on;
ylim([0 max(P1) + 10]);
grid on;
plot(t, P1, 'r-', t, P2, 'b-');
hold off;
title('Output signal from two sensors');
subplot(2, 1, 2);
hold on;
ylim([0 max(P_Dif) + 10]);
grid on;
plot(t, P_Dif, 'g-', 'Displayname', 'Dif');
hold off;
title('Difference of two signals');

DP1 = dispersion1;
DP2 = dispersion2;
E_e_t = systematic1 - systematic2;
D_e_t = dispersion1 - dispersion2;

end

function [P1, P2] = generate_sensor_data(n, systematic1, dispersion1, systematic2, dispersion2)

t = linspace(0, 240, n)';
P_Factual_t = 100 * (1 - (exp(1).^((-1) * t./50)));
E1 = 0;
E2 = 0;
random1 = randn(n, 1) * dispersion1 + E1;
random2 = randn(n, 1) * dispersion2 + E2;
P1 = P_Factual_t + systematic1 + random1;
P2 = P_Factual_t + systematic2 + random2;

EP1 = P_Factual_t + systematic1;
EP2 = P_Factual_t + systematic2;

end