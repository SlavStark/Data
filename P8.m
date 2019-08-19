function prac8()
close all;
clear;
clc;

n = 15;
N = 10000;
systematic1 = 60;
systematic2 = 58;
dispersion1 = 4;
dispersion2 = 0;
[P1, P2] = generate_sensor_data(n, N, systematic1, dispersion1, systematic2, dispersion2); %hypothesis1
[P1_0, P2_0] = generate_sensor_data(n, N, 0, dispersion1, 0, dispersion2); %hypothesis0

%part_3
%e statistics
t_x = -2:0.1:4;
SKO_1 = sqrt((dispersion1 + dispersion2)/n);
e_H0 = e_statistic(P1_0, P2_0);
e_H1 = e_statistic(P1, P2);
g1 = normpdf(t_x, 0, SKO_1); %theoretical H0; E and D: (0) (D1+D2)/n
g2 = normpdf(t_x, 2, SKO_1); %theoretical H1; E and D:(systematic1-systematic2) (D1+D2)/n
[y1, x1] = hist_density(e_H0,20); %hist_density H0
[y2, x2] = hist_density(e_H1,20); %hist_density H1
[h_e_H0, border_low, border_high] = test_equal_means_e(e_H0, dispersion1, dispersion2, n, 0.02); %border
figure(1);
subplot(3,1,1);
title('e statistics');
hold on;
grid on;
plot(t_x, g1, t_x, g2, x1, y1, x2, y2);
plot(border_low, 0,'o');
plot(border_high, 0, 'o');
legend('Theoretical e, H0','Theoretical e, H1','HistDensity e, H0','HistDensity e, H1','low','high');
hold off;

%z statistics
t_x = -4:0.1:8;
z_H0 = z_statistic(P1_0, P2_0, dispersion1, dispersion2);
z_H1 = z_statistic(P1, P2, dispersion1, dispersion2);
g3 = normpdf(t_x, 0, 1); %theoretical H0; E and D: 0 1
E_z = (systematic1-systematic2)/sqrt((dispersion1 + dispersion2)/n);
g4 = normpdf(t_x, E_z, 1); %theoretical H1; E and D: (systematic1-systematic2)/sqrt((D1+D2)/n) 1
[y3, x3] = hist_density(z_H0,20); %hist_density H0
[y4, x4] = hist_density(z_H1,20); %hist_density H1
[h_z_H0, border_low, border_high] = test_equal_means_z(z_H0, 0.02); %border
figure(1);
subplot(3,1,2);
title('z statistics');
hold on;
grid on;
plot(t_x, g3, t_x, g4, x3, y3, x4, y4);
plot(border_low, 0,'o');
plot(border_high, 0, 'o');
legend('Theoretical z, H0','Theoretical z, H1','HistDensity z, H0','HistDensity z, H1','low','high');
hold off;

%t statistics
t_x = -4:0.1:8;
t_H0 = t_statistic(P1_0, P2_0);
t_H1 = t_statistic(P1, P2);
g5 = tpdf(t_x, n-1);
[y5, x5] = hist_density(t_H0,20); %hist_density H0
[y6, x6] = hist_density(t_H1,20); %hist_density H1
[h_t_H0, border_low, border_high] = test_equal_means_t(t_H0, n, 0.02); %border
figure(1);
subplot(3,1,3);
title('t statistics');
hold on;
grid on;
plot(t_x, g5, x5, y5, x6, y6);
plot(border_low, 0,'o');
plot(border_high, 0, 'o');
legend('Theoretical t, H0','HistDensity t, H0','HistDensity t, H1','low','high');
hold off;

%part_4
%n = 15;
%fot type I errors = a errors = false positives
error_e_1 = 1 - h_e_H0
error_z_1 = 1 - h_z_H0
%for type II errors = b errors = false negatives
[h_e_H1, border_low_e_H1, border_high_e_H1] = test_equal_means_e(e_H1, dispersion1, dispersion2, n, 0.02);
[h_z_H1, border_low_z_H1, border_high_z_H1] = test_equal_means_z(z_H1, 0.02);
error_e_2 = h_e_H1
error_z_2 = h_z_H1

%n*2 = 30;
n2 = 30;
%fot type I errors
[P1_n2, P2_n2] = generate_sensor_data(n2, N, systematic1, dispersion1, systematic2, dispersion2); %hypothesis1_n2
[P1_0_n2, P2_0_n2] = generate_sensor_data(n2, N, 0, dispersion1, 0, dispersion2); %hypothesis0_n2
e_H0_n2 = e_statistic(P1_0_n2, P2_0_n2);
e_H1_n2 = e_statistic(P1_n2, P2_n2);
z_H0_n2 = z_statistic(P1_0_n2, P2_0_n2, dispersion1, dispersion2);
z_H1_n2 = z_statistic(P1_n2, P2_n2, dispersion1, dispersion2);

[h_e_H0_n2, border_low_e_H0_n2, border_high_e_H0_n2] = test_equal_means_e(e_H0_n2, dispersion1, dispersion2, n2, 0.02);
[h_z_H0_n2, border_low_z_H0_n2, border_high_z_H0_n2] = test_equal_means_z(z_H0_n2, 0.02);
[h_e_H1_n2, border_low_e_H1_n2, border_high_e_H1_n2] = test_equal_means_e(e_H1_n2, dispersion1, dispersion2, n2, 0.02);
[h_z_H1_n2, border_low_z_H1_n2, border_high_z_H1_n2] = test_equal_means_z(z_H1_n2, 0.02);

error_e_1_n2 = 1 - h_e_H0_n2
error_z_1_n2 = 1 - h_z_H0_n2
%for type II errors
error_e_2_n2 = h_e_H1_n2
error_z_2_n2 = h_z_H1_n2

%part_5
part5();

end

function [P1, P2] = generate_sensor_data(n, N, systematic1, dispersion1, systematic2, dispersion2)

t = linspace(0, 240, n)';
P_Factual_t = 3 + 100 * (1 - (exp(1).^((-1) * t./50)));
P_Factual_t = repmat(P_Factual_t,1,N);
E1 = 0;
E2 = 0;
random1 = randn(n, N) * sqrt(dispersion1) + E1;
random2 = randn(n, N) * sqrt(dispersion2) + E2;
P1 = P_Factual_t + systematic1 + random1;
P2 = P_Factual_t + systematic2 + random2;

end

function e = e_statistic(P1, P2)
e1 = P1 - P2;
e = mean(e1);
end

function z = z_statistic(P1, P2, dispersion1, dispersion2)
e = e_statistic(P1, P2);
n = size(P1,1);
z = e./(sqrt((dispersion1 + dispersion2)/n));
end

function [h, border_low, border_high] = test_equal_means_e(e, dispersion1, dispersion2, n, alpha)
SKO = sqrt((dispersion1 + dispersion2)/n);
border = norminv([alpha, 1-alpha], 0, SKO);
border_low = border(1);
border_high = border(2);
f = find(e > border_low &  e < border_high);
l1 = length(f);
l2 = length(e);
h = l1/l2; %probability to exceed the threshold value
end

function [h, border_low, border_high] = test_equal_means_z(z, alpha)
border = norminv([alpha, 1-alpha], 0, 1);
border_low = border(1);
border_high = border(2);
f = find(z > border_low &  z < border_high);
l1 = length(f);
l2 = length(z);
h = l1/l2;
end

function t = t_statistic(P1, P2)
e1 = P1 - P2;
e = mean(e1);
Se = std(e1,0);
n = size(P2,1);
t = e./(sqrt((Se.^2)/n));
end

function [h, border_low, border_high] = test_equal_means_t(t, n, alpha)
border = tinv([alpha, 1-alpha], n-1);
border_low = border(1);
border_high = border(2);
f = find(t > border_low &  t < border_high);
l1 = length(f);
l2 = length(t);
h = l1/l2;
end


function varargout = hist_density(x, bin_count)
n1 = length(x); min_x = min(x); max_x = max(x);
dx = (max_x - min_x) / bin_count;
[counts, centers] = hist(x, bin_count);
density = (counts/n1)/dx;
if (nargout == 0)
    plot(centers, density);
else if (nargout == 2)
        varargout{1} = density;
        varargout{2} = centers;
    end
end
end

function part5()
%for n diversification
N = 10000;
systematic0 = 0;
systematic1 = 60;
systematic2 = 58;
dispersion1 = 4;
dispersion2 = 0;
alpha = 0.02;

for n = 6:60
    [P1_0, P2_0] = generate_sensor_data(n, N, systematic0, dispersion1, systematic0, dispersion2);
    [P1, P2] = generate_sensor_data(n, N, systematic1, dispersion1, systematic2, dispersion2);
    z_H0 = z_statistic(P1_0, P2_0, dispersion1, dispersion2);
    z_H1 = z_statistic(P1, P2, dispersion1, dispersion2);
    t_H0 = t_statistic(P1_0, P2_0);
    t_H1 = t_statistic(P1, P2);
    [h_z_H0(1,n-5)] = test_equal_means_z(z_H0, alpha);
    [h_z_H1(1,n-5)] = test_equal_means_z(z_H1, alpha);
    [h_t_H0(1,n-5)] = test_equal_means_t(t_H0, n, alpha);
    [h_t_H1(1,n-5)] = test_equal_means_t(t_H1, n, alpha);
end

error_z_1 = 1 - h_z_H0;
error_t_1 = 1 - h_t_H0;
error_z_2 = h_z_H1;
error_t_2 = h_t_H1;
figure(2);
title('1,2 type errors for z,t by diversifying n');
hold on;
grid on;
n1=6:1:60;
plot(n1,error_z_1,n1,error_t_1,n1,error_z_2,n1,error_t_2);
legend('error-z-1','error-t-1','error-z-2','error-t-2');
hold off;

%for systematic diversification
%systematic1: 58--70
n = 15;
for systematic1 = 58:70
    [P1_0_s, P2_0_s] = generate_sensor_data(n, N, systematic0, dispersion1, systematic0, dispersion2);
    [P1_s, P2_s] = generate_sensor_data(n, N, systematic1, dispersion1, systematic2, dispersion2);
    z_H0_s = z_statistic(P1_0_s, P2_0_s, dispersion1, dispersion2);
    z_H1_s = z_statistic(P1_s, P2_s, dispersion1, dispersion2);
    t_H0_s = t_statistic(P1_0_s, P2_0_s);
    t_H1_s = t_statistic(P1_s, P2_s);
    [h_z_H0_s(1,systematic1-57)] = test_equal_means_z(z_H0_s, alpha);
    [h_z_H1_s(1,systematic1-57)] = test_equal_means_z(z_H1_s, alpha);
    [h_t_H0_s(1,systematic1-57)] = test_equal_means_t(t_H0_s, n, alpha);
    [h_t_H1_s(1,systematic1-57)] = test_equal_means_t(t_H1_s, n, alpha);
end

error_z_1_s = 1 - h_z_H0_s;
error_t_1_s = 1 - h_t_H0_s;
error_z_2_s = h_z_H1_s;
error_t_2_s = h_t_H1_s;
figure(3);
title('1,2 type errors for z,t by diversifying systematic1');
hold on;
grid on;
s1=58:1:70;
plot(s1,error_z_1_s,s1,error_t_1_s,s1,error_z_2_s,s1,error_t_2_s);
legend('error-z-1-s','error-t-1-s','error-z-2-s','error-t-2-s');
hold off;

end