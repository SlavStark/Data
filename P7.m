function prac7()
clear;
clc;
close all;
N = 5000;
n_1 = 4;
Dispersion = 2;
NumFig_1 = 1;
MainPart(n_1, N, Dispersion, NumFig_1);
n_2 = 20; %task number 8
NumFig_2 = 2;
MainPart(n_2, N, Dispersion, NumFig_2);
end

function MainPart(n, N, Dispersion, NumFig)
sv = generator(n, N, Dispersion);
sv_mean = mean(sv);
S_K_O = sqrt(Dispersion/n);
z_standard = sv_mean/S_K_O;
sko_eval = std(sv);
t_student = sv_mean./(sko_eval/sqrt(n));
x1 = -3:0.01:3;
sv_teor_mean = normpdf(x1, 0, S_K_O);
[p2, x2] = hist_density(sv_mean, 20);
z_standard_teor = normpdf(x1, 0, 1);
[p4, x4] = hist_density(z_standard, 20);
t_student_teor = tpdf(x1, (n-1));
[p6, x6] = hist_density(t_student, 20);

figure(NumFig);
subplot(3, 1, 1);
plot(x1, sv_teor_mean, 'r-', x2, p2, 'b--');
grid on;
legend('normpdf', 'hist-density');
title('For Xn - sample mean');        
subplot(3, 1, 2);
plot(x1, z_standard_teor, 'r-', x4, p4, 'b--');
grid on;
legend('normpdf', 'hist-density');
title('For z - standardized average sample');        
subplot(3, 1, 3);
plot(x1, t_student_teor, 'r-', x6, p6, 'b--', x1, z_standard_teor, 'm:');
grid on;
legend('tpdf','hist-density', 'normpdf(z)');
title('For t - Student average sample');    
end

function sv = generator(n, N, Dispersion)
SIGMA = sqrt(Dispersion);
sv = SIGMA * randn(n, N);
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