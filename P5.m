function prac5()
clear;
clc;

%part_1
n = 1000;
k = 5;
x = randchi2(k, n);
[p1, x1] = hist_density(x, 20);
X = 0:0.01:30;
YYY = chi2pdf(X, k);
figure(1);
hold on;
plot(x1, p1, 'r-', X, YYY, 'b--');
hold off;
title('Range of relative frequencies, theoretical distribution density of xi^2');
legend('Range of relative frequencies','Theoretical distribution density of xi^2');

%part_2
MK_D = 10;
MK_E = 50;
N = 100000;
n = 5;
MK_Xn = randn(n, N) * sqrt(MK_D) + MK_E;

SnQuadra = std(MK_Xn, 1).^2; %biased
SnIskomoe = std(MK_Xn).^2; %non-biased

[p2, x2] = hist_density(SnQuadra, 100);
[p3, x3] = hist_density(SnIskomoe, 100);
figure(2);
hold on;
plot(x2, p2, 'r-', x3, p3, 'b--');
hold off;
title(['Range of relative frequencies of biased and non-biased estimation of variance (D), Sample n = ',num2str(n)]);
legend('Biased','Non-biased');

ZZ4 = (n - 1) * SnIskomoe/MK_D;
[p4, x4] = hist_density(ZZ4, 100);
XR = 0:0.01:30;
YYY1 = chi2pdf(XR, n-1);
figure(3);
hold on;
plot(x4, p4, 'r-', XR, YYY1, 'b--');
hold off;
title('Range of relative frequencies, theoretical distribution density');
legend('Range of relative frequencies','Theoretical distribution density');
end

function x = randchi2(k, n)
x = sum(((randn(n, k)).^2),2);
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