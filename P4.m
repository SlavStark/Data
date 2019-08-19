function prac4_Starkov_rev2()
clear;
clc;

%part_1
E = 50;
D = 10;
SIGMA = sqrt(D);
x=30:0.1:70;
y = normpdf(x, E, SIGMA); %theoretical distribution density of a given random variable
N = 1000;
X = randn(N, 1) * sqrt(D) + E;
figure;
hist(X); %hist of sample of random variable
title('Sampling histogram');
[p1, x1] = hist_density(X, 20); %normalized histogram with the same scale as in normpdf
figure;
plot(x, y, 'b--', x1, p1, 'r-'); %theoretical distribution density, normalized histogram
title('Theoretical distribution density, normalized histogram');
legend('Theoretical distribution density','Normalized histogram');

%part_2
n = 14 + 10;
Xi = randn(n, N) * sqrt(D) + E;
Xn = mean(Xi, 1);
y2 = normpdf(x, E, SIGMA / sqrt(n)); %theoretical distribution density for average sample
figure;
[p2,x2] = hist_density(Xn, 20); %hist_density on a sample of average
plot(x, y, 'b--', x, y2, 'g--', x2, p2, 'r-'); %theoretical distribution density, theoretical distribution density for average sample, hist_density on a sample of average 
title('Theoretical distribution density, theoretical distribution density for average sample, normalized histogram');
legend('Theoretical distribution density','Theoretical distribution density for average sample','Normalized histogram');
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