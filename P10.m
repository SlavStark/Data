function prac10()
close all;
clear;
clc;

[u, y] = load_tf_ident_data(6);
n1=length(u);
n2=length(y);
%t1u=linspace(1,100,n1);
%t2y=linspace(1,100,n2);
t1u=1:length(u);
t2y=1:length(y);
figure(1);
subplot(1, 2, 1)
plot(t1u,u);
legend('managing impact u');
subplot(1, 2, 2)
plot(t2y,y);
legend('feedback y');

%for 2 order
n=2;
[TransFunc, M, Y, betta] = ident_tf(u, y, n);
figure(2);
subplot(3,1,1);
hold on;
title('Actual control action u(t)');
plot(t1u, u);
hold off;
subplot(3,1,2);
hold on;
title('The actual response y(t),response of the identified system Y(t)');
du=u-u(1);
tt=1:length(y);
dy=lsim(TransFunc,du,tt);
plot(tt,y,tt,dy+y(1));
legend('y(t)','Y(t)');
hold off;
subplot(3,1,3);
hold on;
title('Remains on response e(t)');
e=y-dy;
plot(t2y,e);
hold off;

figure(3);
subplot(3,1,1);
hold on;
title('Actual control action u(t)');
plot(t1u,u);
hold off;
subplot(3,1,2);
hold on;
title('Actual leading derivative (y^n)(t), forecast (Y^n)(t)');
Yn=M*betta;
t1u_n=linspace(1,100,n1-n);
plot(t1u_n, Y,t1u_n,Yn);
legend('(y^n)(t)','(Y^n)(t)');
hold off;
subplot(3,1,3);
hold on;
title('Residuals by derivative e(t)');
e=Y-Yn;
plot(t1u_n,e);
hold off;

end

%part_1
function [u, y] = load_tf_ident_data(var_number)

fname = ['tf_ident_data/' num2str(var_number) '.csv'];
data = load(fname);

u = data(2:end, 1);
y = data(2:end, 2);

end

%part_2 - Identification
%plan matrix
function [M, MY] = build_plan_matrix(u, y, n)

vector(:,1)= y(1:length(y));
k=n+1;
MY(:,k) = vector(1:length(y)-n,1);
for i=1:n
    vector(1:length(y)-i,i+1) = diff(vector(:,i),i);
    k=k-1;   
    MY(:,k) = vector(1:length(y)-n,i+1);
end
M = [- MY(1:length(y)-n,2:end),u(1:length(y)-n)];

end

function yd = diff(y, j)
yd = zeros(length(y)-j, 1);
 for i=1:length(y)-j
     yd(i)=(y(i+1)-y(i))/1;
 end
end

%identification of transfer function
function [TransFunc, M, Y, betta] = ident_tf(u, y, n)

u=u-u(1);
y=y-y(1);
[M, MY] = build_plan_matrix(u, y, n);
Y=MY(:,1);
betta=(((M')*M)^(-1))*(M')*Y;
TransFunc=tf(betta(3),[1 betta(1) betta(2)]);

end
