function prac9()
close all;
clear;
clc;

n=15;
Qmin=3000;
Qmax=5000;
betta_truth=[264.6178; 0; 7.9095*(10^(-5))];
k=length(betta_truth)-1;
De=10;
N=5000;

part_1(n,Qmin,Qmax,betta_truth,k,De);
part_2(n,Qmin,Qmax,betta_truth,k,De,N);

end

%part_1
function Qplan = build_plan(n,Qmin,Qmax)
Qplan = linspace(Qmin,Qmax,n)';
end

function X = build_plan_matrix(Qplan,k)
for i = 1:(k+1)
    X(:,i) = Qplan.^(i-1);
end    
end

function [Y,X,Y_truth]=generate_process_data(Qplan,betta_truth,De)
k=length(betta_truth)-1;
n=length(Qplan);
X = build_plan_matrix(Qplan,k);
for j=1:(k+1)
    for i=1:n
        Y_bord(i,j)=X(i,j).*betta_truth(j);
    end
end
Y_truth=sum(Y_bord,2);
Y=Y_truth+randn(n,1).*sqrt(De);
end

function part_1(n,Qmin,Qmax,betta_truth,k,De)

Qplan = build_plan(n,Qmin,Qmax);
[Y,X]=generate_process_data(Qplan,betta_truth,De);

betta_optimal = (((X') * X)^(-1)) * (X') * Y;
betta_truth
betta_optimal

%I-theoretical points for the interval E-theoretical points without any errors
%in measurment V-truth that we get by measuring with the mistake
QplanI = build_plan(n,0.8*Qmin,1.2*Qmax);
[YI,XI,Y_truthI]=generate_process_data(QplanI,betta_truth,De);
[YV,XV,Y_truthV]=generate_process_data(QplanI,betta_optimal,De);
figure(1);
title('part-1');
hold on;
grid on;
plot(QplanI,Y_truthI,QplanI,Y_truthV);
plot(Qplan,Y,'o');
legend('y-truth','y-recovered','dots');
hold off;
end

%part_2
function betta_optimal_MK = MonteKarloGen(N,Qplan, betta_truth, De)
for i =1:N
    [Y_MK(:,i), X_MK] = generate_process_data(Qplan, betta_truth, De);
    betta_optimal_MK_1 = (((X_MK') * X_MK)^(-1)) * (X_MK') * Y_MK(:,i);
    betta_optimal_MK(:,i)=betta_optimal_MK_1;
end
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

function part_2(n,Qmin,Qmax,betta_truth,k,De,N)
%2.1
Qplan = build_plan(n,Qmin,Qmax);
[Y,X,Y_truth]=generate_process_data(Qplan,betta_truth,De);
K=De*(((X')*X)^(-1));
DOKR=diag(K); %?????????? ?????? ????????????? ?????????

betta_optimal_MK = MonteKarloGen(N,Qplan, betta_truth, De); %variations of betta_optimal

figure(2);
num_of_betta_truth=length(betta_truth);
alpha=0.02;
for i=1:num_of_betta_truth
    [y,x]=hist_density(betta_optimal_MK(i,:),20);
    y1=normpdf(x,betta_truth(i),sqrt(DOKR(i)));
    interval=norminv([alpha,1-alpha],betta_truth(i),sqrt(DOKR(i)));
    subplot(3,1,i);
    hold on;
    grid on;
    plot(x,y,x,y1);
    plot(betta_truth(i),0,'*');
    plot(interval(1), 0,'o');
    plot(interval(2), 0,'o');
    legend('Hist-Density 2.1','Normpdf','truth','low','high');
    hold off;
end

%2.2
n22=length(Qplan);
Q22=[0.8*Qplan(1) 1.2*Qplan(n22)];
[Y_Q22,X_Q22,Y_truth_Q22]=generate_process_data(Q22,betta_truth,De);
X22=build_plan_matrix(Qplan,k);
K22=De*(((X22')*X22)^(-1));
for i=1:2
    D22(i)=X_Q22(i,:)*K22*((X_Q22(i,:))');
end

for i=1:N
[Y_N22(:,i),X_N22,Y_truth_N22(:,i)]=generate_process_data(Q22,betta_optimal_MK(:,i),De);
end

figure(3);
for i=1:2
    [y,x]=hist_density(Y_truth_N22(i,:),20);
    y1=normpdf(x,Y_truth_Q22(i),sqrt(D22(i)));
    interval=norminv([alpha,1-alpha],Y_truth_Q22(i),sqrt(D22(i)));
    subplot(2,1,i);
    hold on;
    grid on;
    plot(x,y,x,y1);
    plot(Y_truth_Q22(i),0,'*');
    plot(interval(1), 0,'o');
    plot(interval(2), 0,'o');
    legend('Hist-Density 2.2','Normpdf','truth','low','high');
    hold off;
end

%2.3
for i=1:N
    e_truth(:,i)=Y_truth_N22(:,i)-Y_truth_Q22;
end
figure(4);
for i=1:2
    [y,x]=hist_density(e_truth(i,:),20);
    y1=normpdf(x,0,sqrt(D22(i)));
    interval=norminv([alpha,1-alpha],0,sqrt(D22(i)));
    subplot(2,1,i);
    hold on;
    grid on;
    plot(x,y,x,y1);
    plot(0,0,'*');
    plot(interval(1), 0,'o');
    plot(interval(2), 0,'o');
    legend('Hist-Density 2.3','Normpdf','truth','low','high');
    hold off;
end

%2.4
for i=1:N
    e_prediction(:,i)=Y_truth_N22(:,i)-Y_Q22;
end
figure(5);
for i=1:2
    [y,x]=hist_density(e_prediction(i,:),20);
    y1=normpdf(x,0,sqrt(D22(i)+De));
    interval=norminv([alpha,1-alpha],0,sqrt(D22(i)+De));
    subplot(2,1,i);
    hold on;
    grid on;
    plot(x,y,x,y1);
    plot(0,0,'*');
    plot(interval(1), 0,'o');
    plot(interval(2), 0,'o');
    legend('Hist-Density 2.4','Normpdf','truth','low','high');
    hold off;
end

end