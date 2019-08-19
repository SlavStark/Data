function prac13()
clc;
close all;
clear all;
warning off;

%2
[Q, H, N] = read_pump_data(1, 1);

[beta_estimation_QH, aprox_QH, S2_QH, TSS_QH, ESS_QH, RSS_QH, R2_QH, R2Sovpadenie_QH, S2_Aprox_QH, TSS_Aprox_QH, ESS_Aprox_QH, RSS_Aprox_QH, R2_Aprox_QH, R2Sovpadenie_Aprox_QH, F_QH, F_Proverka_QH, significant_model_QH, p_significant_Proverka_QH, F_Aprox_QH, F_Proverka_Aprox_QH, significant_model_Aprox_QH, p_significant_Proverka_Aprox_QH, random_residuals_QH, z_QH, Ns_QH, h_test_QH,z_test_QH,Nruns_test_QH, random_residuals_Aprox_QH, z_Aprox_QH, Ns_Aprox_QH, h_test_Aprox_QH, z_test_Aprox_QH, Nruns_test_Aprox_QH, significant_coeff_matrix_QH, significant_coeff_alt_Aprox_QH] = get_beta(Q, H, 'QH');
[beta_estimation_QN, aprox_QN, S2_QN, TSS_QN, ESS_QN, RSS_QN, R2_QN, R2Sovpadenie_QN, S2_Aprox_QN, TSS_Aprox_QN, ESS_Aprox_QN, RSS_Aprox_QN, R2_Aprox_QN, R2Sovpadenie_Aprox_QN, F_QN, F_Proverka_QN, significant_model_QN, p_significant_Proverka_QN, F_Aprox_QN, F_Proverka_Aprox_QN, significant_model_Aprox_QN, p_significant_Proverka_Aprox_QN, random_residuals_QN, z_QN, Ns_QN, h_test_QN,z_test_QN,Nruns_test_QN, random_residuals_Aprox_QN, z_Aprox_QN, Ns_Aprox_QN, h_test_Aprox_QN, z_test_Aprox_QN, Nruns_test_Aprox_QN, significant_coeff_matrix_QN, significant_coeff_alt_Aprox_QN] = get_beta(Q, N, 'QN');

significant_model_QH(1)=NaN;
significant_model_QN(1)=NaN;

% Âûâîä áåç excel - ×àñòü 3
disp('×àñòü 3')
S2_QH
R2_QH
R2Sovpadenie_QH
S2_Aprox_QH
R2_Aprox_QH
R2Sovpadenie_Aprox_QH
S2_QN
R2_QN
R2Sovpadenie_QN
S2_Aprox_QN
R2_Aprox_QN
R2Sovpadenie_Aprox_QN

% Âûâîä áåç excel - ×àñòü 5
disp('×àñòü 5')
F_QH 
F_Proverka_QH 
significant_model_QH
F_QN 
F_Proverka_QN 
significant_model_QN
F_Aprox_QH 
F_Proverka_Aprox_QH 
significant_model_Aprox_QH
F_Aprox_QN 
F_Proverka_Aprox_QN 
significant_model_Aprox_QN

% Âûâîä áåç excel - ×àñòü 4
disp('×àñòü 4')
random_residuals_QH
z_QH
Ns_QH
h_test_QH
z_test_QH
Nruns_test_QH
random_residuals_QN
z_QN
Ns_QN
h_test_QN
z_test_QN
Nruns_test_QN
% äëÿ ëèòåðàòóðû (Aprox)
random_residuals_Aprox_QH
z_Aprox_QH
Ns_Aprox_QH
h_test_Aprox_QH
z_test_Aprox_QH
Nruns_test_Aprox_QH
random_residuals_Aprox_QN
z_Aprox_QN
Ns_Aprox_QN
h_test_Aprox_QN
z_test_Aprox_QN
Nruns_test_Aprox_QN

% Âûâîä áåç excel - ×àñòü 6
disp('×àñòü 6')
significant_coeff_matrix_QH
significant_coeff_matrix_QN
% H(Q) = a - b * Q^2
% N(Q) = k1 * Q - k2 * Q^2

HQSig=[significant_coeff_alt_Aprox_QH(1) 0 significant_coeff_alt_Aprox_QH(2)];
NQSig=[0 significant_coeff_alt_Aprox_QN(1) significant_coeff_alt_Aprox_QN(2)];

HQSig
NQSig

%3_a
figure;
plot (Q, H, 'o');
grid on;
title('Íàïîðíàÿ õàðàêòåðèñòèêà QH íàñîñíûõ àãðåãàòîâ');

%3_b
Q_plan=(0:10:1.5*max(Q))';

Q_plan_matrix=[Q_plan.^0 Q_plan.^1 Q_plan.^2 Q_plan.^3];
for i=1:4
H_estimate(:, i)=Q_plan_matrix(:, 1:i)*beta_estimation_QH(1:i, i);
end
hold on;
plot (Q_plan, H_estimate (:, 1), 'r', Q_plan, H_estimate (:, 2), 'g', Q_plan, H_estimate (:, 3), 'b', Q_plan, H_estimate (:, 4), 'm');
H_aprox=[Q_plan_matrix(:,1) -Q_plan_matrix(:,3)]*aprox_QH;
hold on;
plot (Q_plan, H_aprox, '--k');

%3_c
hleg = legend('Ðåàëüíûå äàííûå', 'Ïîëèíîì 0-ãî ïîðÿäêà', 'Ïîëèíîì 1-ãî ïîðÿäêà','Ïîëèíîì 2-ãî ïîðÿäêà', 'Ïîëèíîì 3-ãî ïîðÿäêà', 'Àïïðîêñèìàöèÿ èç ëèò-ðû');
ylim([0, 1.5 * max(H)]);

%4_a
figure;
plot (Q, N, 'o');
grid on;
title('Õàðàêòåðèñòèêà ìîùíîñòè QN íàñîñíûõ àãðåãàòîâ');

%4_b
for i=1:4
N_estimate(:, i)=Q_plan_matrix(:, 1:i)*beta_estimation_QN(1:i, i);
end
hold on;
plot (Q_plan, N_estimate (:, 1), 'r', Q_plan, N_estimate (:, 2), 'g', Q_plan, N_estimate (:, 3), 'b', Q_plan, N_estimate (:, 4), 'm');
N_aprox=[Q_plan_matrix(:,2) -Q_plan_matrix(:,3)]*aprox_QN;
hold on;
plot (Q_plan, N_aprox, '--k');

%4_c
ylim([0, 1.5 * max(N)]);

%4_d
hleg = legend('Ðåàëüíûå äàííûå', 'Ïîëèíîì 0-ãî ïîðÿäêà', 'Ïîëèíîì 1-ãî ïîðÿäêà','Ïîëèíîì 2-ãî ïîðÿäêà', 'Ïîëèíîì 3-ãî ïîðÿäêà', 'Àïïðîêñèìàöèÿ èç ëèò-ðû');

%5_a äëÿ QH
Q_ist=[Q.^0 Q.^1 Q.^2 Q.^3];
% èçìåíåíèå - â îòñîðòèðîâàííûé âèä ïåðåõîä ×àñòü 4
[sorted_Q_new,ind] = sort(Q); % Q çàìåíÿåì íà ind
for i=1:4
    H_estimate2(:, i)=Q_ist(:, 1:i)*beta_estimation_QH(1:i, i);
end
rest_QH=[H H H H]-H_estimate2;
figure;

subplot (5, 1, 1);
plot(ind,rest_QH (:, 1), 'o');
title('Ãðàôèê ðåãðåññèîííûõ îñòàòêîâ 0-ãî ïîðÿäêà äëÿ QH');

subplot (5, 1, 2);
plot(ind,rest_QH (:, 2), 'o');
title('Ãðàôèê ðåãðåññèîííûõ îñòàòêîâ 1-ãî ïîðÿäêà äëÿ QH');

subplot (5, 1, 3);
plot(ind,rest_QH (:, 3), 'o');
title('Ãðàôèê ðåãðåññèîííûõ îñòàòêîâ 2-ãî ïîðÿäêà äëÿ QH');

subplot (5, 1, 4);
plot(ind,rest_QH (:, 4), 'o');
title('Ãðàôèê ðåãðåññèîííûõ îñòàòêîâ 3-ãî ïîðÿäêà äëÿ QH');

H_aprox2=[Q_ist(:,1) -Q_ist(:,3)]*aprox_QH;
rest_QH2=H-H_aprox2;
subplot (5, 1, 5);
plot (ind, rest_QH2, 'o');
title('Ãðàôèê ðåãðåññèîííûõ îñòàòêîâ äëÿ QH');

%5_a äëÿ QN
for i=1:4
    N_estimate2(:, i)=Q_ist(:, 1:i)*beta_estimation_QN(1:i, i);
end
rest_QN=[N N N N]-N_estimate2;
figure;

subplot (5, 1, 1);
plot(ind,rest_QN (:, 1), 'o');
title('Ãðàôèê ðåãðåññèîííûõ îñòàòêîâ 0-ãî ïîðÿäêà äëÿ QN');

subplot (5, 1, 2);
plot(ind,rest_QN (:, 2), 'o');
title('Ãðàôèê ðåãðåññèîííûõ îñòàòêîâ 1-ãî ïîðÿäêà äëÿ QN');

subplot (5, 1, 3);
plot(ind,rest_QN (:, 3), 'o');
title('Ãðàôèê ðåãðåññèîííûõ îñòàòêîâ 2-ãî ïîðÿäêà äëÿ QN');

subplot (5, 1, 4);
plot(ind,rest_QN (:, 4), 'o');
title('Ãðàôèê ðåãðåññèîííûõ îñòàòêîâ 3-ãî ïîðÿäêà äëÿ QN');

N_aprox2=[Q_ist(:,2) -Q_ist(:,3)]*aprox_QN;
rest_QN2=N-N_aprox2;
subplot (5, 1, 5);
plot (ind, rest_QN2, 'o');
title('Ãðàôèê ðåãðåññèîííûõ îñòàòêîâ äëÿ QN');

%5_b äëÿ QH
figure;
subplot (5, 1, 1);
hist(rest_QH (:, 1), 20);
title('Ãèñòîãðàììà ðåãðåññèîííûõ îñòàòêîâ 0-ãî ïîðÿäêà äëÿ QH');

subplot (5, 1, 2);
hist(rest_QH (:, 2), 20);
title('Ãèñòîãðàììà ðåãðåññèîííûõ îñòàòêîâ 1-ãî ïîðÿäêà äëÿ QH');

subplot (5, 1, 3);
hist(rest_QH (:, 3), 20);
title('Ãèñòîãðàììà ðåãðåññèîííûõ îñòàòêîâ 2-ãî ïîðÿäêà äëÿ QH');

subplot (5, 1, 4);
hist(rest_QH (:, 4), 20);
title('Ãèñòîãðàììà ðåãðåññèîííûõ îñòàòêîâ 3-ãî ïîðÿäêà äëÿ QH');

subplot (5, 1, 5);
hist(rest_QH2, 20);
title('Ãèñòîãðàììà ðåãðåññèîííûõ îñòàòêîâ äëÿ QH');

%5_b äëÿ QN
figure;
subplot (5, 1, 1);
hist(rest_QN (:, 1), 20);
title('Ãèñòîãðàììà ðåãðåññèîííûõ îñòàòêîâ 0-ãî ïîðÿäêà äëÿ QN');

subplot (5, 1, 2);
hist(rest_QN (:, 2), 20);
title('Ãèñòîãðàììà ðåãðåññèîííûõ îñòàòêîâ 1-ãî ïîðÿäêà äëÿ QN');

subplot (5, 1, 3);
hist(rest_QN (:, 3), 20);
title('Ãèñòîãðàììà ðåãðåññèîííûõ îñòàòêîâ 2-ãî ïîðÿäêà äëÿ QN');

subplot (5, 1, 4);
hist(rest_QN (:, 4), 20);
title('Ãèñòîãðàììà ðåãðåññèîííûõ îñòàòêîâ 3-ãî ïîðÿäêà äëÿ QN');

subplot (5, 1, 5);
hist(rest_QN2, 20);
title('Ãèñòîãðàììà ðåãðåññèîííûõ îñòàòêîâ äëÿ QN');

% Âûâîä â excel:
Label={'Ïîëèíîì 0','Ïîëèíîì 1','Ïîëèíîì 2','Ïîëèíîì 3','Ëèòåðàòóðà'}';
Names={'Íîìåð ìîäåëè','Betta 0','Betta 1','Betta 2','Betta 3','R2','R2 Ïðîâåðêà','Àíàëèç îñòàòêîâ','Îöåíêà äèñïåðñèè øóìà','Çíà÷èìîñòü ìîäåëè','F-Ñòàòèñòèêà','F-Ïðîâåðêà','Çíà÷èìîñòü êîýôôèöèåíòîâ'};

filename='models_Berdalieva.xls';
xlswrite(filename,Names,1,'A1')
xlswrite(filename,Label,1,'A2')
xlswrite(filename,Label,1,'A8')

% ×àñòü 3
xlswrite(filename,R2_QH,1,'F2')
xlswrite(filename,R2Sovpadenie_QH,1,'G2')
xlswrite(filename,S2_QH,1,'I2')
xlswrite(filename,R2_Aprox_QH,1,'F6')
xlswrite(filename,R2Sovpadenie_Aprox_QH,1,'G6')
xlswrite(filename,S2_Aprox_QH,1,'I6')

xlswrite(filename,R2_QN,1,'F8')
xlswrite(filename,R2Sovpadenie_QN,1,'G8')
xlswrite(filename,S2_QN,1,'I8')
xlswrite(filename,R2_Aprox_QN,1,'F12')
xlswrite(filename,R2Sovpadenie_Aprox_QN,1,'G12')
xlswrite(filename,S2_Aprox_QN,1,'I12')

% Âûâîä áýòò
BettaFirst_H=[beta_estimation_QH(1, :) aprox_QH(1)]';
BettaSecond_H=[beta_estimation_QH(2, :) 0]';
BettaThird_H=[beta_estimation_QH(3, :) aprox_QH(2)]';
BettaFourth_H=[beta_estimation_QH(4, :) 0]';
BettaFirst_N=[beta_estimation_QN(1, :) 0]';
BettaSecond_N=[beta_estimation_QN(2, :) aprox_QN(1)]';
BettaThird_N=[beta_estimation_QN(3, :) aprox_QN(2)]';
BettaFourth_N=[beta_estimation_QN(4, :) 0]';
% H(Q) = a - b * Q^2
% N(Q) = k1 * Q - k2 * Q^2
xlswrite(filename,BettaFirst_H,1,'B2')
xlswrite(filename,BettaSecond_H,1,'C2')
xlswrite(filename,BettaThird_H,1,'D2')
xlswrite(filename,BettaFourth_H,1,'E2')
xlswrite(filename,BettaFirst_N,1,'B8')
xlswrite(filename,BettaSecond_N,1,'C8')
xlswrite(filename,BettaThird_N,1,'D8')
xlswrite(filename,BettaFourth_N,1,'E8')

% ×àñòü 5
xlswrite(filename,significant_model_QH,1,'J2')
xlswrite(filename,significant_model_QN,1,'J8')
xlswrite(filename,F_QH,1,'K2')
xlswrite(filename,F_QN,1,'K8')
xlswrite(filename,F_Proverka_QH,1,'L2')
xlswrite(filename,F_Proverka_QN,1,'L8')
xlswrite(filename,significant_model_Aprox_QH,1,'J6')
xlswrite(filename,significant_model_Aprox_QN,1,'J12')
xlswrite(filename,F_Aprox_QH,1,'K6')
xlswrite(filename,F_Aprox_QN,1,'K12')
xlswrite(filename,F_Proverka_Aprox_QH,1,'L6')
xlswrite(filename,F_Proverka_Aprox_QN,1,'L12')

% ×àñòü 4
xlswrite(filename,random_residuals_QH,1,'H2')
xlswrite(filename,random_residuals_QN,1,'H8')
xlswrite(filename,random_residuals_Aprox_QH,1,'H6')
xlswrite(filename,random_residuals_Aprox_QN,1,'H12')
% âûâîä äëÿ ïðîâåðêè
LabelPart4={'Ïðîâåðêà ñëó÷àéíîñòè','z-ñòàòèñòèêà','Êîë ñåðèé','Ïðîâåðêà ïðîâåðêè ñëó÷àéíîñòè','z-ïðîâåðêà','Êîë ñåðèé ïðîâåðêà'};
LabelPart4Title={'Ïðîâåðêà ×àñòü 4'};
xlswrite(filename,LabelPart4,1,'H14')
xlswrite(filename,LabelPart4Title,1,'G14')
xlswrite(filename,random_residuals_QH,1,'H15')
xlswrite(filename,random_residuals_QN,1,'H21')
xlswrite(filename,random_residuals_Aprox_QH,1,'H19')
xlswrite(filename,random_residuals_Aprox_QN,1,'H25')
xlswrite(filename,z_QH,1,'I15')
xlswrite(filename,z_QN,1,'I21')
xlswrite(filename,z_Aprox_QH,1,'I19')
xlswrite(filename,z_Aprox_QN,1,'I25')
xlswrite(filename,Ns_QH,1,'J15')
xlswrite(filename,Ns_QN,1,'J21')
xlswrite(filename,Ns_Aprox_QH,1,'J19')
xlswrite(filename,Ns_Aprox_QN,1,'J25')
% test
xlswrite(filename,h_test_QH,1,'K15')
xlswrite(filename,h_test_QN,1,'K21')
xlswrite(filename,h_test_Aprox_QH,1,'K19')
xlswrite(filename,h_test_Aprox_QN,1,'K25')
xlswrite(filename,z_test_QH,1,'L15')
xlswrite(filename,z_test_QN,1,'L21')
xlswrite(filename,z_test_Aprox_QH,1,'L19')
xlswrite(filename,z_test_Aprox_QN,1,'L25')
xlswrite(filename,Nruns_test_QH,1,'M15')
xlswrite(filename,Nruns_test_QN,1,'M21')
xlswrite(filename,Nruns_test_Aprox_QH,1,'M19')
xlswrite(filename,Nruns_test_Aprox_QN,1,'M25')

% ×àñòü 6
xlswrite(filename,significant_coeff_matrix_QH,1,'M2')
xlswrite(filename,significant_coeff_matrix_QN,1,'M8')
xlswrite(filename,HQSig,1,'M6')
xlswrite(filename,NQSig,1,'M12')

end

function [beta_estimation, aprox, S2, TSS, ESS, RSS, R2, R2Sovpadenie, S2_Aprox, TSS_Aprox, ESS_Aprox, RSS_Aprox, R2_Aprox, R2Sovpadenie_Aprox, F, F_Proverka, significant_model, p_significant_Proverka, F_Aprox, F_Proverka_Aprox, significant_model_Aprox, p_significant_Proverka_Aprox, random_residuals, z, Ns, h_test, z_test, Nruns_test, random_residuals_Aprox, z_Aprox, Ns_Aprox, h_test_Aprox, z_test_Aprox, Nruns_test_Aprox, significant_coeff_matrix, significant_coeff_alt_Aprox] = get_beta(x, y, param)

X = [x.^0 x.^1 x.^2 x.^3];  
for i = 1 : 4
    beta_estimation(1:i, i) = (X(:, 1:i)'*X(:, 1:i))^-1*X(:, 1:i)'*y;
    B_reg(1 : i, i) = regress(y, X(:, 1 : i));
    
    Y_New(:, i) = X(:, 1:i) * beta_estimation(1:i, i);
    % ×àñòü 3
    S2(i, 1) = calc_residual_dispersion(y, Y_New(:, i), i);
    TSS(i, 1) = calc_TSS(y);
    ESS(i, 1) = calc_ESS(y, Y_New(:, i));
    RSS(i, 1) = calc_RSS(y, Y_New(:, i));
    R2(i, 1) = calc_R2(y, Y_New(:, i));
    R2Sovpadenie(i, 1) = RegressProverkaNaSovpadenie(y, X(:, 1:i));
    
    % ×àñòü 5
    [F(i, 1),n1,n2]=model_Fstat(y,Y_New(:, i),i);
    significant_model(i, 1)=test_model_significance(F(i, 1),n1,n2,0.05);
    [~,~,~,~,stats5] = regress(y,X(:, 1:i));
    F_Proverka(i, 1)=stats5(2);
    p_significant_Proverka(i, 1)=stats5(3);
    
    % ×àñòü 4
    E(:, i) = y - Y_New(:, i);
    [random_residuals(i, 1),z(i, 1),Ns(i, 1),h_test(i, 1),z_test(i, 1),Nruns_test(i, 1)]=test_residuals_randomness(E(:, i));
    
    % ×àñòü 6
    K=S2(i, 1)*(X(:, 1:i)'*X(:, 1:i))^(-1);
    K_est=diag(K);
    n_len = length(x);
    [beta_st,beta_bounds]=studentize_beta_est(beta_estimation(1:i, i),K_est,n_len,0.05,i);
    significant_coeff=test_coefficient_significance(beta_st, n_len, 0.05,i);
    [~,bint,~,~,~]=regress(y,X(:, 1:i),0.05);
    beta_bounds=bint;
    significant_coeff_alt=test_coefficient_significance_alt(beta_bounds,i);
    r=i;
    if r ==1
        A=significant_coeff_alt;
    elseif r==2
        B=significant_coeff_alt;
    elseif r==3
        C=significant_coeff_alt;
    elseif r==4
        D=significant_coeff_alt;
    end
    
end

significant_coeff_matrix=[A NaN NaN NaN;B NaN NaN;C NaN;D];

if param == 'QH'
x_aprox = [X(:, 1) -X(:,3)];
end
if param == 'QN'
x_aprox = [X(:, 2) -X(:,3)];
end
aprox(:, 1) = inv(x_aprox' * x_aprox) * x_aprox' * y;

Y_New_Aprox = x_aprox * aprox(:, 1);
% ×àñòü 3
S2_Aprox = calc_residual_dispersion(y, Y_New_Aprox, 3);
TSS_Aprox = calc_TSS(y);
ESS_Aprox = calc_ESS(y, Y_New_Aprox);
RSS_Aprox = calc_RSS(y, Y_New_Aprox);
R2_Aprox = calc_R2(y, Y_New_Aprox);
R2Sovpadenie_Aprox = RegressProverkaNaSovpadenie(y, x_aprox);
% ×àñòü 5
[F_Aprox,n1_Aprox,n2_Aprox]=model_Fstat(y,Y_New_Aprox,2);
significant_model_Aprox=test_model_significance(F_Aprox,n1_Aprox,n2_Aprox,0.05);
[~,~,~,~,stats5_Aprox] = regress(y,x_aprox);
F_Proverka_Aprox=stats5_Aprox(2);
p_significant_Proverka_Aprox=stats5_Aprox(3);
% ×àñòü 4
E_Aprox = y - Y_New_Aprox;
[random_residuals_Aprox,z_Aprox,Ns_Aprox,h_test_Aprox,z_test_Aprox,Nruns_test_Aprox]=test_residuals_randomness(E_Aprox);
% ×àñòü 6
K_Aprox=S2_Aprox*(x_aprox'*x_aprox)^(-1);
K_est_Aprox=diag(K_Aprox);
n_len_Aprox = length(x);
[beta_st_Aprox,beta_bounds_Aprox]=studentize_beta_est(aprox(:, 1),K_est_Aprox,n_len_Aprox,0.05,2);
significant_coeff_Aprox=test_coefficient_significance(beta_st_Aprox, n_len_Aprox, 0.05,2);
[~,bint_Aprox,~,~,~]=regress(y,x_aprox,0.05);
beta_bounds_Aprox=bint_Aprox;
significant_coeff_alt_Aprox=test_coefficient_significance_alt(beta_bounds_Aprox,2);

end

% Ñ÷èòûâàíèå äàííûõ
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

% ×ÀÑÒÜ ¹3
function S2 = calc_residual_dispersion(y, y_prediction, k)
  n = length(y);
  S2 = sum((y - y_prediction).^2) * (1 / (n - k));
end

function TSS = calc_TSS(y)
  y_sred = mean(y);
  TSS = sum((y - y_sred).^2);
end

function ESS = calc_ESS(y, y_prediction)
  ESS = sum((y_prediction - mean(y)).^2);
end

function RSS = calc_RSS(y, y_prediction)
  RSS = sum((y - y_prediction).^2);
end

function R2 = calc_R2(y, y_prediction)
  RSS = calc_RSS(y, y_prediction);
  TSS = calc_TSS(y);
  R2 = 1 - (RSS / TSS);
end

function [R2Sovpadenie] = RegressProverkaNaSovpadenie(y, X)
  [~, ~, ~, ~, stats] = regress(y, X);
  % Âûíîñèì R^2
  R2Sovpadenie = stats(1);
end

% ×àñòü 4
function Ns=get_series_count(e)
Ns=0;
for i=1:(length(e)-1)
    if ((e(i)>0 && e(i+1)<0) || (e(i)<0 && e(i+1)>0))
        Ns = Ns+1;
    end
end
end

function [Nplus,Nminus]=get_signed_values_count(e)
Nplus=0;
for i=1:length(e)
    if e(i,1)>=0
        Nplus=Nplus+1;
    end
end
Nminus=length(e)-Nplus;
end

function [random_residuals,z,Ns,h_test,z_test,Nruns_test]=test_residuals_randomness(e)
[Nplus,Nminus]=get_signed_values_count(e);
Ns=get_series_count(e);
ENs=((2*Nplus*Nminus)/(Nplus+Nminus))+1;
DNs=((2*Nplus*Nminus)/((Nplus+Nminus)^2))*((2*Nplus*Nminus-(Nplus+Nminus))/(Nplus+Nminus-1));
z=(Ns-ENs)/(sqrt(DNs));
[~,~,stats]=runstest(e);
% disp(stats.nruns);
% disp(z);
% disp(stats.z);
pdf = normpdf(z,0,1);
random_residuals=1;
if (pdf > 0.05/2)
    random_residuals=0;   
end

% ïðîâåðêà 
[h_test,~,stats_test]=runstest(e);
Nruns_test = stats_test.nruns;
z_test = stats_test.z;

end

% ×àñòü 5

function [F,n1,n2]=model_Fstat(y,y_prediction,k)
n1=k-1;
n2=length(y)-k;
ESS=calc_ESS(y,y_prediction);
RSS=calc_RSS(y,y_prediction);
F=(ESS/n1)/(RSS/n2);
end

function significant_model=test_model_significance(F,n1,n2,alpha)
X=finv(1-alpha,n1,n2);
significant_model=1; % 1 - çíà÷èò çíà÷èìûé
if F < X
    significant_model=0; % 0 - çíà÷èò íåçíà÷èìûé
end
end

% ×àñòü 6

function [beta_st, beta_bounds]=studentize_beta_est(beta_est,K_est,n,alpha,l)

beta_st=beta_est./(sqrt(K_est));
t_porog=tinv(1-alpha,(n-l-1));
beta_bound_1=beta_est-t_porog;
beta_bound_2=beta_est+t_porog;
beta_bounds=[beta_bound_1 beta_bound_2];

end

function significant_coeff=test_coefficient_significance(beta_st,n,alpha,l)

border_low(l) = tinv(alpha,(n-l-1));
border_high(l) = tinv(1-alpha,(n-l-1));
significant_coeff(l)=1;
if beta_st(l) > border_low(l) && beta_st(l) < border_high(l)
    significant_coeff(l)=0;
end

end

function significant_coeff_alt=test_coefficient_significance_alt(beta_bounds,l)

for i=1:l
    significant_coeff_alt(i)=1;
    if (beta_bounds(i,1)<0) && (beta_bounds(i,2)>0) 
       significant_coeff_alt(i)=0;
    end
end

end
