function prac12()
close all;
clear;
clc;
warning('off');

station_number = 1;
pump_number = 1;
[Q_data, H_data, N_data] = read_pump_data(station_number, pump_number);
[betta_data,X_data] = approx(Q_data, H_data);
H_estimated = X_data*betta_data;

e_of_model = H_data - H_estimated; % calculating e remains
n = length(H_data);
D_of_e = sum(e_of_model.^2)/n; % formula for dispersion of remains
D_of_e
k = 3; % difficulty of the model - connected with the number of koeff - a,b,c..
S_for_e = sum(e_of_model.^2)/(n - k);
S_for_e

% leave_one_out = LOO - kicking out one point of data and finding the dispersion of noise
counter = 0;
M_data=X_data;
H_loo=H_data;
for i = 1:n
    M_data(i,:)=[];
    H_loo(i,:)=[];
    % calculating new betta for new data without i-element
    betta_loo = (((M_data') * M_data)^(-1)) * (M_data') * H_loo;
    e_loo(i)=H_data(i)-X_data(i,:)*betta_loo;
    M_data=X_data;
    H_loo=H_data;
end
D_e_loo=sum(e_loo.^2)/n;
D_e_loo

for i=1:n
    D_of_e_n=sqrt(S_for_e*(X_data(i,:)*((X_data'*X_data)^(-1))*X_data(i,:)'+1));
    t_stat=e_loo(i)/D_of_e_n;
    bounds = tinv([0.05/2, 1-0.05/2], n-3);
    if (t_stat>=bounds(1) && t_stat<=bounds(2)) 
        counter=counter+1;
    end
end

% number of hits
counter
% percentage of hits
percentage = (counter/n) * 100;
percentage

end

function [Q, H, N] = read_pump_data(station_number, pump_number)
filename = ['.\123\' get_pump_string(station_number, pump_number) '.csv'];
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

function [betta_data,X] = approx(Q_data, H_data)

for i=0:2
    X(:,i+1)=Q_data.^(i);
end
betta_data = (((X') * X)^(-1)) * (X') * H_data; % optimal

end
