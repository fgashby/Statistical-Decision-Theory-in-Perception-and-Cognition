clear all; close all; clc;  %Clear all workspaces

x0 = [40,10];         % Starting guess at the solution
lb = [-Inf; 0];       % lower bounds for all parameters
ub = [Inf; Inf];      % upper bounds for all parameters       
[param,fval] = fmincon(@OneD,x0,[],[],[],[],lb,ub)   % call the minimization routine

% create the function that computes 1D model predictions
function BIC = OneD(param)
filename = 'Data.txt';
D = importdata(filename);        % input data matrix

%Use the following line to fit a 1D rule on dimension 1 in which the A 
%response region is on the left and the B region is on the right
Z=inv(param(2))*(param(1)*ones(length(D),1)-D(:,2));

%Use the following line to fit a 1D rule on dimension 2 in which the A
%response region is above the bound and the B region is below
%Z=inv(param(2))*(D(:,3)-param(1)*ones(length(D),1));

PC=normcdf(Z,0,1);
for i=1:length(D)
    if D(i,4)>1.5
        PC(i)=1-PC(i);
    end
end

lnL=sum(log(PC));              % compute log L* (Eq. 7.13)
BIC=2*log(length(D))-2*lnL;    % compute BIC

end