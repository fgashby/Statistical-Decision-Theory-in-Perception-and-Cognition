clear all; close all; clc;      % clear all workspaces

x0 = [1; -1; 0; 10];        % starting guesses for all parameters, which
                                % assume the A response region is above
                                % the bound and the B region is below
lb = [-Inf; -Inf; -Inf; 0];     % lower bounds for all parameters
ub = [Inf; Inf; Inf; Inf];      % upper bounds for all parameters       
[param,fval] = fmincon(@GLC,x0,[],[],[],[],lb,ub)  % call the minimization routine

% create the function that computes GLC predictions
function BIC = GLC(param)
filename = 'Data.txt'; 
D = importdata(filename);       % input data matrix

Z=inv(param(4))*(param(1)*D(:,2)+param(2)*D(:,3)+param(3)*ones(length(D),1));
PC=normcdf(Z,0,1);
for i=1:length(D)
    if D(i,4)>1.5
        PC(i)=1-PC(i);
    end
end

lnL=sum(log(PC));               % compute log L* (Eq. 7.13)
BIC=3*log(length(D))-2*lnL;     % compute BIC

end