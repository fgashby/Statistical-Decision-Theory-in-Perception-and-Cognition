clear all; close all; clc;  %Clear all workspaces

x0 = [.5];         % Initial estimate of probability of guessing response A
lb = [0];          % lower bound on guessing probability
ub = [1];          % upper bound on guessing probability       
[param,fval] = fmincon(@Guess,x0,[],[],[],[],lb,ub)   % call the minimization routine

% create the function that computes predictions of the guessing model
function BIC = Guess(param)
filename = 'Data15.txt';
D = importdata(filename);        % input data matrix

PC=param(1)*ones(length(D),1);
for i=1:length(D)
    if D(i,4)>1.5
        PC(i)=1-PC(i);
    end
end

lnL=sum(log(PC));              % compute log L* (Eq. 7.13)
BIC=log(length(D))-2*lnL;    % compute BIC

end