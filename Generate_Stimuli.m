clear all; close all; clc;  %Clear all workspaces

N = 330;    % Number of samples to draw for each category
N_final = 300;  % Final number of stimuli in each category
mu_RB_A = [50; 65];   % Mean of category A in the RB condition
mu_RB_B = [50; 35];   % Mean of category B in the RB condition
mu_II_A = [39.4; 60.6];   % Mean of category A in the II condition
mu_II_B = [60.6; 39.4];   % Mean of category B in the II condition
Sig_RB = [225 0; 0 41.6];    % Variance-covariance matrix in the RB condition
Sig_II = [133.3 91.7; 91.7 133.3];    % Variance-covariance matrix in the II condition
RB_A = mvnrnd(mu_RB_A,Sig_RB,N);    % Category A raw random samples in RB condition
RB_B = mvnrnd(mu_RB_B,Sig_RB,N);    % Category B raw random samples in RB condition

% Discard outliers
j=1; k=1;
for i=1:N
    dA=sqrt((RB_A(i,:)'-mu_RB_A)'*inv(Sig_RB)*(RB_A(i,:)'-mu_RB_A));
    if dA>3
        Out_A(j)=i;
        j=j+1;
    end
    dB=sqrt((RB_B(i,:)'-mu_RB_B)'*inv(Sig_RB)*(RB_B(i,:)'-mu_RB_B));
    if dB>3
        Out_B(k)=i;
        k=k+1;
    end
end

if length(Out_A)>.5         % make sure there is at least 1 outlier
    RB_A(Out_A,:) = [];
end
if length(Out_B)>.5         % make sure there is at least 1 outlier
    RB_B(Out_B,:) = [];
end

% Trim to exactly N_final stimuli in each category
RB_A = RB_A(1:N_final,:);
RB_B = RB_B(1:N_final,:);

% Transform samples to have population mean and variance-covariance matrix
S_A = cov(RB_A); % Sample variance-covariance matrix of category A stimuli
S_B = cov(RB_B); % Sample variance-covariance matrix of category B stimuli
P_RB = (chol(Sig_RB))';  % Cholesky matrix for RB condition
Q_A = (chol(S_A))';  % Q matrix for category A
Q_B = (chol(S_B))';  % Q matrix for category B
for i=1:N_final
    RB_A_pop(i,:) = (P_RB*inv(Q_A)*(RB_A(i,:)'-mean(RB_A)'))' + mu_RB_A';
    RB_B_pop(i,:) = (P_RB*inv(Q_B)*(RB_B(i,:)'-mean(RB_B)'))' + mu_RB_B';
end

% Plot the RB stimuli
scatter(RB_A_pop(:,1),RB_A_pop(:,2)); 
axis([0 100 0 100]); axis square;
hold on;
scatter(RB_B_pop(:,1),RB_B_pop(:,2)); 

DataA=[ones(N_final,1) RB_A_pop];
DataB=[2*ones(N_final,1) RB_B_pop];
Data(1:N_final,:)=DataA;
Data(N_final+1:2*N_final,:)=DataB;

dlmwrite('DataRB.txt',Data,'delimiter', '\t');  %write the RB stimuli in a file

% Rotate the RB stimuli 45 degrees to create the II stimuli
RB_A_zero=RB_A_pop-(ones(300,1)*mean(RB_A_pop));
RB_B_zero=RB_B_pop-(ones(300,1)*mean(RB_B_pop));

R=[cos(pi/4) -sin(pi/4); sin(pi/4) cos(pi/4)];
II_A_zero = (R * RB_A_zero')';
II_B_zero = (R * RB_B_zero')';

II_A_pop=II_A_zero + (ones(300,1)*mu_II_A');
II_B_pop=II_B_zero + (ones(300,1)*mu_II_B');

DataIIA=[ones(N_final,1) II_A_pop];
DataIIB=[2*ones(N_final,1) II_B_pop];
DataII(1:N_final,:)=DataIIA;
DataII(N_final+1:2*N_final,:)=DataIIB;

dlmwrite('DataII.txt',DataII,'delimiter', '\t');  %write the II stimuli in a file

% Plot the II stimuli
figure;
scatter(II_A_pop(:,1),II_A_pop(:,2)); 
axis([0 100 0 100]); axis square;
hold on;
scatter(II_B_pop(:,1),II_B_pop(:,2)); 