%% Kalman interpol, Task A:
load proj18.mat
% rain_org_EG = getfield(ElGeneina,'rain_org'); %actual rain data
% rain_org_t_EG = getfield(ElGeneina,'rain_org_t');
% rain_EG = getfield(ElGeneina,'rain');
rain_org_EG = getfield(Kassala,'rain_org'); %actual rain data
rain_org_t_EG = getfield(Kassala,'rain_org_t');
rain_EG = getfield(Kassala,'rain');
%%
u = rain_org_EG(1:end);
% S = 12;
% A12 = [1 zeros(1,S-1) -1];
% A36 = [1 zeros(1,35) -1];
% us = filter(A12,1,u);
% us = us(13:end);
%% Kalman filter! 
% u = us; 
N = length(u);
% Make a guess from ARMA
Re = eye(1)*1.15; % model error 
Rw = 0.2; % measurement error
% Guess a1
a1 = -0.22; 
A = -a1; 
Rxx = eye(1); % Initial variance
xtt = 0;  % Initial guess state 
C = 1; 
Zsave=zeros(3,N);
Usave=zeros(1,N);

for k=1:N
    
  Rxxtt = Rxx-Kt*Ryy10*Kt'; 
  Rxx = A*Rxxtt*A'+Re;
    x = [0; 0; 0];

    for j = 1:3
        Ryy10 = C*Rxx*C'+Rw;
        Kt = (Rxx*C')/Ryy10;
        x(j) = C*A^(j-1)*xtt; %? 
    end
    % the sum of these should add up to u.. 
    error = u(k)-sum(x); 
    % there will be an error.. split this on the three and redo estimate
    for j = 1:3
        x(j) = xtt+Kt*error; 
        xtt = A*x(j); 
    end
    
  % Save state
  Zsave(:,k)=x;
  Usave(k)=[1 1 1]*x; 
end


%% PLOT HERE
u_linear = rain_EG(1:end);
% u_linear = u_linear - mean(u_linear);
% u_linear = filter(A36,1,u_linear);
% u_linear = u_linear(37:end);
plot(u_linear) % fr?ga om detta..
hold on 
plot(Zsave(:))

ucomp = zeros(length(u),3); 

ucomp(:,2) = u; 
ucomp(:,1) = NaN(length(u),1); 
ucomp(:,3) = NaN(length(u),1); 
ucomp = ucomp'; 
ucomp = ucomp(:); 
ucomp(ucomp==0)=NaN; 

plot(ucomp, 'or')
sum(u)
sum(Usave)

rain_Kalman = Zsave(:);
rain_Kalman(rain_Kalman<0) = 0;