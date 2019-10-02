%% Task A
load proj18.mat
rain_org_EG = getfield(ElGeneina,'rain_org'); %actual rain data
rain_org_t_EG = getfield(ElGeneina,'rain_org_t');
rain_EG = getfield(ElGeneina,'rain');
u = rain_org_EG(1:end);


%%
a1 = -0.5;
A = (-a1);
Re = 1.5; % State noise variance - CHANGE: The higher V1 -> higher variance in result
Rw = 0; % Measure variance - CHANGE: The higher V2-> lower variance in result
N = length(u);
% Vector to store values in
Zsave=zeros(3,N);
Ysave=zeros(1,N);
Rxxt = 1; 
Xtt = zeros(3,1);
C = [1 1 1];
Ryyt = C*Rxxt*C' + Rw;
% Kalman filter. Start from k=3, since we need old values of y.

for k=1:N
    for i = 1:3
        xtt = Xtt(1,1);
        Rxx = Rxxt;
        xtt1 = A*xtt;
        Kt = Ryyt\(Rxxt*C');
        Rxxt = A*Rxx*A' + Re;
        Ryyt = C*Rxxt*C' + Rw;
        Xtt(i,1) = xtt1;
    end
  xtt = Xtt+Kt*(u(k)-C*Xtt);
  xtt(xtt<0) = 0;
  % Save state
  Xsave(:,k)= xtt;
  Usave(k) = C*xtt;
end

% Zsave = Zsave(:);
plot(Xsave(:))

sum(rain_org_EG)
sum(Xsave(:))
rain_Kalman = Xsave(:);

% QUESTIONS: 
% DO WE TAKE AWAY THE MA(12) PART ASWELL, WHY ONLY AR(1)?
% IF WE TAKE THE SQRT DOES THE EQUALITY HOLD, Y = X1 + X2 + X3?