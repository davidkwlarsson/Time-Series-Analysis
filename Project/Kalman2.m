%% Task A
load proj18.mat
rain_org_EG = getfield(ElGeneina,'rain_org'); %actual rain data
rain_org_t_EG = getfield(ElGeneina,'rain_org_t');
rain_EG = getfield(ElGeneina,'rain');
% rain_org_EG = getfield(Kassala,'rain_org'); %actual rain data
% rain_org_t_EG = getfield(Kassala,'rain_org_t');
% rain_EG = getfield(Kassala,'rain');
u = rain_org_EG(1:end);


%%
a1 = -0.3;
A = (-a1);
Re = eye(3)*50; % State noise variance - CHANGE: The higher V1 -> higher variance in result
Rw = 0; % Measure variance - CHANGE: The higher V2-> lower variance in result
N = length(u);
% Vector to store values in
Xsave=zeros(3,N);
Usave=zeros(1,N);
Rxxt = 1; 
Xtt = zeros(3,1);
C = [1 1 1];
Ryyt = C*Rxxt*C' + Rw;
% Kalman filter. Start from k=3, since we need old values of y.

for k=1:N
    for i = 1:2
        if (i == 1)
            xtt = Xtt(3,1);
        else 
            xtt = Xtt(1,1);
        end
        Rxx = Rxxt;
        xtt1 = A*xtt;
        Kt = Ryyt\(Rxxt*C');
        Rxxt = A*Rxx*A' + Re;
        Ryyt = C*Rxxt*C' + Rw;
        Xtt(i,1) = xtt1;
    end
    xtt = Xtt(2,1);
    Rxx = Rxxt - Kt*Ryyt*Kt';
    xtt1 = A*xtt;
    Kt = Ryyt\(Rxxt*C');
    Rxxt = A*Rxx*A' + Re;
    Ryyt = C*Rxxt*C' + Rw;
    Xtt(3,1) = xtt1;
    Xtt = Xtt+Kt*(u(k)-C*Xtt);
    Xtt(Xtt<0) = 0;
  % Save state
  Xsave(:,k)= Xtt;
  Usave(k) = C*Xtt;
end

% Zsave = Zsave(:);
% plot(Xsave(:))
plot(Usave, 'b')
hold on
plot(rain_org_EG,'*r')
sum(rain_org_EG)
sum(Xsave(:))
rain_Kalman = Xsave(:);

% QUESTIONS: 
% DO WE TAKE AWAY THE MA(12) PART ASWELL, WHY ONLY AR(1)?
% IF WE TAKE THE SQRT DOES THE EQUALITY HOLD, Y = X1 + X2 + X3?