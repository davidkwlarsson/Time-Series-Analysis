%% Computer Exercise 3

load tar2.dat
load thx.dat
%%
subplot(311)
plot(tar2)
ylabel('tar2')
subplot(312)
plot(thx(:,1))
subplot(313)
plot(thx(:,2))
%% Q1
lambda = 0.95;
X = recursiveAR(2,'ForgettingFactor', lambda); %AR(2)
z = iddata(tar2);
Aest = zeros(3,500);
Yhat = zeros(1,500);
for kk = 1:numel(z.y)
    [A,yhat] = step(X,z.y(kk));
    Aest(:,kk) = A(1,1:3)';
    Yhat(kk) = yhat;
end
sysX = idpoly(X); %A(z) = 1 + 1.041z^-1 + 0.2948z^-2

subplot(221)
plot(thx(:,1))
subplot(222)
plot(thx(:,2))
subplot(223)
plot(Aest(2,:))
subplot(224)
plot(Aest(3,:))
% Lambda will make the values previous in time matter less to the new est.
% i.e the lower the lambda the faster past values become insignificant,
% faster change but higher variance.

%% Q.2

n = 100;
lambdaline=linspace(0.85,1,n);
ls2 = zeros(n,1);
for i=1:length(lambdaline)
    X = recursiveAR(2,'ForgettingFactor', lambdaline(i)); %AR(2)
    for kk = 1:numel(z.y)
        [A,yhat] = step(X,z.y(kk));
        Aest(:,kk) = A(1,1:3)';
        Yhat(kk) = yhat;
    end
    ls2(i)=sum((tar2-Yhat').^2);
end

%%
figure
plot(lambdaline,ls2) %ls2 contains the squared error for the predictions
                     %from using the different lambdas. The minimum shows
                     %the optimal, least squares estimate of lambda.
min_index = find(ls2 == min(ls2));
lambda_min = lambdaline(min_index); %This is 0.9424
%% Q.3
% Length of process
y = tar2;
N = length(tar2);
% Set parameters
A = [1 0;
     0 1];
V1 = [1e-4 0;
      0 1e-4]; % State noise variance
V2 = 1; % Measure variance
%usually C should be set here to, but in this case C is a function of
%time.
% Set initial values
Vtt = 1*eye(2); % Initial variance
% Ztt = [mean(Aest(2,:)); mean(Aest(3,:))]; % Initial state (Initial guess?)
Ztt = [1; 1];
% Vector to store values in
Zsave=zeros(2,N);

% Kalman filter. Start from k=3, since we need old values of y.
for k=3:N
  % C is a function of time.
  C = [-y(k-1) -y(k-2)];
  % Time update
  %Ztt_1 = Ztt + Kt*(y(k)-C*Ztt);
  Ztt_1 = A*Ztt;
  Vtt_1 = A*Vtt*A' + V1; %R_xx for t+1
  % Measure update
  Vtt = C*Vtt_1*C' + V2; %R_yy for t+1
  Kt = (Vtt_1*C')/Vtt;
  Ztt = Ztt_1 + Kt*(y(k)-C*Ztt_1);
  % Save the state
  Zsave(:,k)=Ztt;
  Yhat(:,k) = C*Ztt;
end
hold on
plot(thx(:,1), 'b')
plot(Zsave(1,:), 'r')
figure
hold on
plot(thx(:,2),'b')
plot(Zsave(2,:),'r')

%Lower R_e will give lower variance for the parameters
%Small R_w will give us higher variance (this is becuase the Kalman gain
%will be lower)

%% Q.4
b = 20;
y = tar2 + b*states(1:500)';
u = states(1:500);
simgae_2 = 1;
simgav_2 = 4;
% Length of process
N = length(tar2);
% Set parameters
A = [1 0;
     0 1];
V1 = [1 0;
      0 1]; % State noise variance - CHANGE: The higher V1 -> higher variance in result
V2 = 4; % Measure variance - CHANGE: The higher V2-> lower variance in result
%usually C should be set here to, but in this case C is a function of
%time.
% Set initial values
Vtt = eye(2); % Initial variance
Ztt = [1 1]'; % Initial state (Initial guess?)
% Vtt = A*Vtt*A' + V1;
% Vtty = C*Vtt*C' + V2;
% Kt = (Vtt*C')/Vtty;
% Vector to store values in
Zsave=zeros(2,N);
Ysave=zeros(1,N);

% Kalman filter. Start from k=3, since we need old values of y.

for k=3:N
  % C is a function of time.
  %C = [1 0];
  C = [1 -u(k)];
  % Update Kalman gain 
  Vtt = A*Vtt*A' + V1; 
  Vtty = C*Vtt*C' + V2;
  Kt = (Vtt*C')/Vtty;

  % Time update and Measure ?? Which is what? 
  Ztt_1 = A*Ztt;
  Ztt = Ztt_1 + Kt*(y(k) - C*Ztt_1); 
 
  % Save state
  Zsave(:,k)=Ztt;
  Ysave(k) = C*Ztt;
end
% figure
subplot(211)
plot(Ysave)
subplot(212)
plot(y,'r')
MSE = sum((y(1:end)-Ysave(1:end)').^2)
%Why does A become the identity matrix?
%R_e is the error for the model, noise, hence we chose 1(given in task).
%R_v is the error of the measuremnt which is given as 4.
% plot(Zsave(3,:))
figure
hold on
plot(Zsave(1,:), 'r');
plot(tar2,'b');
%% Q.6
load svedala94.mat
T = linspace(datenum(1994,1,1),datenum(1994,12,31),length(svedala94));
plot(T,svedala94) ;
datetick('x');
%%
y = svedala94;
Season = [1 0 0 0 0 0 -1];
% ys = y;
ys = filter(Season, 1, y);
% ys = ys(7:end);
th = armax(ys,[2 2]);
th_winter = armax(ys(1:540),[2 2]);
th_summer = armax(ys(907:1458),[2 2]); %Different coeff.
%% Q.7
th0 = [th_winter.A(2:end) th_winter.C(2:end)];
[thr,yhat] = rarmax(ys,[2 2],'ff',0.99,th0);
subplot(311)
plot(T,svedala94);
datetick('x')
subplot(312)
plot(thr(:,1:2))
hold on
plot(repmat(th_winter.A(2:end),[length(thr) 1]),'b:');
plot(repmat(th_summer.A(2:end), [length(thr) 1]),'r:');
axis tight
hold off

subplot(313)
plot(thr(:,3:end))
hold on
plot(repmat(th_winter.C(2:end),[length(thr) 1]),'b:');
plot(repmat(th_summer.C(2:end),[length(thr) 1]),'r:');
axis tight
hold off
%ANS : This seems to coincide with the summer process.
%% Q.8
load svedala94
y = svedala94(850:1100);
y = y-mean(y);
t=(1:length(y))';
U=[sin(2*pi*t/6) cos(2*pi*t/6)];
Z = iddata(y,U);
model = [3 [1 1] 4 [0 0]];
thx = armax(Z,model); %Different season here gives a worse estimate of B poly
Sf = U*cell2mat(thx.b)'; %less part of y is described by b*U
hold on
plot(y)
plot(Sf,'r') 
%Can one add more b coeffs to better model the season?

%% Q.9
% y = svedala94(850:1100);
y = svedala94;
y = y - y(1);
% y = y - mean(y);
t=(1:length(y))';
U = [sin(2*pi*t/6) cos(2*pi*t/6) ones(size(t))];
Z = iddata(y,U);
A = thx.A(2:end);
B = thx.B;
C = thx.C(2:end);
m0 = [A B 0 C];
m0 = cell2mat(m0);
Re = diag([0 0 0 0 0 0.7 0 0 0 0]);
model = [3 [1 1 1] 4 0 [0 0 0] [1 1 1]];
[thr, yhat] = rpem(Z,model,'kf', Re, m0);

m = thr(:,6);
a = thr(end,4);
b = thr(end,5);
y_mean = m + a*U(:,1) + b*U(:,2);
y_mean = [0; y_mean(1:end-1)];
hold on
plot(y_mean, 'r')
plot(y,'b')

%Seems to be a good estimate...
