%% Computer Exercise1 3.5 Estimation of Real Data
load svedala
%%
y = svedala;
subplot(3,1,1)
acf(y,100,0.05, 1);
ylabel('ACF') 
subplot(3,1,2)
pacf(y,100, 0.05, 1);
ylabel('PACF')
subplot(3,1,3)
normplot(y)
%% Remove the trend
A1 = [1 -1];
y = filter(A1,1,y);
y = y(2:end);
% subplot(3,1,1)
% acf(u,100,0.05, 1);
% ylabel('ACF') 
% subplot(3,1,2)
% pacf(u,100, 0.05, 1);
% ylabel('PACF')
% subplot(3,1,3)
% normplot(u)
%% Remove the seasonality
A24 = [1 zeros(1,23) -1];
z = filter(A24, 1, y);
z = z(25:end);
subplot(3,1,1)
acf(z,100,0.05, 1);
ylabel('ACF') 
subplot(3,1,2)
pacf(z,100, 0.05, 1);
ylabel('PACF')
subplot(3,1,3)
normplot(z)
%% Remove AR(1) param
data = iddata(z);
model_init = idpoly([1 0],[],[]);
model_armax = pem(data, model_init);
z1 = filter(model_armax.a, model_armax.c,z);
subplot(3,1,1)
acf(z1,100,0.05, 1);
ylabel('ACF') 
subplot(3,1,2)
pacf(z1,100, 0.05, 1);
ylabel('PACF')
subplot(3,1,3)
normplot(z1)

%% Remove the MA(24)
data = iddata(z);
model_init = idpoly([1 0],[],[1 zeros(1,24)]);
model_init.Structure.c.Free = [zeros(1,24) 1];
model_armax = pem(data, model_init);
z2 = filter(model_armax.a, model_armax.c,z);
subplot(3,1,1)
acf(z2,100,0.05, 1);
ylabel('ACF') 
subplot(3,1,2)
pacf(z2,100, 0.05, 1);
ylabel('PACF')
subplot(3,1,3)
normplot(z2)  %Param with model_armax in command window.

%Lower fit and FPE when not removing the seasonality.
%Without Season and just AR(1) removed : FPE = 0.4185, 25.91%
%Without season and AR and MA: FPE 0.4106 , 26.66%
%Removed everything: FPE 0.3761, 29.46%
%Didn't remove trend: FPE 0.521, 70.99%

%Didn't get a better estimate with not differentiating.
%What SARIMA PROCESS IS THIS...?
