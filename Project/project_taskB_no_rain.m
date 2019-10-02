%% Script to run the project. 
load proj18.mat
%% B (ElGeneina):
nvdi_EG = getfield(ElGeneina,'nvdi'); %this seems to be AR(1)
nvdi_t_EG = getfield(ElGeneina,'nvdi_t'); 
rain_EG = getfield(ElGeneina,'rain'); %This is the interpolated data
rain_org_EG = getfield(ElGeneina,'rain_org'); %actual rain data
rain_org_t_EG = getfield(ElGeneina,'rain_org_t');
nvdi_EG_2 = getfield(ElGeneina,'nvdi'); %this seems to be AR(1)

%% Divide into Model,Validation and test data:
nvdi_EG = (nvdi_EG/255)*2 - 1; % [0,255] ->[-1,1]
y = nvdi_EG(1:450);
y_valid = nvdi_EG(451:end);
y_test = nvdi_EG(581:end);
u = rain_Kalman(793:1242);
u_valid = rain_Kalman(1243:end);
u_test = rain_Kalman(1373:end);

%% Without percipation:
mu_y = mean(y);
%y0 = y; 
y0 = y - mu_y;
%% trend - we assume no trend! Maybe we should make a test or linear regression and check significance
acfpacf(y0)
%% seasonality
%data = iddata(y0);
S = 36; 
A36 = [ 1 zeros(1,S-1) -1];
% model_init = idpoly(A36 ,[], []) ;
% model_init.Structure.a.Free = [zeros(1,S) 1];
% model = pem(data, model_init);
% present(model)
% A36 = model.a;
%ys = filter(model.a, model.c, yt);
ys = filter(A36, 1, y0); 
ys = ys(37:end);
acfpacf(ys,100)
%%
subplot(131)
plot(y)
title('Given vegetation data, rescaled to [-1,1]')
xlim([0 450])
subplot(132)
plot(linspace(37, 450, length(ys) ),ys)
xlim([0 450])
title('Removed seasonality and made zero mean')
subplot(133)
normplot(y)

%%
n = 50; 
subplot(2,2,1)
acf(ys,n,0.05, 1);
ylabel('ACF')
title('Before')
subplot(2,2,3)
pacf(ys,n, 0.05, 1);
ylabel('PACF')

%%
data = iddata(ys); 
A = [1 0];
C = [1 zeros(1,S)];
%C = []; 
model_init = idpoly(A ,[], C) ;
model_init.Structure.c.Free = [ zeros(1,S) 1];
model = pem(data, model_init);
present(model)
yar = filter(model.a, model.c, ys);
yar = yar(2:end); 
acfpacf(yar)

subplot(2,2,2)
acf(yar,n,0.05, 2);
ylabel('ACF')
title('After')
subplot(2,2,4)
pacf(yar,n, 0.05, 1);
ylabel('PACF')

%% k-step without rain:
k = 5;
y_validm = y_valid - mu_y; 
A = conv(model.a,A36);
C = [1];
[CS,AS] = equalLength(C,A);
[Fk,Gk] = deconv(conv([1,zeros(1,k-1)],CS),AS);
yhat_k = filter(Gk, C, y_validm);
yhat_k = yhat_k(length(Gk)+10:end);
y_validm = y_validm(length(Gk)+10:end); 
res = y_validm(k:end)-yhat_k(k:end);
subplot(2,2,4)
acf(res,40,0.05, 1);
 %this should be a MA(k-1)
%whitenessTest(res) %run this for k=1 to check if white
std_norain = std(res)
subplot(2,2,2)
hold on
val_plot =y_valid(length(Gk)+10:end);
plot(val_plot(k:end), 'r')
plot(yhat_k(k:end)+mu_y)
title('One-step prediction')
hold off
mse = sum(res.^2)/length(res)