%% Script to run the project.
load proj18.mat
%Tested removing season 36, didnt work
%Tested an AR(3) process, didnt work


%% B (ElGeneina)
nvdi_EG = getfield(ElGeneina,'nvdi'); %this seems to be AR(1)
nvdi_t_EG = getfield(ElGeneina,'nvdi_t'); 
rain_EG = getfield(ElGeneina,'rain'); %This is the interpolated data
rain_org_EG = getfield(ElGeneina,'rain_org'); %actual rain data
rain_org_t_EG = getfield(ElGeneina,'rain_org_t');
%acfpacf(nvdi_EG,20) %Season at 36 and AR(1)?
nvdi_EG = nvdi_EG./255*2 - 1;
nvdi_model_EG = nvdi_EG(1:549);
nvdi_valid_EG = nvdi_EG(550:end);
%% Without percipiation: I.E ARMA
%Remove seasonality: Dosent make an improvemnt
z_EG = nvdi_model_EG;
% A36 = [1 zeros(1,35) -1];
% z_EG = filter(A36, 1, nvdi_model_EG);
% z_EG = z_EG(25:end);
data = iddata(z_EG);
model_init = idpoly([1 0 0],[],[1 0 0]); %ARMA(2,2) makes the best white
model_init.Structure.c.Free = [0 0 1];   %white noise
model_armax = pem(data, model_init);
z_EG = filter(model_armax.a, model_armax.c,z_EG); 
whitenessTest(z_EG);
%z_EG = filter(model_armax.a, model_armax.c,nvdi_EG); 
%% 1 step method (without rain):
k = 1;
A = model_armax.a;
C = model_armax.c;
[CS,AS] = equalLength(C,A);
[Fk,Gk] = deconv(conv([1,zeros(1,k-1)],CS),AS);
yhat_k = filter(Gk, C, nvdi_valid_EG);
yhat_k = yhat_k(10:end);
res = nvdi_valid_EG(10:end) - yhat_k;
acfpacf(res)
whitenessTest(res)

%% k-step method (without rain):
k = 6;
A = model_armax.a;
C = model_armax.c;
[CS,AS] = equalLength(C,A);
[Fk,Gk] = deconv(conv([1,zeros(1,k-1)],CS),AS);
yhat_k = filter(Gk, C, nvdi_valid_EG);
yhat_k = yhat_k(10:end);
res = nvdi_valid_EG(10:end) - yhat_k;
acfpacf(res)  %this is an MA(7), should be MA(5)

%% With percipation: I.E ARMAX
%Start with modeling the percipation, squareroot.
% rain_model_EG = rain_EG(1:end-100);
rain_model_EG = sqrt(rain_EG(1:end-100));
y = sqrt(nvdi_EG(1:549));
% Trend = [1 -1];
% rain_model_EG = filter(Trend, 1, rain_model_EG);
% data = iddata(rain_model_EG);
% A36 = [1 zeros(1,36)];
% model_init = idpoly(A36,[],[1]);
% model_init.Structure.a.Free = [zeros(1,36) 1];
% model = pem(data, model_init);
% A36 = model.a;
% rain_model_EG = filter(A36, model.c, rain_model_EG); %Take away period
% rain_model_EG = rain_model_EG(36:end);
data = iddata(rain_model_EG);
model_init = idpoly([1 0 0],[],[1 0 0 0]);
% model_init = idpoly([1 0 0],[],[1 zeros(1,36)]);
model_init.Structure.c.Free = [0 0 0 1];
model_armax = pem(data, model_init);
upw = filter(model_armax.a, model_armax.c,rain_model_EG);
% Apw = conv(conv(A36,Trend),model_armax.a);
ypw = filter(model_armax.a, model_armax.c,y);
whitenessTest(upw) %this clears two tests
upw = upw(792:end); %21*36 = 756
%% Check the crosscorrelation, upw and ypw
M = 40;
stem(-M:M,crosscorr(upw,ypw,M)); % d = 4, s = 0, r = 0,1 or 2
title('Crosscorrelationfunction'),xlabel('Lag') 
%%
% vhat = ypw - upw; 
% whitenessTest(vhat) %this becomes white?
%%