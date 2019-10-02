%% Script to run the project. 
load proj18.mat
%% B (ElGeneina):
nvdi_EG = getfield(ElGeneina,'nvdi'); %this seems to be AR(1)
nvdi_t_EG = getfield(ElGeneina,'nvdi_t'); 
rain_EG = getfield(ElGeneina,'rain'); %This is the interpolated data
rain_org_EG = getfield(ElGeneina,'rain_org'); %actual rain data
rain_org_t_EG = getfield(ElGeneina,'rain_org_t');
%% Divide into Model,Validation and test data:
nvdi_EG = nvdi_EG./255*2 - 1; % [0,255] ->[-1,1]
y = nvdi_EG(1:500);
y_valid = nvdi_EG(501:580);
y_test = nvdi_EG(581:end);
u = rain_EG(793:1292);
u_valid = rain_EG(1293:1372);
u_test = rain_EG(1373:end);
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
%% k-step without rain:
k = 5;
y_validm = y_valid - mu_y; 
A = conv(model.a,A36);
C = [1];
[CS,AS] = equalLength(C,A);
[Fk,Gk] = deconv(conv([1,zeros(1,k-1)],CS),AS);
yhat_k = filter(Gk, C, y_validm);
yhat_k = yhat_k(38+k:end);
res = y_validm(38+k:end)-yhat_k;
acfpacf(res,20) %this should be a MA(k-1)
%whitenessTest(res) %run this for k=1 to check if white
V_norain = std(res)^2;
figure
hold on
plot(y_validm(38+k:end), 'r')
plot(yhat_k)
%% Now using the rain: (using the square root of rain)
uf = sqrt(u);
uf = uf - mean(uf);
S = 36; 
A36 = [1 zeros(1,S-1) -1];
us = filter(A36, [1],uf);
us = us(37:end);
data = iddata(us);
C3 = [1 zeros(1,S)];
model_init = idpoly([1 0 0],[],C3);
model_init.Structure.c.Free = [0 0 0 1 zeros(1,S-4) 1];
model_armax = pem(data, model_init);
upw = filter(model_armax.a, model_armax.c,us); 
upw = upw(3:end);
ypw = filter(model_armax.a, model_armax.c,ys); %the no season and zero mean process y
ypw = ypw(3:end);
crosscorr(upw,ypw,40)
%%
A2 = [1 0 0]; %to -r 
B = [0 0 0 0 0 0]; %-d to -(d+s)
Mi = idpoly([1],[B],[],[],[A2]);
Mi.Structure.b.Free = [0 0 0 0 0 1];
zpw = iddata(ypw,upw);
Mba2 = pem(zpw,Mi); 
present(Mba2)
vhat = resid(Mba2,zpw,'corr');
v_res = vhat.OutputData;
crosscorr(upw,v_res)
A2 = Mba2.f;
B = Mba2.b; %This B already contains the delay
%%
x = ys - filter(B,A2,us);
x = x(3:end);
crosscorr(x,upw,40) %this shows that all dependecy of u in the remaining process was removed

%% Reestimate all param
A1 = [1 0]; %model.d
A2 = [1 0 0]; %model.f
B = [0 0 0 0 0 0]; %model.b
C = [1 zeros(1,36)]; %model.c
Mi = idpoly(1,B,C,A1,A2);
Mi.Structure.b.Free = [0 0 0 0 0 1];
Mi.Structure.c.Free = [zeros(1,36) 1];
z = iddata(ys,us);
model = pem(z,Mi);
% present(MboxJ)
Au = conv(model.d,model.b);
Cu = conv(model.f,model.c);
e_res = filter(model.d,model.c,ys) - filter(Au,Cu,us);
e_res = e_res(8:end);
% ehat=resid(MboxJ,z);
% e_res = ehat.OutputData;
whitenessTest(e_res)
crosscorr(e_res,us(8:end));
%% k-step prediction with rain:
u_valid1 = sqrt(u_valid);
u_validm = u_valid1 - mean(u_valid1);
k = 4;
A = conv(model.d,model.f);
C = conv(model.f,model.c);
B = conv(model.d,model.b);
A = conv(A36,A);
[CS,AS] = equalLength(C,A);
[Fk,Gk] = deconv(conv([1,zeros(1,k-1)],CS),AS);
yhat_k = filter(Gk,C,y_validm);
yhat_k = yhat_k(44:end);
B = conv(B,A36);
BF = conv(B,Fk);
[BFS,CS] = equalLength(BF,C);
[Fku,Gku] = deconv(conv([1,zeros(1,k-1)],BFS),CS);
uhat_k = filter(Gku,C,u_validm);
uhat_k = uhat_k(44:end);

Yhat = yhat_k + uhat_k;
y_naive = y_validm(44-k:end-k);
hold on
plot(Yhat(k:end),'b')
plot(y_validm(44:end-k),'r')

plot(yhat_k,'--r')
plot(uhat_k,'--b')
res_naive = (y_validm(44:end)-y_naive);
res = (Yhat - y_validm(44:end));