%% B-with KALMAN RAIN:
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
y = nvdi_EG(1:450);
y_valid = nvdi_EG(451:end);
y_test = nvdi_EG(581:end);
u = rain_Kalman(793:1242);
u_valid = rain_Kalman(1243:end);
u_test = rain_Kalman(1373:end);

mu_y = mean(y);
y0 = y - mu_y;
S = 36; 
A36 = [ 1 zeros(1,S-1) -1];
ys = filter(A36, 1, y0); 
ys = ys(37:end);
y_validm = y_valid - mu_y; 
% u_valid = u_valid - mean(u_valid);
u_valid1 = sqrt(u_valid);
%%
uf = sqrt(u);
u_validm = u_valid1 - mean(uf);
uf = uf - mean(uf);

S = 36; 
A36 = [1 zeros(1,S-1) -1];
us = filter(A36, 1,uf);
us = us(37:end);
data = iddata(us);
S = 36;
C3 = [1 zeros(1,S)];
model_init = idpoly([1 0 0],[],C3);
model_init.Structure.c.Free = [0 0 0 1 zeros(1,S-4) 1];
model_armax = pem(data, model_init);
upw = filter(model_armax.a, model_armax.c,us); 
upw = upw(3:end);
ypw = filter(model_armax.a, model_armax.c,ys); %the no season and zero mean process y
ypw = ypw(3:end);
crosscorr(upw,ypw,50)
%%
A2 = [1 0 0]; %to -r 
B = [0 0 0]; %-d to -(d+s)
Mi = idpoly([1],[B],[],[],[A2]);
Mi.Structure.b.Free = [0 0 1];
zpw = iddata(ypw,upw);
Mba2 = pem(zpw,Mi); 
% present(Mba2)
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
C = [1]; %model.c
Mi = idpoly(1,B,C,A1,A2);
Mi.Structure.b.Free = [0 0 1 0 0 0];
% Mi.Structure.c.Free = [zeros(1,36) 1];
z = iddata(ys,us);
model = pem(z,Mi);
Au = conv(model.d,model.b);
Cu = conv(model.f,model.c);
e_res = filter(model.d,model.c,ys) - filter(Au,Cu,us);
e_res = e_res(8:end);
whitenessTest(e_res)
% crosscorr(e_res,us(8:end));
%% k-step prediction with rain:
k = 5;
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

Yhat = yhat_k + uhat_k + mu_y;
% 
plot(Yhat(k:end),'b')
hold on
plot(y_valid(43+k:end),'r')

% plot(yhat_k(k:end),'--r')
% plot(uhat_k(k:end),'--b')
res_naive = (y_valid(44:end)-y_valid(44-k:end-k));
MSE_naive = sum(res_naive.^2)/length(res_naive)
STD_naive = std(res_naive)
res = (Yhat(k:end) - y_valid(43+k:end));
% acf(res,50,0.05, 1);
MSE_kalman = sum(res.^2)/length(res)
STD_kalman = std(res)

% %% Apply the model on Kassala:
% % 
% nvdi_K = getfield(Kassala,'nvdi'); %this seems to be AR(1)
% nvdi_t_K = getfield(Kassala,'nvdi_t'); 
% nvdi_K = nvdi_K./255*2 - 1; % [0,255] ->[-1,1]
% y_valid = nvdi_K(451:end);
% u = rain_Kalman(793:1242);
% u_valid = rain_Kalman(1243:end);
% u_test = rain_Kalman(1373:end);
% 
% mu_y = mean(y);
% y0 = y - mu_y;
% S = 36; 
% A36 = [ 1 zeros(1,S-1) -1];
% ys = filter(A36, 1, y0); 
% ys = ys(37:end);
% y_validm = y_valid - mu_y; 
% u_valid1 = sqrt(u_valid);
% u_validm = u_valid1 - mean(u_valid1);
% %%
% % figure
% k = 5;
% A = conv(model.d,model.f);
% C = conv(model.f,model.c);
% B = conv(model.d,model.b);
% A = conv(A36,A);
% [CS,AS] = equalLength(C,A);
% [Fk,Gk] = deconv(conv([1,zeros(1,k-1)],CS),AS);
% yhat_k = filter(Gk,C,y_validm);
% yhat_k = yhat_k(44:end);
% B = conv(B,A36);
% BF = conv(B,Fk);
% [BFS,CS] = equalLength(BF,C);
% [Fku,Gku] = deconv(conv([1,zeros(1,k-1)],BFS),CS);
% uhat_k = filter(Gku,C,u_validm);
% uhat_k = uhat_k(44:end);
% 
% Yhat = yhat_k + uhat_k + mu_y;
% y_naive = y_valid(44-k:end-k);
% 
% % plot(Yhat(k:end),'b')
% % hold on
% % plot(y_valid(43+k:end),'r')
% 
% % plot(yhat_k(k:end),'--r')
% % plot(uhat_k(k:end),'--b')
% res_naive = (y_valid(44:end)-y_naive);
% MSE_naive = sum(res_naive.^2)/length(res_naive)
% STD_naive = std(res_naive)
% res = (Yhat(k:end) - y_valid(43+k:end));
% % acf(res,50,0.5,1);
% MSE_EGKassala = sum(res.^2)/length(res)
% STD_EG_Kassala = std(res)