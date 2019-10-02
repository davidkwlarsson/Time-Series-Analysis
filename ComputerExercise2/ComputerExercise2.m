%% Computer Exercise 2

n = 500;
A1 = [1 -.65]; A2 = [1 .90 .78];
C = 1 ; 
B = [0 0 0 0 .4];
e = sqrt(1.5)*randn(n+100,1);
w = sqrt(2)*randn(n+200,1);
A3 = [1 .5];
C3 = [1 -.3 .2];
u = filter(C3,A3,w);
u = u(101:end);
y = filter(C,A1,e)+filter(B,A2,u);
u = u (101:end);
y = y (101:end);
% clear A1, A2, C, B, e, w, A3, C3

%% Q.1 ARMA model for input u 
data = iddata(u);
model_init = idpoly([1 0],[],[]); %AR(1), white!
model_armax = pem(data, model_init);
upw = filter(model_armax.a, model_armax.c,u);
A3 = model_armax.a;
C3 = model_armax.c;
whitenessTest(upw)
ypw = filter(model_armax.a, model_armax.c,y);
%%
M = 40;
stem(-M:M,crosscorr(upw,ypw,M)); %d=4 r=2 s=1
title('Crosscorrelationfunction'),xlabel('Lag') 
hold on
plot(-M:M,2/sqrt(length(u))*ones(1,2*M+1),'--')
plot(-M:M,-2/sqrt(length(u))*ones(1,2*M+1),'--')
hold off
%% Q.2
A2 = [1 0 0]; %to -r 
B = [0 0 0 0 0]; %-2 to -4
Mi = idpoly([1],[B],[],[],[A2]);
Mi.Structure.b.Free = [0 0 0 0 1];
zpw = iddata(ypw,upw);
Mba2 = pem(zpw,Mi); 
present(Mba2)
vhat = resid(Mba2,zpw,'corr');
v_res = vhat.OutputData; %Not white noise and we dont think it should be
stem(-M:M,crosscorr(upw,v_res,M));  %Look at the definition for v(t)...
%% Q.3
A2 = Mba2.f;
B = Mba2.b; %This B already contains the delay
x = y - filter(B,A2,u);
stem(-M:M,crosscorr(x,u,M));
CCF_xu = crosscorr(u,x,M); %This is not white.
%%
A1 = [1 0 0 0 0]; 
C1 = [1 0 0];
Mi = idpoly(A1,[],C1);
Mi.Structure.a.Free = [0 0 0 1 1];
Mi.Structure.c.Free = [0 1 1];
xdata = iddata(x);
model_armax = pem(xdata,Mi);
x_res = filter(model_armax.a, model_armax.c, x);
whitenessTest(x_res)

%% Q.4
A1 = [1 0 0 0 0];
A2 = [1 0 0];
B = [0 0 0 0 0];
C = [1 0 0];
Mi = idpoly(1,B,C,A1,A2);
%Mi.Structure.b.Free = [0 0 1 1 1];
Mi.Structure.b.Free = [0 0 0 0 1];
% Mi.Structure.f.Free = [0 0 1]; %f is A2
Mi.Structure.d.Free = [0 0 0 1 1]; %d is A1
Mi.Structure.c.Free = [0 1 1];
z = iddata(y,u);
MboxJ = pem(z,Mi);
present(MboxJ)
ehat=resid(MboxJ,z);
e_res = ehat.OutputData;
whitenessTest(e_res)
crosscorr(u,e_res);

%The residual is now white noise, however it dosent seem completely
%uncorrelated with the input signal. Some parameters in B seems to be
%very close to zero and the 4:th in D,(A1) is also kind of small
%Try setting the parameters to zero, that are already close to zero.
%How do we know if the parameters are significantly different from zero?

%% Q.5
load tork.dat
tork = tork - repmat(mean(tork),length(tork),1);
y = tork(:,1); u = tork(:,2);
% u=u(101:end);
% y = y(101:end);
z = iddata(y,u);
plot(z(1:300))
A = [1 0];
Mi = idpoly(A,[],[]);
model_armax = pem(z.inputData,Mi);
upw = filter(model_armax.a, model_armax.c,u);
ypw = filter(model_armax.a, model_armax.c,y);
upw = upw(5:end);
ypw = ypw(5:end);
%%
M = 40;
stem(-M:M,crosscorr(upw,ypw,M)) % d=3 s=2 r=2
hold on
plot(-M:M,2/sqrt(length(w))*ones(1,2*M+1),'--')
plot(-M:M,-2/sqrt(length(w))*ones(1,2*M+1),'--')
hold off
%%
A2 = [1 0 0]; %to -r 
B = [0 0 0 0 0 0]; %-3 to -5
Mi = idpoly([1],[B],[],[],[A2]);
Mi.Structure.b.Free = [0 0 0 1 1 1];
zpw = iddata(ypw,upw);
Mba2 = pem(zpw,Mi); 
present(Mba2)
vhat = resid(Mba2,zpw,'corr');
v_res = vhat.OutputData; %Not white noise and we dont think it should be
% v_res = v_res(10:end);
% upw = upw(10:end);
stem(-M:M,crosscorr(upw,v_res,M));
%%
x = y - filter(Mba2.b,Mba2.f,u);
crosscorr(x,u)
%%
A1 = [1 0]; 
Mi = idpoly(A1,[],[]);
xdata = iddata(x);
model_armax = pem(xdata,Mi);
x_res = filter(model_armax.a, model_armax.c, x);
whitenessTest(x_res)
%% Reestimate all param
A1 = [1 0];
A2 = [1 0 0];
B = [0 0 0 0 0 0]; %from 0 to -s or -d to -(d+s)
Mi = idpoly(1,B,[1],A1,A2);
Mi.Structure.b.Free = [0 0 0 1 1 1];
z = iddata(y,u);
MboxJ = pem(z,Mi);
present(MboxJ)
resid(MboxJ,z,'corr')
ehat=resid(MboxJ,z,'corr');
figure
e_res = ehat.OutputData;
whitenessTest(e_res)
%Nearly uncorrelated and the residual is white.

%Verify that u to x is white? Does this have to be the case?
%Why is this weird?
%%
load svedala
y = svedala;
k = 3;
A = [1 -1.79 0.84];
C = [1 -0.18 -0.11];
[CS,AS] = equalLength(C,A);
[Fk,Gk] = deconv(conv([1,zeros(1,k-1)],CS),AS);
yhat_k = filter(Gk,C,y);
yhat_k = yhat_k(30:end);
e_res = svedala(30:end) - yhat_k; 
%k=1 gives white noise
%the variance will be the coefficients squared up to k-1
%mean(e_res) % = 0.0810, k =3 or 0.2480, k = 26
%Expectation of F(z)e_(t+k) = 0
N_var = Fk*Fk'*sigmae;
[r,tau] = kovarians(Gk,C,2);
r(1) %% 33.8668
%33.6 for k = 26 and 7.3243 for k = 3
%difference between theo and estimated.
conf = 1.96*sqrt(N_var);
hold on
sum(e_res>conf) % zero outside confidence intervall
%For k = 3 we get an ma2 process but k = 26 just gives gives ar process
%% 3.4
load sturup
u = sturup;
A = [1 -1.49 0.57];
B = [0 0 0 0.28 -0.26]; %number of zeros in b = delay
C = [1];
[CS,AS] = equalLength(C,A);
[Fk,Gk] = deconv(conv([1,zeros(1,k-1)],CS),AS);
yhat_k = filter(Gk,C,y);
yhat_k = yhat_k(101:end);
BF = conv(B,Fk);
[BFS,CS] = equalLength(BF,C);
[Fku,Gku] = deconv(conv([1,zeros(1,k-1)],BFS),CS);
uhat_k = filter(Gku,C,u);
uhat_k = uhat_k(101:end);
Yhat_k = uhat_k + yhat_k;
pred_e = svedala(30:end)-Yhat_k;
e_res = svedala(30:end)-yhat_k;
% k = 3 gives the same variance as for the normal arma
%since k =< d and therefore all x-values will be in the past
%not contributing to errors.
ARMAX_var = Fk*Fk'*std(pred_e)^2;
hold on
plot(Yhat_k, 'r')
plot(svedala(101:end))
ylabel('armax')
figure
hold on
plot(yhat_k, 'r')
plot(svedala(101:end))
ylabel('without input')
%The prediction error has a lower norm value
%The theoretical is fron the F ploynomial while the estimated is from STD^2
%% 3.5 Q.14
%period = 24
S = 24;
AS = [1 zeros(1,S-1) -1];
Mi = idpoly(AS,[],[]);
Mi.Structure.a.Free = [zeros(1,S) 1];
data = iddata(svedala);
model_armax = pem(data,Mi);
ys = filter(model_armax.a,model_armax.c,svedala);
%%
A = [1 0 0];
C = [1 zeros(1,24)];
Mi = idpoly(A,[],C);
Mi.Structure.c.Free = [zeros(1,24) 1];   
data = iddata(ys);
model_armax = pem(data,Mi);
w = filter(model_armax.a,model_armax.c,ys);
acfpacf(w)
whitenessTest(w)
%%
k = 3;
Atot = conv(A,AS);
[CS,AS] = equalLength(C,Atot);
[Fk,Gk] = deconv(conv([1,zeros(1,k-1)],CS),AS);
yhat_k = filter(Gk,C,y);
yhat_k = yhat_k(101:end);
res = svedala(101:end) - yhat_k;
hold on
plot(svedala(101:end),'b')
plot(yhat_k,'r')
e_var = std(res)^2;