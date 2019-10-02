%% ComputerExercise1
% A = [1 a1 0 a3];
S = 5;
% C = [1 c1 zeros(1,10) c12];

model_init = idpoly([1 0 0 0], [],[1 zeros(1,12)]);
model_init.Structure.a.Free = [0 1 0 1];
model_init.Structure.c.Free = [0 1 zeros(1,10) 1];
model_armax = pem(data,model_init);

% model_init = idpoly([1 0], [],[]);
% model_armax = pem(data,model_init);
% model_init = idpoly([1 armax(2) 0 0], [],[]);
% model_init.Structure.a.Free = [0 1 0 1];
% model_armax = pem(data,model_init);
%% Q1
N = 200;
A1 = [1 -1.79 0.84];
C1 = [1 -0.18 -0.11];
A2 = [1 -1.79];
C2 = [1 -0.18 -0.11];
ARMApoly1 = idpoly(A1,[],C1);
ARMApoly2 = idpoly(A2,[],C2);
%pzmap(ARMApoly1)
%pzmap(ARMApoly2)
sigma2 = 1.5;
e = sqrt(sigma2) * randn(N,1);
y1 = filter(ARMApoly1.c,ARMApoly1.a,e);
y2 = filter(ARMApoly2.c,ARMApoly2.a,e); %pole outside unit circle --> divergens
subplot(211)
plot(y1)
subplot(212)
plot(y2)
%% Q2
m=80;
rtheo = kovarians(ARMApoly1.c,ARMApoly1.a,m);
stem(0:m,rtheo * sigma2)
hold on
rest = covf(y1,m+1);
stem(0:m,rest,'r')
%% Q3
% subplot(3,1,1)
% normplot(y1)
% subplot(3,1,2)
% acf(y1,20, 0.05, 1);
% subplot(3,1,3)
% pacf(y1,20, 0.05, 1);
data = iddata(y1);
ar_model = arx(data, 3); %Slightly lower than 2, maybe because of error
                       %and cost for extra parameters may not be high
                       %enough. This has the same FPE as the arma below.
arma_model = armax(data,[2 0]); %This is the best for arma, adding any
                       %any ma part will increase more because of cost of
                       %extra parameters rather then decreasing because of
                       %accuracy. Therefore adding extra parameters is not
                       %worth it according to FPE
ehat=filter(ar_model.a,ar_model.c,y1);
%%
subplot(3,1,1) %These look really good with the aproximation of
normplot(ehat) %the arma [2 0] model computed above.
subplot(3,1,2) %Alot of models looked pretty good, i.e.
acf(ehat,20, 0.05, 1);  %The residual was white noise
subplot(3,1,3)
pacf(ehat,20, 0.05, 1);
%% (3.2) Q4
n = 500;
A = [1 -1.35 0.43];
sigma2 = 4;
noise = sqrt(sigma2)*randn(n+100,1);
y = filter(1,A,noise);
y= y(101:end); %Why discard 100? dont we just need 2 previous values to
subplot(211)   %compute the firt y value of the AR-process
plot(y)
subplot(212)
plot(noise)
%%
n_est=floor(2/3*n);
y_est = iddata(y(1:n_est));
y_val = iddata(y(n_est+1:end));
NN = [1:10]';
V = arxstruc(y_est,y_val,NN);
n_order = selstruc(V,0);
n_aic = selstruc(V,'aic');
%%
n_order = zeros(1,100);
n_aic = zeros(1,100);
for i = 1:100
    noise = sqrt(sigma2)*randn(n+100,1);
    y = filter(1,A,noise);
    y= y(101:end);
    y_est = iddata(y(1:n_est));
    y_val = iddata(y(n_est+1:end));
    V = arxstruc(y_est,y_val,NN);
    n_order(i) = selstruc(V,0);
    n_aic(i) = selstruc(V,'aic');
end
subplot(211)
hist(n_order) %order seems to be 2 for this one but more spread then AIC
ylabel('order') %this seems to be f-distributed for many iterations
subplot(212)    %the n_order also has more procentage on 2 with many iter.
hist(n_aic) %Order to with very few values different from 2.
ylabel('AIC') 
%%
ar_model = arx(y,n_order(end));
ar_model.NoiseVariance
ar_model.CovarianceMatrix
present(ar_model)
%% (3.3) Q5
load data.dat
load noise.dat
data = iddata(data);
ar1_model = arx(data, 1);
ar2_model = arx(data, 2);
rar1 = resid(ar1_model, data); %plot with just resid(...) in commandwindow
rar2 = resid(ar2_model, data); %rar1 seems to have acf in 1 + the small in
                               %rar2.
                               %rar2 have just a small correlation left.
subplot(311)
plot(rar1.y)
ylabel('rar1')
subplot(312)
plot(rar2.y)
ylabel('rar2')
subplot(313)
plot(noise)
ylabl('noise')
% present(ar2_model) 
%The ar1_model estimates a1 better but the ar2 has lower FPE
%but the estimate for a1 is much worse. We think this is due to 
%ar2 trying to approximate the error while the ar1 just focus on the first 
%term, ignoring the c1 since it doesn't allow for that estimate anyway.
%ar2 has lower FPE so it's better in theory, we think (also better fit).
%% Q6
am11_model = armax(data,[1 1]); %Best one according to FPE, reasonably 
                                %since this is the order of the process.
am22_model = armax(data,[2 2]); %
resid(am11_model, data)
figure
resid(am22_model, data) %From these two and using KISS we use am11
%% (3.4) Q7
A = [ 1 -1.5 0.7];
C = [ 1 zeros(1,11) -0.5];
A12 = [ 1 zeros(1,11) -1];
A_star= conv(A,A12);
e = randn(600,1);
y = filter(C,A_star,e);
y=y(100:end);


subplot(3,1,1)
acf(y,50,0.05, 1);    %Slightly decaying Sinus curve with period 12
ylabel('ACF')
subplot(3,1,2)
pacf(y,20, 0.05, 1);  %Lots of values above the confidence   
ylabel('PACF')
subplot(3,1,3)
normplot(y)
%% Q8
y_s = filter(A12,1,y);
data = iddata(y_s);
subplot(3,1,1)
acf(y_s,20,0.05, 1);
ylabel('ACF') %Decaying sinus curve (not slow as before)
subplot(3,1,2)
pacf(y_s,20, 0.05, 1);  %AR(2)
ylabel('PACF')
subplot(3,1,3)
normplot(y_s)
%% Q9
model_init = idpoly([1 0 0],[],[]);
model2_armax = pem(data, model_init);
% model_init = idpoly([1 0],[],[]);
% model1_armax = pem(data, model_init); %AR(1) approx, not good
y_s1 = filter(model2_armax.a,model2_armax.c,y); %.c is a one
subplot(3,1,1)
acf(y_s1,20,0.05, 1); %Peak at 12
ylabel('ACF') 
subplot(3,1,2)
pacf(y_s1,20, 0.05, 1);  %Peak at 12
ylabel('PACF')
subplot(3,1,3)
normplot(y_s1)
%% Q10
model_init=idpoly([1 0 0],[],[1 zeros(1,12)]);
model_init.Structure.c.Free = [zeros(1,12) 1];
model_armax = pem(data,model_init);
y_s1 = filter(model_armax.a,model_armax.c,y_s);
data = iddata(y_s1);
subplot(3,1,1)
acf(y_s1,20,0.05, 1);
ylabel('ACF') 
subplot(3,1,2)
pacf(y_s1,20, 0.05, 1);
ylabel('PACF')
subplot(3,1,3)
normplot(y_s1) %We have white noise, parameters in model armax

