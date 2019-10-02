%% Computer Exercise 3
%Prepp 1:

p11 = 7/8; p22 = 7/8; 
p12 = 1/8; p21 = 1/8; 
states = zeros(1,1000);
state = 1;
n = 1;
while(n<1001) 
    if (rand(1) > 7/8)
        if (state == 1) 
             state = 0; 
        else
            state = 1; 
        end
    end
    states(1,n) = state; 
    n = n+1; 
end
%% Prepp 2:
% Example of Kalman filter

% Simulate process
y = ?;

% Length of process
N = length(y);

% Set parameters
A = [1 0;
     0 1];
V1 = [1e-4 0;
      0 1e-4]; % State noise variance
V2 = 1.25; % Measure variance
%usually C should be set here to, but in this case C is a function of
%time.

% Set initial values
Vtt = var(1)*eye(2); % Initial variance
Ztt = [? ?]'; % Initial state

% Vector to store values in
Zsave=zeros(2,N);

% Kalman filter. Start from k=3, since we need old values of y.
for k=3:N
  % C is a function of time.
  C = [? ?];

  % Time update
  Ztt_1 = ?;
  Vtt_1 = ?;

  % Measure update
  Vtt = ?;
  Kt = ?;
  Ztt = ?;

  % Save the state
  Zsave(:,k)=Ztt;
end;


