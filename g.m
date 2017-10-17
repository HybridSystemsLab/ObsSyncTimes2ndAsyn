function xplus = g(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file               Author: Torstein Ingebrigtsen Bø
%
% Project: Simulation of a hybrid system (bouncing ball)
%
% Description: Jump map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global G A H1 H2 K11 K12 K21 K22 T1 T2 W mm

% states
N = length(A);
xp      = x(1:N);               % state of the consensus variable
xo1     = x(N+1:2*N);
xo2     = x(2*N+1:3*N);
eta1    = x(3*N+1:4*N);
eta2    = x(4*N+1:5*N);
timer  = x(end-1:end);


xpplus  = xp; 
xo1plus = xo1;
xo2plus = xo2;
eta1plus = eta1;
eta2plus = eta2;
timerplus = timer;



m1 = mm*rand;
m2 = mm*rand;
c1x = mm*rand*ones(2,1);
c2x = mm*rand*ones(2,1);
c1y = mm*rand;
c2y = mm*rand;

% measurement
yy1 = H1*xp + m1;
yy2 = H2*xp + m2;

xo2t = xo2 + c1x;
xo1t = xo1 + c2x;
yy2t = yy2 + c1y;
yy1t = yy1 + c2y;

% eta1plus = K11*(H1*xo1 - yy1) + K12*(H2*xo2t - yy2t) + W*(xo1 - xo2t);
% eta2plus = K21*(H1*xo1t - yy1t) + K22*(H2*xo2 - yy2) + W*(xo2 - xo1t);

theta = 0.9;

dropoutrate = 0;
if rand <= dropoutrate
    drop = 0;
else
    drop = 1;
end

if timer(1) <= 0
%     eta1plus = K11*sign(H1*xo1 - yy1).*abs(H1*xo1 - yy1).^theta ...
%            + K12*sign(H2*xo2t - yy2t).*abs(H2*xo2t - yy2t).^theta ...
%            + W*sign(xo1 - xo2t).*abs(xo1 - xo2t).^theta;
    eta1plus = K11*(H1*xo1 - yy1) + drop*K12*(H2*xo2t - yy2t) + drop*W*(xo1 - xo2t);
    timerplus(1) = rand*(T2 - T1) + T1;
end       
     
if timer(2) <= 0
%     eta2plus = K21*sign(H1*xo1t - yy1t).*abs(H1*xo1t - yy1t).^theta ...
%            + K22*sign(H2*xo2 - yy2).*abs(H2*xo2 - yy2).^theta ...
%            + W*sign(xo2 - xo1t).*abs(xo2 - xo1t).^theta;
    eta2plus = K21*(H1*xo1t - yy1t) + drop*K22*(H2*xo2 - yy2) + drop*W*(xo2 - xo1t);
    timerplus(2) = rand*(T2 - T1) + T1;   
end

% timerplus = ones(2,1)*(rand*(T2 - T1) + T1);

xplus = [xpplus;...
         xo1plus;...
         xo2plus;...
         eta1plus;...
         eta2plus;...
         timerplus];
             
end