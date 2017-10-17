function xdot = f(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file                Author: 
%
% Project: Simulation of a hybrid system (bouncing ball)
%
% Description: Flow map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% global data -----------
global A varsigma

% states
N = length(A);
xp      = x(1:N);               % state of the consensus variable
xo1     = x(N+1:2*N);
xo2     = x(2*N+1:3*N);
eta1    = x(3*N+1:4*N);
eta2    = x(4*N+1:5*N);
% timer  = x(end-1:end);

xpdot  = A*xp;
xo1dot = A*xo1 + eta1;
xo2dot = A*xo2 + eta2;
eta1dot = zeros(2,1);
eta2dot = zeros(2,1);
timerdot = -ones(N, 1) + varsigma*ones(N, 1);

xdot = [xpdot;...
        xo1dot;...
        xo2dot;...
        eta1dot;...
        eta2dot;...
        timerdot]; 
    
end