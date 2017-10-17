function inside= D(x) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file                Author: Torstein Ingebrigtsen Bø
%
% Description: Jump set
% Return 0 if outside of D, and 1 inside of D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% global data -----------

% states
% xp      = x(1:2);               % state of plant
timer      = x(end-1:end);               % state of timer

if find(timer <= 0)
    inside = 1;
else
    inside = 0; %jump map is empty, never jumps
end