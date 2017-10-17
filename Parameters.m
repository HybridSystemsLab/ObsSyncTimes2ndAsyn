%%% finding feasible parameters for oscillatory plant
clear all 
clc

A = [0 1;-1 0];
H1 = [1 0];
H2 = [1 0.4];

%%% gain matrices 

%%% aux params
delta = -10.52;
% beta = 0.2;
% eps = 2000;

%%% timer bounds 
T1 = 0.0002; T2 = 0.0002;
nu = [T1, T2];

N = length(A);
ZE = zeros(N,N);
IN = eye(N);

Af = [A, ZE, IN;...
      ZE, A, IN;...
      ZE, ZE, ZE];
Bf = [ZE, ZE, ZE;...
      ZE, ZE, -IN;...
      ZE, ZE, ZE];

%%% optimization variable
K1 = sdpvar(2,1,'full');
K2 = sdpvar(2,1,'full');
% P  = sdpvar(3*N,3*N);
beta = sdpvar(1,1);

Ag1 = [IN+K1*H1, ZE, ZE;...
       K1*H1, IN, ZE;...
       ZE, K1*H1, ZE];
Ag2 = [IN+K1*H1, ZE, ZE;...
       K2*H2, IN, ZE;...
       ZE, K2*H2, ZE];

A1 = exp(delta*nu(1)/2)*expm(Af*nu(1))*Ag1;
A2 = exp(delta*nu(2)/2)*expm(Af*nu(2))*Ag2;


% F = [P >= 1e1*eye(size(Af))];
F = [A1'*A1 - eye(size(Af)) <= -beta*eye(size(Af))];
F = [F; A2'*A2 - eye(size(Af)) <= -beta*eye(size(Af))];

% obj = trace(PP1);
obj = -beta;
yalpen = export(F,obj,sdpsettings('solver','penbmi'),[],[],1);

bmi = yalmip2bmi(yalpen);
penm = bmi_define(bmi);
prob = penlab(penm);
solve(prob);
% Pvec = prob.x;
% Pvar = zeros(N);
% inddd = find(tril(ones(N)));
% Pvar(inddd) = Pvec;
% P = Pvar + tril(Pvar,-1)';

K1 = prob.x(1:2);
K2 = prob.x(3:4); 
beta = prob.x(5);

beta

% P = 1e2*eye(size(Af));

Ag1 = [IN, ZE, ZE;ZE, IN, ZE;K1*H1, K1*H1, ZE];
Ag2 = [IN, ZE, ZE;ZE, IN, ZE;K2*H2, K2*H2, ZE];

A1 = exp(delta*nu(1)/2)*expm(Af*nu(1))*Ag1;
A2 = exp(delta*nu(2)/2)*expm(Af*nu(2))*Ag2;

%%% test whether the P satisfies the condition 
flag = 0;
% for ii = nu(1):0.01:nu(2)
%     A1i = expm(Af*ii)*Ag1;
%     A2i = expm(Af*ii)*Ag2;
%     if (find(eig(exp(delta*ii)*A1i'*A1i - eye(size(Af)) + beta*eye(size(Af))) > 0))
%         flag = 1;
%         
%     elseif (find(eig(exp(delta*ii)*A2i'*A2i - eye(size(Af)) + beta*eye(size(Af))) > 0))
%         flag = 1;
%         ii
%     end
% end

display(['flag = ' num2str(flag)])

betamax = 1000;

if (flag == 0) 

alpha1 = 1e20;
alpha2 = 0;

for ii = 0:.001:nu(2)
    tau = ii;
%     alpha_eig = abs(eig(exp(delta*tau)*expm(Af'*tau)*P*expm(Af*tau)));
    alpha_eig = eig(exp(delta*tau)*expm(Af'*tau)*expm(Af*tau));
    if alpha1 > min(alpha_eig);
        alpha1 = min(alpha_eig);
    end
    if alpha2 < max(alpha_eig)
        alpha2 = max(alpha_eig);
    end
end
        
rho = 0;

for ii = 0:.01:nu(2)
    rho_tem = norm(exp(delta*ii)*Bf'*expm(Af'*ii)*expm(Af*ii)*Bf,2);
    if rho < rho_tem
        rho = rho_tem;
        iter = ii;
    end
end

display(['rho = ' num2str(rho)])
display(['alpha1 = ' num2str(alpha1)])
display(['alpha2 = ' num2str(alpha2)])

% rhsout = exp((1/eps - delta + eps*rho/alpha1)*nu(2));
rhsout = exp(2*sqrt(rho/alpha1)*nu(2) - delta*nu(2));
lhsout = 1/(1 - beta);
        
if lhsout > rhsout
    fprintf('Stability Condition Satisfied: LHS > RHS; %3.2f > %3.2f\n ',lhsout,rhsout )
else
    fprintf('Stability Condition NOT Satisfied (LHS >= RHS); %3.2f !=> %3.2f\n ',lhsout,rhsout )
end

end