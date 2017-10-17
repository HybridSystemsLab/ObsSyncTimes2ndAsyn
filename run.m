%% simulation of event based static consensus

%%% global data -----------
clear all
global G A H1 H2 K11 K12 K21 K22 T1 T2 W mm varsigma

% plant information 
A  = [0 -1;1 0];
% H1 = [1 0]; 
H1 = [0 0];
H2 = [0 1];
T1 = 0.1; T2 = 0.3;
K11 = -[0.5 0.2]';
K12 = -[0.2 0.2]';
K21 = [0.2 0.3]';
K22 = [-0.1 -0.5]';
W   = -0.4;

% skew on the dynamics of clock
varsigma = 0.0;

% mm = 0.8; %%% noise factor
mm = 0;

%%%% Graph (1) - 2 agents
G = ones(2,2);

%%% -----------------------
% IC for plant states
xp0 = [2 2]';

% IC for agent1;
xo10 = [15 5]';
eta10 = [1 1]';
timer10 = 0;

% IC for agent2;
xo20 = [-1 0]';
eta20 = [-1 -1]';
timer20 = 0.3;

y0 = [xp0; xo10; xo20; eta10; eta20; timer10; timer20]; 

% simulation horizon
TSPAN = [0 10];
JSPAN = [0 20000];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-1,'MaxStep',1e-2);

% simulate
[t y j] = hybridsolver(@f,@g,@C,@D,y0,TSPAN,JSPAN,rule,options,1);

save('y1.mat','y') 
save('t1.mat','t')

%%
figure
hold on 
plot(y(:,1),y(:,2),'r','linewidth',3)
plot(y(:,3),y(:,4),'k--','linewidth',1.5)
plot(y(:,5),y(:,6),'b-.','linewidth',1.5)
grid on
set(gca,'FontSize',20)
legend('x0', 'x1', 'x2')
plot(y(1,1),y(1,2),'*r','linewidth',2)
plot(y(1,3),y(1,4),'*k','linewidth',1.5)
plot(y(1,5),y(1,6),'*b','linewidth',1.5)
axis equal
box on 
%%

%%
figure 
subplot(311)
plot(t, y(:,3) - y(:,1), 'k--','linewidth',1.5)
hold on
plot(t, y(:,5) - y(:,1), 'b-.','linewidth',1.5)
grid on 
set(gca,'FontSize',20)
legend('e11', 'e12')
axis([0, TSPAN(2), -20, 20])
subplot(312)
plot(t, y(:,4) - y(:,2), 'k--','linewidth',1.5)
hold on
plot(t, y(:,6) - y(:,2), 'b-.','linewidth',1.5)
grid on
set(gca,'FontSize',20)
legend('e21', 'e22')
axis([0, TSPAN(2), -20, 20])
subplot(313), hold on, grid on, 
plotHarcT(t, j, y(:,11), 'k')
plotHarcT(t, j, y(:,12), 'b')
axis([0, 4, 0, 0.4]), set(gca,'FontSize',20)
h1 = plot([0,0],[1,3],'k','linewidth',2);
hold on
h2 = plot([0 3],[1,1],'b','linewidth',2);
legend([h1 h2],{'phi1','phi2'})
set(h1,'Visible','off'); set(h2,'Visible','off');
%%

