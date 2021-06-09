% Liam Wigney , April, Jun 2021
% compare scaleddimensionless and dimensional models

clear
close all


options = odeset('abstol',1e-12,'reltol',1e-12) ; %lower if running to slow
odeSolver = @ ode15s;

%%
tSpan = [0 15];

bw0 = [90000; 6000];
r=4;
k=1e6;
alph=4e-6;
s=0.8;
gb = 130;
gw =40;
A = 0.02;
bet = 0.0002;

theta = (gb*A)/r;
eps = s/r;
sigm = (alph*k)/r;
phi=(gw*A)/r;

%%
dbwdt = @(t,bw)  [
    bw(1)-bw(1).^2-(bw(1)*bw(2))-(theta*bw(1));
    (sigm*bw(2)*bw(1))-(eps*bw(2))-(phi*bw(2))];

%tSpanD = tSpan/r;
bwd(1)=bw0(1)/k;
bwd(2)=bw0(2)/(r/bet);
[tnd, pnd] = odeSolver(dbwdt, tSpan, bwd, options); 

pn(:,1) = pnd(:,1)*k;
pn(:,2) = pnd(:,2)*r/bet;
tn = tnd/r;
%figure(1)

p = plot(tn, pn);
%set(p, {'DisplayName'}, {sprintf('initial condition %i and %i', bw0(1), bw0(2))})
hold on

%% 
dbwdt2 = @(t,bw)  [
    ((r*bw(1))*(1-(bw(1)/k)))-(bet*bw(1)*bw(2))-(A*gb*bw(1));
    (alph*bw(2)*bw(1))-(s*bw(2))-(gw*A*bw(2))];

options = odeset('abstol',1e-12,'reltol',1e-12) ; %lower if running to slow
odeSolver = @ ode15s;

[tn2, pnd2] = odeSolver(dbwdt2, tSpan, bw0, options); 

pn2(:,1) = pnd2(:,1);
pn2(:,2) = pnd2(:,2);
%figure(2)

pp = plot(tn2, pn2,'--');
%%


%%
xlabel('Solution $t$ time','interpreter','latex')
ylabel('Solution $p$ population','interpreter','latex')
legend

