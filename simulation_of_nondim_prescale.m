%% Liam Wigney Jun 2021

clear
close all


options = odeset('abstol',1e-12,'reltol',1e-12) ; %lower if running to slow
odeSolver = @ ode15s;

%%
tSpan = [0 150000];

bw0 = [60000; 9000];
r=4;
k=1e6;
alph=4e-6;
s=0.8;
gb = 0;
gw =400;
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
pn(:,2) = pnd(:,2);
tn = tnd;
%figure(1)

p1 = plot(tn, pn(:,1));
hold on
xlabel('Solution $t$ time','interpreter','latex')
ylabel('Solution $p$ population','interpreter','latex')
set(p1, {'DisplayName'}, {sprintf('Bees when A=%.3f', A)})
legend

p2 = plot(tn, pn(:,2));
hold on
xlabel('Solution $t$ time','interpreter','latex')
ylabel('Solution $p$ population','interpreter','latex')
set(p2, {'DisplayName'}, {sprintf('Wasps when A=%.3f', A)})
legend


%%
A=0

dbwdt = @(t,bw)  [
    bw(1)-bw(1).^2-(bw(1)*bw(2))-(theta*bw(1));
    (sigm*bw(2)*bw(1))-(eps*bw(2))-(phi*bw(2))];

%tSpanD = tSpan/r;
bwd(1)=bw0(1)/k;
bwd(2)=bw0(2)/(r/bet);
[tnd4, pnd4] = odeSolver(dbwdt, tSpan, bwd, options); 

pn4(:,1) = pnd4(:,1;
pn4(:,2) = pnd4(:,2);
tn4 = tnd4;
%figure(1)

p14 = plot(tn4, pn4(:,1));
hold on
xlabel('Solution $t$ time','interpreter','latex')
ylabel('Solution $p$ population','interpreter','latex')
set(p14, {'DisplayName'}, {sprintf('Bees when A=%.3f', A)})
legend

p24 = plot(tn4, pn4(:,2));
hold on
xlabel('Solution $t$ time','interpreter','latex')
ylabel('Solution $p$ population','interpreter','latex')
set(p24, {'DisplayName'}, {sprintf('Wasps when A=%.3f', A)})
legend