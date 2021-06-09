clear
syms theta xi phi sigm
j1 = [1-theta 0; 0 -xi-phi]
eig(j1)


syms b w
b=(xi+phi)/sigm
w = (-b-theta+1)


%syms gb A r alph k s r gw

%theta = (gb*A)/r;
%xi = s/r;
%sigm = (alph*k)/r;
%phi=(gw*A)/r;


j2=[1-(2*b)-w-theta, -b; sigm*(w), (sigm*b)-xi-theta]
eig(j2)



det(j2)
trace(j2)

r=4;
k=1e6;
alph=4e-6;
s=0.8;
gb = 130;
gw =20;
A = 0.02;