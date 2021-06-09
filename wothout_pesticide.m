clear
syms eps sigm



syms b w
b=(eps)/sigm
w = -b+1


syms gb A r alph k s r gw

%eps = s/r;
%sigm = (alph*k)/r;


j2=[1-(2*b)-w, -b; sigm*(w), (sigm*b)-eps]
eig(j2)



dd=det(j2)
tt=trace(j2)

r=4;
k=1e6;
alph=4e-6;
s=0.8;
gb = 130;
gw =40;
A = 0.02;

(eps*s*sigm - 2*alph*eps^2*k + alph*eps*k*sigm)/(r*sigm^2)
(alph*eps*k)/(r*sigm) - s/r - eps/sigm

