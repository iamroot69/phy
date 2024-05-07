clc;clear;clf;
h=6.626D-34;c=3D8;k=8.617333262D-5;v=1;s=1;u=-1;s1=0.5;u1=2
T=[100,1000,5000]
e=linspace(0,6,1000)  //Mev
e1=linspace(0,6,1000)
function y=n(e,u,T,A)
    k=8.617333262D-5
    for i=1:length(e)
        y=1./(exp((e-u)./(k.*T))+A)
    end
endfunction
c=(2*s)*(4*%pi*v)/((h^3)*c^3)
g=c*(e.^2);
subplot(2,3,1)
plot(e,g)
subplot(2,3,2)
plot(e,n(e,u,1D3,-1),e,n(e,u,1D4,-1))
subplot(2,3,3)
f=g'*n(e,u,1D4,-1)
f1=g'*n(e,u,2D4,-1)
plot(e',diag(f),e',diag(f1))
c1=(2*s1)*(4*%pi*v)/(h^3*c^3)
g1=c1*(e.^2)
subplot(2,3,4)
plot(e,g)
subplot(2,3,5)
plot(e1,n(e1,u1,1D3,1),e1,n(e1,u1,2D1,1),'--')
subplot(2,3,6)
f2=g1'*n(e1,u1,1D3,1)
f3=g1'*n(e1,u1,2D1,1)
plot(e1',diag(f2),e1',diag(f3),'--')
