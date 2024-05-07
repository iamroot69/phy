clc;clear;clf();
h=6.626e-34; kb=1.38e-23; c=3e8
temp=1200:100:1700; L0=0.5:0.002:30
l=L0*1e-6; temp0=1500;
a=8*%pi*h*c;
style=['b-';'kx-';'g:';'c-.';'m--d';'ro-']
figure(0);scf(0);clf(0);
//subplot(1,3,1)
for i=1:length(temp)
    for j=1:length(L0)
        Ep(i,j)=(a/(l(j)^5))/(exp(h*c/(l(j)*kb*temp(i)))-1)
        Ew(i,j)=(a/(l(j)^5))/(exp(h*c/(l(j)*kb*temp(i))))
        Erj(i,j)=(a/(l(j)^4))*(kb*temp(i)/(h*c))
    end
    [m,k]=max(Ep(i,:));
    plot(L0(1:50:$-1),Ep(i,1:50:$-1),style(i),'linewidth',2)
    plegend(i)=["T="+string(temp(i))]
    Lmax(i)=L0(k);
    if temp(i)==temp0 then q=i end
end
replot([0.5 max(L0)/3;0 %nan])
xlabel("$\lambda\;(\mu m)$",'fontsize',5)
ylabel("$E(\lambda)\;(\times10^{3}Wm^{-2}\mu m^{-1})$",'fontsize',5)
title("Planck''s Radiation law",'fontsize',6)
L=legend(plegend,1)
figure(1);scf(1);clf(1);
plot(L0',[Ep(q,:)' Ew(q,:)' Erj(q,:)'],'linewidth',2)
xlabel("$\lambda\;(\mu m)$",'fontsize',5)
ylabel("$E(\lambda)\;(\times10^{3}Wm^{-2}\mu m^{-1})$",'fontsize',5)
title("Planck''s, Wien''s, Reyleigh-Jeans Distribution(T="+string(temp0)+"K)",'fontsize',6)
L=legend("Planck''s law","Wien''s law","Reyleigh-Jeans law")
L.font_size=4;
replot([0 0.33*max(L0);0 max(Ep(q,:))*1.05])
disp("Temperature(T)  Lambda_max(um)  Lambda_max*T(um K)")
disp([temp' Lmax Lmax.*temp' ])
figure(2);scf(2);clf(2);
plot(1./temp,Lmax','ro-','linewidth',2)
xlabel("$\frac{1}{T}\;(K^{-1})$",'fontsize',5)
ylabel("$\lambda_{max}\;(\mu m)$",'fontsize',5)
title("Plot of max Wavelength vs Temperature inverse",'fontsize',6)
L=legend("max Wavelength in Planck''s plot",4)
