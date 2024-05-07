clc;clear;clf();
kb=1.38e-23;Na=6.022e23;
Elmnt=["Nitrogen";"Oxygen";"Flourine";"Neon"]
m1=[28;32;19;20]
m=m1/(Na*1000);
T=400:400:1600;
v=0:1:3500;
v2=v.^2;iter=1;
style=['k-d';'k-o';'k-x';'k-']

subplot(1,2,1)
mprintf("\t   |\t\t\t\t\t     |\t\t\t\t\t     |\n\t   |\t\t From Formula\t\t     |\t\tFrom distribution graph      |\n\t   |\t\t\t\t\t     |\t\t\t\t\t     |\n")
mprintf("\t   |\t\t\t\t\t     |\t\t\t\t\t     |\n")
mprintf("----------------------------------------------------------------------------------------------\n")
mprintf("\tElement : "+string(Elmnt(1))+"\t\t\t\t\t\t\t\t     |\n")
mprintf("----------------------------------------------------------------------------------------------\n")
mprintf("Temperature|    Vmp     |Vav          |Vrms          |Vmp         |Vav          |Vrms\t     |\n")
for i=1:length(T)
     for j=1:length(v)
          P_m(i,j)=4*%pi*v2(j)*((m(1)/(2*%pi*kb*T(i)))^1.5)*exp(-0.5*m(1)*v2(j)/(kb*T(i)))
 
     end
     [d,k]=max(P_m(i,:))
     vmp_m(i)=sqrt(2*kb*T(i)/m(1)); vav_m(i)=sqrt(8*kb*T(i)/(m(1)*%pi)); vrms_m(i)=sqrt(3*kb*T(i)/m(1))
     Vmp_m(i)=v(k);
     
     
      Vav_m(i)=(4*P_m(i,1:2:$)*v(1:2:$)'+2*P_m(i,2:2:$)*v(2:2:$)')/3
     Vrms_m(i)=sqrt((4*P_m(i,1:2:$)*v2(1:2:$)'+2*P_m(i,2:2:$)*v2(2:2:$)')/3)
     mprintf("%8d K | %5.2f m/s | %7.2f m/s | %7.2f m/s  | %5.2f m/s | %7.2f m/s | %7.2f m/s|\n",T(i),vmp_m(i),vav_m(i),vrms_m(i),Vmp_m(i),Vav_m(i),Vrms_m(i))
     mprintf("\t   |\t\t|\t      |\t\t     |\t          |\t        |\t     |\n")
     plot(v(1:10:$),P_m(i,1:10:$),style(i),'linewidth',2)
     replot([0 2500;0 %nan])
     l(i)=["T="+string(T(i))]
end
title("Maxwell velocity distribution for "+Elmnt(1)+"(mol mass="+string(m1(1))+" g/mol)",'fontsize',6)
xlabel("Velocity(m/s)",'fontsize',5)
ylabel("P(v)",'fontsize',5)
L=legend(l)
L.font_size=4
subplot(1,2,2)
mprintf("-------------------------------------------------------------------------------------------------------------\n")
mprintf("Element    |    Mass    |    Vmp      |Vav           |Vrms         |Vmp          |Vav         |Vrms\t    |\n")
mprintf("-------------------------------------------------------------------------------------------------------------\n")
for i=1:length(m1)
     for j=1:length(v)
          P_t(i,j)=4*%pi*v2(j)*((m(i)/(2*%pi*kb*T(1)))^1.5)*exp(-0.5*m(i)*v2(j)/(kb*T(1)))
     end
     [d,k]=max(P_t(i,:))
     vmp_t(i)=sqrt(2*kb*T(1)/m(i)); vav_t(i)=sqrt(8*kb*T(1)/(m(i)*%pi)); vrms_t(i)=sqrt(3*kb*T(1)/m(i))
     Vmp_t(i)=v(k); Vav_t(i)=(4*P_t(i,2:2:$-1)*v(2:2:$-1)'+2*P_t(i,3:2:$-1)*v(3:2:$-1)')/3
     Vrms_t(i)=sqrt((4*P_t(i,2:2:$-1)*v2(2:2:$-1)'+2*P_t(i,3:2:$-1)*v2(3:2:$-1)')/3)
     mprintf("%8s   |  %2d g/mol  | %5.2f m/s  | %7.2f m/s  |%7.2f m/s  | %5.2f m/s | %7.2f m/s | %7.2f m/s |\n",Elmnt(i),m1(i),vmp_t(i),vav_t(i),vrms_t(i),Vmp_t(i),Vav_t(i),Vrms_t(i))
     mprintf("\t   |\t\t|\t      |\t\t     |\t           |\t        |\t      |\t\t    |\n")
     plot(v(1:10:$),P_t(i,1:10:$),style(i),'linewidth',2)
     replot([0 2000;0 %nan])
     l(i)=["Element:"+Elmnt(i)+", mol mass="+string(m1(i))+" g/mol"]
end
mprintf("-------------------------------------------------------------------------------------------------------------\n")
title("Maxwell velocity distribution for T="+string(T(1))+"K",'fontsize',6)
xlabel("Velocity(m/s)",'fontsize',5)
ylabel("P(v)",'fontsize',5)
L=legend(l)
L.font_size=4
msprintf([Elmnt' vmp_m' Vmp_m'])
