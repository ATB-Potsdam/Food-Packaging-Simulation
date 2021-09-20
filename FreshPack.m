clc,clear
D=2;            %(Experiment duration, days
dur=D*24;       %(Experiment duration, hours)
t=1/3600;       %(Time interval, hours)
NS=dur/t;       %(Number of steps)

% PRODUCT PROPERTIES (related to strawberry MAP)
Wp=0.4;     %(Product mass, kg)
np=14;      %(Product number)
Cp=3722.1;      %(specific heat of fresh produce, J/(kg.K))
Vp=(1.098*Wp*1000-0.39)*10^-6;      %(Fruit volume, m3)
dp=2*(3*Vp/(4*np*pi))^(1/3);            %(Fruit average radius, m)
As=(1.862*Wp*1000+10.650)*10^-4;        %(Overall surface area of fruits, m2)

% PACKAGE PROPERTIES
Wt=0.03;        %(Packaging tray mass, kg)
x=0.001;        %(Tray wall thickness, m)
e=100*10^-6;    %(Film thickness, m)
Po2ref=19*10^-6;%(Coefficient of permeability of film to O2, mg/(m2.h))
Pco2ref=76*10^-6;%(Coefficient of permeability of film to CO2, mg/(m2.h))
WVPref=26*10^-5;%(Permeability of film to water vapour, mg/(m2.h.Pa))
Kt=0.19*3600;   %(thermal conductivity of polyethylene, j/(h.m.K))
Ct=1670;        %(specific heat of plastic, j/(kg.K))
Dp=0.001;       %(Radius of the perforation, m)
Np=1;           %(Number of film perforations)
W=0.105;        %(width package assumed as horizontal plate, m)
L=0.128;        %(length package assumed as horizontal plate, m)
D=0.07;         %(Depth of package assumed as vertical plate)
Ap1=L*W;        %(surface area of top of package, m2)
Ap1=L*W;        %(surface area of bottom of package, m2)
Ap2=D*(W+L+W+L);%(surface area of sides as symmetrical trapezoids, m2)
Ap=2*Ap1+Ap2;   %(surface area of package, m2)
A=Ap1;          %(Breathable film area, m^2)
Vt=W*L*D;       %(Total package volume, m3)
Vf=Vt-Vp;       %(Package free volume, Cubic meter)
 
% ENVIRONMENTAL CONDITIONS
Tout=20;        %(Ambient temperature, C)
RHout=50;       %(Ambient relative humidity, %)
yo2out=20.9;    %(Outside O2 concentration, percent)
yco2out=0.03;   %(Outside CO2 concentration, percent)
%% Initial values
T(1)=Tout(1);   %(Initial headspace temperature, C)
Ts(1)=Tout(1);  %(Initial fruit surface temperature, C)
Tt(1)=Tout(1);  %(Initial package tray wall temperature, C)     
RH(1)=RHout(1);     %(Initial headspace relative humidity, %)
yo2(1)=yo2out(1); %(Initial O2 concentration, percent)
yco2(1)=yco2out(1); %(Initial CO2 concentration, percent)
Mtr(1)=0;       %(Initial total moisture transpiration, mg)
Mfilm(1)=0;     %(Initial total moisture permeation, mg)
Mabsorb(1)=0;   %(Initial total moisture transpiration, mg)
Mloss(1)=0;     %(Initial total oxidative mass loss, mg)
MlossT(1)=0;    %(Initial total fruit mass loss, mg)
Det(1)=0.27;        %(Initial total microbial deterioration, %)
Det_max=13;         %(Maximum acceptable total microbial deterioration, %)

%Constants
%% Module 1: Packaging gas concentration
kmo2=2.63;             %(Michaelis-Menten constant, Percent)
vmo2ref=0.27;          %(Max O2 consumption rate at reference Temp, micromole/kg.S)
Eavmo2=74826;          %(O2 respiration activation energy, J/mol)
RQ=0.91;               %(Respiratory quotient)
kmco2f=0.056;          %(Michaelis-Menten constant, Percent)
vmco2fref=0.50;        %(Max CO2 production rate at reference Temp, micromole/kg.S)
Eavmco2f=57374;        %(CO2 respiration activation energy, J/mol)
Do2=0.073; 		%(Diffusion coefficient of O2 in air) m2/h
Dco2=0.0742; 	%(Diffusion coefficient of CO2 in air, m2/h)
Tref=10;		%(Reference temperature, C)
R=8.314472; 	%(Gas constant, J/mol.K)
%% Module 2: Packaging humidity and condensation
Patm=1.01325*10^5;%(Atmospheric Pressure, Pa)
Rd=287.058;		%(Specific gas constant for dry air, 287.058 J/(kg·K)
Rv=461.495;		%(Specific gas constant for water vapor, 461.495 J/(kg·K)
awi=0.99;  %(Water activity on fruit surface: 0-1)
Ks=13.6*10^-3*3600;  %(Fruit skin mass transfer coefficient, mg/(m2.h.Pa))
MMh2o=0.018016;	%Molar mass of water vapor, 0.018016 kg/mol)
Babsorb=11.4;
MMdair=0.028964;	%(Molar mass of dry air, 0.028964 kg/mol)
%% Module 4: Microbial deterioration and shelf life
yco2max=30;     %Maximal quantity (%) of CO2 withstanding by microorganisms
%% Module 5: Temperature
Mu=1.87e-5;     %(Dynamic viscosity, Pa.s=kg/(m.s2))
k=0.027;        %(Thermal conductivity of air, W/(m.K)=J/(s.m.K))
g=9.81;     %(gravitational acceleration, m/s2)

for i=1:NS
%% Module 1: Packaging gas concentration
% O2
Densityo2(i)=(1.429-0.0049*T(i));
Ro2(i)=yo2(i)/(kmo2+yo2(i))*vmo2ref*exp((Eavmo2/R)*(1/(Tref+273.15)-...
    1/(Ts(i)+273.15)))*(32*3.6);
dVo2p(i)=t*(Po2ref/e*A+(10^6*pi*(Dp/2)^2*Do2*Np/(e+(Dp/2))))*...
    (yo2out/100-yo2(i)/100);
dVo2r(i)=-Wp*t*Ro2(i)/Densityo2(i);
dVo2(i)=dVo2p(i)+dVo2r(i);
dyo2(i)=dVo2(i)*10^-6/Vf*100;
yo2(i+1)=yo2(i)+dyo2(i);
% CO2
Densityco2(i)=(1.977-0.0068*T(i));
Rco2(i)=(44*3.6)*(RQ*Ro2(i)/(32*3.6)+(1/(1+yo2(i)/kmco2f+1))*...
    vmco2fref*exp((Eavmco2f/R)*(1/(Tref+273.15)-1/(Ts(i)+273.15))));
dVco2p(i)=t*(Pco2ref/e*A+(10^6*pi*(Dp/2)^2*Dco2*Np/(e+(Dp/2))))*...
    (yco2out/100-yco2(i)/100);
dVco2r(i)=Wp*t*Rco2(i)/Densityco2(i);
dVco2(i)=dVco2p(i)+dVco2r(i);
dyco2(i)=dVco2(i)*10^-6/Vf*100;
yco2(i+1)=yco2(i)+dyco2(i);
%% Module 2: Packaging humidity and condensation
% Psychrometry
Lv(i)=(2502535.259-2385.76424*T(i))*10^-6;	%(Latent heat of vaporisation, j/mg)
Dair(i)=Patm/(287.05*(T(i)+273.15));	%(Density of air, kg/m3)%      
Psat(i)=610.94*exp(17.625*T(i)/(243.04+T(i)));	%(Water vapour saturation pressure in headspace, Pa)%
Psats(i)=610.94*exp(17.625*Ts(i)/(243.04+Ts(i)));	%(Water vapour saturation pressure, Pa)%
Psatout(i)=610.94*exp(17.625*Tout(1)/(243.04+Tout(1)));	%(Water vapour saturation pressure in ambient, Pa)%
Pin(1)=RH(1)/100*Psat(1);	%(Water vapour partial pressure, Pa)%
Pdair(i)=Patm-Pin(i);	%(Dry air partial pressure, Pa)%
Dairhumid(i)=Pdair(i)/(Rd*(T(i)+273.15))+Pin(i)/(Rv*(T(i)+273.15));%(Density of humid air, kg/m3)
Wa(i)=Vf*Dairhumid(i); 	%(Mass of headspace air, kg)
Mairmax(i)=0.002166*Psat(i)/(T(i)+273.16)*Vf*10^6; 	%(Max. moisture content of headspace air, kg)
Mtotal(1)=0.002166*Pin(1)/(T(1)+273.16)*Vf*10^6; 	%(initial moisture in packaging free space kg)
Mair(1)=Mtotal(1); 	%(initial moisture content of headspace air, kg)
% Moisture condensation and relative humidity
if Mtotal(i)>Mairmax(i)
    Mcond(i+1)=Mtotal(i)-Mairmax(i);Mair(i)=Mairmax(i);
else
    Mcond(i+1)=0;Mair(i)=Mtotal(i);
end
dMcond(i)=Mcond(i+1)-Mcond(i);
RH(i)=100*(Patm/Psat(i))*(Mair(i)*10^-3/(0.622*(Dairhumid(i)*1000*Vf)...
    +Mair(i)*10^-3));
% Moisture transpiration  
Dv(i)=9.1*10^-9*((T(i)+273.15)^2.5)/((T(i)+273.15)+245.18)*3600;        %(Diffusivity of water vapor in air (m2/h))
Ka(i)=2*Dv(i)/((T(i)+273.15)*dp*(R/MMh2o))*10^6*3600/1;       %(Thin air layer mass transfer coefficient, mg/(m2.h.Pa))
Ps(i)=awi*Psats(i);     %(Water vapor pressure on fruit surface, Pa)
VPD(i)=Ps(i)-Pin(i);   %(Vapor Pressure Deficit, Pa)
dMtr(i)=t*As*VPD(i)/(1/Ks+1/Ka(i));
if dMtr(i)<0
dMtr(i)=0;
end
Mtr(i+1)=Mtr(i)+dMtr(i);
% Moisture absorption   
Mabsorbeq(i)=0.057*exp(0.057*RH(i))-Mabsorb(i)/1000;
dMabsorb(i)=1*Mabsorbeq(i)*(1-exp((-t/24)/Babsorb))*1000;
Mabsorb(i+1)=Mabsorb(i)+dMabsorb(i);   
% Moisture permeation
Pin(i)=(Psat(i)*RH(i)/100);
Pout(i)=(Psatout(i)*RHout(1)/100);
Dwv(i)=Dair(i)*(1+Mair(i)*10^-3/Wa(i))/(1+MMdair/MMh2o*Mair(i)*10^-3/Wa(i));
dMfilm(i)=Dwv(i)*t*(WVPref/e*A+(10^6*pi*(Dp/2)^2*Dv(i)*Np/(e+(Dp/2))))*Dwv(i)*(RH(i)/100-RHout(1)/100)*Psat(i)/(R/MMh2o*(T(i)+273.15));         
Mfilm(i+1)=Mfilm(i)+dMfilm(i);
% Final moistures
dMtotal(i)=dMtr(i)-dMabsorb(i)-dMfilm(i);
Mtotal(i+1)=Mtotal(i)+dMtotal(i);
if Mtotal(i+1)>Mairmax(i)
    Mair(i+1)=Mairmax(i);
else
    Mair(i+1)=Mtotal(i);
end
Pin(i+1)=Mair(i+1)/(Vf*10^6)*(T(i)+273.16)/0.002166;

%% Module 3: Fruit mass loss
% Oxidative mass loss
dMloss(i)=Rco2(i)*Wp*t*(180-108)/264;    
Mloss(i+1)=Mloss(i)+dMloss(i);
% Total mass loss
MlossT(i)=Mloss(i)+Mtr(i);

%% Module 4: Microbial deterioration and shelf life
%(Deterioration, %)
Ksp(i)=(3e-6*log(T(i))+2e-6)*3600;        %(Spoilage constant  h^-1)
co2_rel(i)=1-yco2(i)/yco2max;           %(Deterioration inhibition parameter)
dDet(i)=t*Ksp(i)*Det(i)*((100-Det(i))/100)*co2_rel(i); 
Det(i+1)=Det(i)+dDet(i);
%(Keeping quality (shelf life under given conditions))
if Det(i)<=Det_max 
    KQ(i)=i;
else
    KQ(i)=0;
end
%% Module 5: Temperature
% convective heat transfer coefficient for packaging tray
Ca(i)=(1.005+1.82*Mair(i)*10^-6/(Vf*Dairhumid(i)))*10^3;            %(Humid heat of air, J/kg.K) 
beta1(i)=1/((T(i)+Tt(i))/2+273.15);
%(Grashof number for top, bottom and side wall)
Grtop(i)=L^3*Dairhumid(i)^2*g*abs(T(i)-Tt(i))*beta1(i)/Mu^2;
Grbott(i)=L^3*Dairhumid(i)^2*g*abs(T(i)-Tt(i))*beta1(i)/Mu^2;
Grside(i)=D^3*Dairhumid(i)^2*g*abs(T(i)-Tt(i))*beta1(i)/Mu^2;
Pr(i)=Mu*Ca(i)/k;       %Prandtl number
%(Rayleigh number for top, bottom and side wall)
Ratop(i)=Grtop(i)*Pr(i);
Rabott(i)=Grbott(i)*Pr(i);
Raside(i)=Grside(i)*Pr(i);
%(Nusselt number for top, bottom and side wall)
if Tt(i)<T(i)
    Nutop(i)=0.54*Ratop(i)^0.25;Nubott(i)=0.27*Rabott(i)^(1/3);
else
    Nutop(i)=0.27*Ratop(i)^(1/3);Nubott(i)=0.54*Rabott(i)^0.25;
end
Nuside(i)=(0.825+(0.387*Raside(i)^(1/6)/((1+(0.492/Pr(i))^(9/16))^(8/27))))^2;        %(For all values of Ra)
%(heat transfer coefficient (W/(m2.K)for top, bottom and side wall)
hptop(i)=3600*k/L*Nutop(i);
hpbott(i)=3600*k/L*Nubott(i);
hpside(i)=3600*k/D*Nuside(i);
hp(i)=1*(Ap1*hptop(i)+Ap1*hpbott(i)+Ap2*hpside(i))/Ap;      % (Average heat transfer coefficient (W/(m2.K))
% convective heat transfer coefficient for fruit
beta2(i)=1/((T(i)+Ts(i))/2+273.15);
Gr(i)=dp^3*Dairhumid(i)^2*g*abs(T(i)-Ts(i))*beta2(i)/Mu^2;      %Grashof number)
Pr(i)=Mu*Ca(i)/k;       %(Prandtl number)
Ra(i)=Gr(i)*Pr(i);      %(Rayleigh number)
Nu(i)=(2+0.589*Ra(i)^0.25/(1+(0.469/Pr(i))^(9/16))^(4/9));      %(Nusselt number)
hs(i)=1*3600*k*Nu(i)/dp;        %(heat transfer coefficient, W/(m2.K))
% final temperatures
%(Tray wall)
dTt(i)=(t*((Kt*Ap*(Tout(1)-Tt(i))/x)+hp(i)*Ap*(T(i)-Tt(i))+(dMabsorb(i)-dMfilm(i)+dMcond(i))*Lv(i)))/(Wt*Ct);
Tt(i+1)=Tt(i)+dTt(i);
%(Fruit surface)
dTs(i)=(t*(Wp*Rco2(i)*6.21+hs(i)*As*(T(i)-Ts(i)))-(dMtr(i))*Lv(i))/(Wp*Cp);
Ts(i+1)=Ts(i)+dTs(i);
%(Headspace air)
dT(i)=(t*(hp(i)*Ap*(Tt(i)-T(i))-hs(i)*As*(T(i)-Ts(i)))+(dMtotal(i))*Lv(i))/(Wa(i)*Ca(i));
T(i+1)=T(i)+dT(i);
    end

                            %OUTPUT VISUALIZATION
%_____________________________________________________________________________
xaxisdiv=48;xaxislable='Time (h)';

figure(1)
plot(yo2,'LineWidth',2),hold on,plot(yco2,'g--','LineWidth',2)
legend('O_2','CO_2');
xlabel(xaxislable,'fontsize',14);ylabel('Gas Concentration(%)','fontsize',14);
set(gca,'XTick',0:xaxisdiv/t:NS,'XTickLabel',{0:xaxisdiv:NS*t},'fontsize',16,'XLim',[0 NS]);

figure(2)
plot(Mtr*10^-3,'LineWidth',2),hold on,plot(Mabsorb*10^-3,'r','LineWidth',2),plot(Mcond*10^-3,'g','LineWidth',2),plot(Mfilm*10^-3,'k','LineWidth',2)
legend('Moisture loss','Absorption','Condensation','Permeation','Total mass loss','Ox. mass loss','Location','northwest')
xlabel(xaxislable,'fontsize',14);ylabel('Mass (g)','fontsize',14);
set(gca,'XTick',0:xaxisdiv/t:NS,'XTickLabel',{0:xaxisdiv:NS*t},'fontsize',16,'XLim',[0 NS],'YLim',[0 4]);

figure(3)
plot(RH,'LineWidth',2)
xlabel(xaxislable,'fontsize',14);ylabel('Relative Humidity (%)','fontsize',14);
set(gca,'XTick',0:xaxisdiv/t:NS);set(gca,'XTickLabel',{0:xaxisdiv:NS*t},'fontsize',16);
set(gca,'XTick',0:xaxisdiv/t:NS,'XTickLabel',{0:xaxisdiv:NS*t},'fontsize',16,'XLim',[0 NS]);

figure(4)
plot(MlossT*10^-3/(Wp*10^3)*100,'LineWidth',2),hold on,plot(Mloss*10^-3/(Wp*10^3)*100,'r','LineWidth',2),plot(Mtr*10^-3/(Wp*10^3)*100,'g','LineWidth',2)
legend('Total mass loss','Ox. mass loss','Moisture loss','Location','northwest')
xlabel(xaxislable,'fontsize',14);ylabel('Loss (%)','fontsize',14);
set(gca,'XTick',0:xaxisdiv/t:NS,'XTickLabel',{0:xaxisdiv:NS*t},'fontsize',16,'XLim',[0 NS],'YLim',[0 1]);

figure(5)
plot(Det,'LineWidth',2),hold on
xlabel(xaxislable,'fontsize',14);ylabel('Deterioration (%)','fontsize',14);
MaxKQ=max(KQ);
line([0 NS],[Det_max Det_max],'color',[0 0 0],'LineStyle','-','LineWidth',0.5)
line([MaxKQ MaxKQ],[0 Det_max],'color',[0 0 0],'LineStyle','--','LineWidth',0.5)
annotation(figure(5),'textbox',[0.19 0.67 0.28 0.05],'String','MAD','LineStyle','none','FontSize',14,'FitBoxToText','off');
annotation(figure(5),'textbox',[0.67 0.24 0.16 0.05],'String','Shelf life','LineStyle','none','FontSize',14,'FitBoxToText','off');
set(gca,'XTick',0:xaxisdiv/t:NS,'XTickLabel',{0:xaxisdiv:NS*t},'fontsize',16,'XLim',[0 NS],'YLim',[0 20]);

figure(6)
plot(T,'LineWidth',2),hold on,plot(Ts,'r','LineWidth',2),plot(Tt,'g','LineWidth',2),plot(Tout,'k','LineWidth',2)
legend('T_{in}','T_{fruit surface}','T_{tray wall}','T_{out}','Location','northwest')
xlabel(xaxislable,'fontsize',14);ylabel('Temperature (°C)','fontsize',14);
set(gca,'XTick',0:xaxisdiv/t:NS);set(gca,'XTickLabel',{0:xaxisdiv:NS*t},'fontsize',16);
set(gca,'XTick',0:xaxisdiv/t:NS,'XTickLabel',{0:xaxisdiv:NS*t},'fontsize',16,'XLim',[0 NS],'YLim',[Tout(1) Tout(1)+2]);
