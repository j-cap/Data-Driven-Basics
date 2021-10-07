% Example on the calculation of an Inverting Amplifier with OPA
% Linear System - Superposition!

U1 = [0 1e-3]           % Input Signal
R1=unc(1000,10,'R1')    % Resistors for amplification V=-10
R2=unc(10000,100,'R2')  % Standard uncertainty 1/100 R
B = 1000                % equivalent noise bandwidth
T = 300                 % Temperature

% random, yet static deviations from ideal OPA

IBm=unc(10e-12,3e-12,'IBm')
IBp=unc(10e-12,3e-12,'IBp')
Uofs=unc(0,1e-4,'Uofs')

k=unc(1.38064852e-23,0.00000079e-23,'k');   % Boltzmann Konstante dzt. SI System!

% Noise - Varies for each measurement!
Un_t=unc([0 0],10e-9*sqrt(B)*[1 1],{'Un','Un'})
UnR1_t=unc([0 0],sqrt(gmv(4*k*T*R1*B))*[1 1],{'UnR1','UnR1'})
UnR2_t=unc([0 0],sqrt(gmv(4*k*T*R2*B))*[1 1],{'UnR2','UNR2'})

Inm_t=unc([0 0],10e-11*sqrt(B)*[1 1],{'Inm','Inm'})
Inp_t=unc([0 0],10e-11*sqrt(B)*[1 1],{'Inp','Inm'})


% Contributions on the outpur of the amplifier

UA_U1=-U1*R2/R1;
UA_Uofs=Uofs*(1+R2/R1);
UA_IBm=-IBm*R2;
UA_IBp=0;

UA_Un=Un_t*(1+R2/R1)
UA_UnR1=-UnR1_t*R2/R1
UA_UnR2=-UnR2_t;
UA_Inm=-Inm_t*R2;
UA_Inp=0;


% Superposition of results (linear system!)

UA=UA_U1+UA_Uofs+UA_IBm+UA_IBp+UA_Un+UA_UnR1+UA_UnR2+UA_Inm+UA_Inp
disp_contribution(UA(2))

% difference method (aka: correlated double sampling, autozeroing, ...)

UAd=UA(2)-UA(1)
disp_contribution(UAd)

