% Example for uncertainty calculation for a Maxwell-Wien Bridge.
% We compute the value of the impedance Zx to find the unknown 
% inductance (Lx) and resistance (Rx) using compensation Method

Ug = unc(0,0.01,'Ug');              %  Voltage shown at balance-instrument
R1 = unc(1.001e3,1,'R1');           %  resistor value - arm 1
R2 = unc(1.001e3,3,'R2');           %  resistor value - arm 2
R3 = unc(1.001e3,3,'R3');           %  resistor value - arm 3
C1 = unc(15.e-6,1e-6,'C1');         %  capacitor value - arm 1
Urms = unc(10,5,'Urms');            %  Voltage (RMS) of the power supply
f = 50;                             %  linear frequency of the power supply
w = 2*pi*f;                         %  angular frequency of the power supply

% equivalent impedance of arm 1 (Z1) formed by R1||C1. Z1 = R1 /( 1+j*w*R1*C1 ) 
realZ1 = R1 /(1+(w*R1*C1)^2);         % real part of the impedance Z1
imagZ1 = - w*C1*R1^2/(1+(w*R1*C1)^2); % img part of the impedance Z1 
Z1 = complex(realZ1, imagZ1);         % defining the impedance Z1

Z2 = R2;                            % equivalent impedance of arm 2 formed by R2
Z3 = R3;                            % equivalent impedance of arm 3 formed by R3


% actual caluclations
I1rms = Urms/(Z1+Z2);               % current through Z1 and Z2
I2rms = (I1rms*Z1+Ug)/Z3;           % current through Z3 (and thus Zx)

% calculate the impedance Zx, this arm is formed by the resistor Rx and the
% inductor Lx connected in series. This measurement model considers 
% the influences introduced by the choice of Urms and the accuracy 
% of the instrument that measures Ug.  
Zx = ((Urms*Z3*(Z1+Z2))/(Urms*Z1 + Ug*(Z1+Z2))) - Z3                                    

Rx = real(Zx)                       % resistor Rx is the real part of Zx
X_L = imag(Zx)                      % reactance XL is the img part of Zx
Lx = X_L/w   

disp_contribution(Rx)
disp_contribution(Lx)

