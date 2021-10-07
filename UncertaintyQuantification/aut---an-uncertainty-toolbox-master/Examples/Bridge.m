% Example for uncertainty calculation for a Wheatstone Bridge
% using compensation Method

Ug=0;           % Voltage shown at balance-instrument
R2=1e3;         % nominal resistor value
R3=1e3;         % nominal resistor value
R4=1.001e3;     % nominal resistor value to balance the bridge
U0=5;           % nominal supply voltage for the bridge

u_R2 = 5;       % standard uncertainty of resistor
u_R3 = 5;       % standard uncertainty of resistor
u_R4 = 1;       % standard uncertainty of resistor

u_U0 = 0.1;     % standard uncertainty of supply
u_Ug = 0.01;    % standard uncertainty of balance instrument


% assign the values to the variables
R2 = unc(R2,u_R2,'R2');
R3 = unc(R3,u_R3,'R3');
R4 = unc(R4,u_R4,'R4');
U0 = unc(U0,u_U0,'U0');
Ug = unc(Ug,u_Ug,'Ug');

% actual calculations
I2=U0/(R3+R4);          % current through R3 and R4
I1=(I2*R4+Ug)/R2;       % current through R2 (and thus R1)
R1 = U0/I1-R2           % calculate result for R1
disp_contribution(R1)  % show uncertainty contributions to R1

