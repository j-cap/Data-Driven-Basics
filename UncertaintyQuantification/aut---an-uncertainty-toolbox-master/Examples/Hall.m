% Simple example for calculation of the uncertainty of a theoretical Hall
% Sensor
% 
% Please refer to 
%
% Edward Ramdsen: "Hall-Effect Sensors - Theory and Applications",
% Elsevier, Oxford, UK, 2006
%
% for more details on the model.

% create some simulated measurement data
t=0:pi/20:4*pi;
Un=(sin(t)+2)*0.01;
Ugem=unc(Un,Un.*0+0.0001,'ADC');

% assign parameters and corresponding standard uncertainties
I=unc(1e-3,1e-5,'Bias Current');
U_oIV=unc(0,1e-3,'Amplifier Offset');
U_oH=unc(0,1e-3,'Hall Sensor Offset');
V=unc(100,1,'IV Amplification');
S=unc(2e-3,1e-4,'Hall Sensor Sensitvity');
B_K=unc(0,3,'Calibration Flux Density');


% calculate the Flux Density from the measurements
% just as if there were no uncertainty...
B=(Ugem/V-U_oIV-U_oH)/I/S;


% plot results
figure;
plot(gmv(B));
hold on;
plot(gmu(B),'g-.');
legend('Magn. Flux Density','Standard Uncertainty')
ylabel('B [\muT]');
xlabel('Sample');
title('Uncalibrated Results');

% show contributions to the uncertainty of the 10th sample
legend('Magn. Flux Density','Standard Uncertainty')
disp([13,13,13,'Uncertainty contributions to results without calibration']);
disp_contribution(B(10))

% now apply calibration
BK=B(:)-B(1)+B_K;

% plot results
figure
plot(gmv(BK));
hold on;
plot(gmu(BK),'g-.');
ylabel('BK [\muT]');
xlabel('Sample');
title('Calibrated Results');

% show contributions to the uncertainty of the 10th sample
disp([13,13,13,'Uncertainty contributions to results after calibration']);
disp_contribution(BK(10))

