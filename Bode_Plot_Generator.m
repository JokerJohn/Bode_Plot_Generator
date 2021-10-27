% Bode Plot Generator
% Created by Evan Schimmel on 10/26/2021
% Free for editing and distribution
% If you encounter issues, email schimmer@rose-hulman.edu

close all
clear all
clc

% Specify stepsize, 0.1 is plenty small
stepsize = 0.1;
i=0;

% Loop over frequencies
for omega = 0.1:stepsize:1000
    i=i+1;
    s = omega*1i;
    
    % **** ADD YOUR TRANSFER FUNCTION HERE ****
    TF_numer = s+10;
    TF_denom = s^2+0.6*s+1;
    % **** ADD YOUR TRANSFER FUNCTION HERE ****
    
    % Calculations for plots
    freq(i) = omega;
    re_numer = real(TF_numer);
    im_numer = imag(TF_numer);
    re_denom = real(TF_denom);
    im_denom = imag(TF_denom);
    FRF_mag = (sqrt((re_numer)^2+(im_numer)^2)) / (sqrt((re_denom)^2+(im_denom)^2));
    DB_mag(i) = 20*log10(abs(FRF_mag));
    numer_phase = atand(im_numer/re_numer);
    denom_phase = atand(im_denom/re_denom);
    if im_numer > 0 && re_numer < 0
        numer_phase = 180 + numer_phase;
    elseif im_denom > 0 && re_denom < 0
        denom_phase = 180 + denom_phase;
    elseif im_numer < 0 && re_numer < 0
        numer_phase = 180 + numer_phase;
    elseif im_denom < 0 && re_denom < 0
        denom_phase = 180 + denom_phase;
    end
    FRF_phase(i) = numer_phase - denom_phase;
end

% Create magnitude and phase plots
figure
set(gcf, 'color', 'w')
subplot(2,1,1)
plot(freq,DB_mag,'k-')
set(gca, 'XScale','log');
xlabel('Frequency [rad/s]')
ylabel('Magnitude [dB]')
axis([0.1 1000 -40 40])
xticks([10e-1 10e0 10e1 10e2 10e3]);
yticks([-40 -30 -20 -10 0 10 20 30 40]);
grid
subplot(2,1,2)
plot(freq,FRF_phase,'k-')
set(gca, 'XScale','log');
xlabel('Frequency [rad/s]')
ylabel('Phase [deg]')
axis([0.1 1000 -180 180])
xticks([10e-1 10e0 10e1 10e2 10e3]);
yticks([-180 -135 -90 -45 0 45 90 135 180]);
grid