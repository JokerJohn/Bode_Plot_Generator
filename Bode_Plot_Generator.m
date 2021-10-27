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
for omega = 10e-2:stepsize:10e2
    i=i+1;
    s = omega*1i;
    
    % **** ADD YOUR TRANSFER FUNCTION HERE ****
    numer_TF = 12250*(s+4);
    denom_TF = (s^2)*((s^2)+28*s+4900);
    % **** ADD YOUR TRANSFER FUNCTION HERE ****
    
    % Calculations for plots
    omega_val(i) = omega;
    numer_re = real(numer_TF);
    numer_im = imag(numer_TF);
    denom_re = real(denom_TF);
    denom_im = imag(denom_TF);
    FRF_mag = (sqrt((numer_re)^2+(numer_im)^2)) / (sqrt((denom_re)^2+(denom_im)^2));
    DB_mag(i) = 20*log10(abs(FRF_mag));
    numer_phase_rad = atan2(numer_im,numer_re);
    denom_phase_rad = atan2(denom_im,denom_re);
    FRF_phase_rad(i) = numer_phase_rad - denom_phase_rad;
end

FRF_phase_deg = (360/(2*pi)) * unwrap(FRF_phase_rad);

if max(FRF_phase_deg) > 90
    FRF_phase_deg = FRF_phase_deg - 180;
elseif min(FRF_phase_deg) < -180
    FRF_phase_deg = FRF_phase_deg + 180;
end

% Create magnitude and phase plots
figure
set(gcf, 'color', 'w')

subplot(2,1,1)
plot(omega_val,DB_mag,'k-')
set(gca, 'XScale','log');
xlabel('Frequency [rad/s]')
ylabel('Magnitude [dB]')
xlim([10e-2 10e2])
ymin = min(DB_mag)-20;
ymax = max(DB_mag)+20;
ylim([ymin ymax])
xticks([10e-2 10e-1 10e0 10e1 10e2]);
yticks([-1000:20:1000]);
grid

subplot(2,1,2)
plot(omega_val,FRF_phase_deg,'k-')
set(gca, 'XScale','log');
xlabel('Frequency [rad/s]')
ylabel('Phase [deg]')
xlim([10e-2 10e2])
ymin = min(FRF_phase_deg)-22.5;
ymax = max(FRF_phase_deg)+22.5;
ylim([ymin ymax])
xticks([10e-2 10e-1 10e0 10e1 10e2]);
yticks([-3600:45:3600]);
grid