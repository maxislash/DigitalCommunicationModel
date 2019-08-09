clear all, close all, clc

nIterations = 2e3;
Fs          = 1e6;   % Nominal sampling frequency
y_ppm       = 50;    % Fractional frequency offset in ppm
Fc_tx       = 1e3;   % Nominal clock frequency
pn_var      = 1e-9;  % Phase noise variance
Kp          = 0.05;  % Proportional Constant
Ki          = 0.01;  % Integral Constant

nIterations = 100;
Fs = 10000;
Fc_tx = 400;
Kp = 0.5;

%% Derived parameters
y        = y_ppm * 1e-6;  % Fractional frequency offset in ppm
f_offset = y * Fc_tx;     % Absolute freq error (in Hz)
F_offset = f_offset/Fs;   % Normalized absolute frequency offset

% Nominal Phase Increment (used in the loop phase accumulator)
delta_phi_c = 2*pi*Fc_tx/Fs;

% PLL pull range
% Normalize frequency offset must be smaller than Kp*pi, namely:
%   F_offset < Kp*pi
if (F_offset > Kp*pi)
    warning('Frequency offset is not within the pull range\n');
end

delta_phi_in = 2*pi*(Fc_tx + f_offset)/Fs;
phase_noise  = sqrt(pn_var)*randn(nIterations, 1);
phi_in       = (0:nIterations-1).' * delta_phi_in + phase_noise;
% Finally, generate the exponential:
s_in         = exp(1j * phi_in);

NCO           = zeros(nIterations, 1);
% Initialize the loop DDS to a random phase
NCO(1)  = sqrt(pi)*randn;
Vi = 0;
V = zeros(nIterations, 1);
s_loop             = zeros(nIterations, 1); 
dds_mult           = zeros(nIterations, 1);
phi_error          = zeros(nIterations, 1);

for i = 1:nIterations

  %% Loop DDS:
  s_loop(i) = exp(1j*NCO(i));
  dds_mult(i) = s_in(i) * conj(s_loop(i));
  phi_error(i) = imag(dds_mult(i));

  %% Loop Filter
  Vp = phi_error(i)*Kp;
  Vi = Vi + phi_error(i)*Ki;
  V(i) = Vp + Vi;

  %% NCO
  NCO(i+1) = NCO(i) + delta_phi_c + V(i);

end

%% Performance
% Input vs. Output Instantaneous Phase
figure
plot(phi_in)
hold on
plot(NCO, 'r')
title('Instantaneous Phase')
xlabel('Sample')
ylabel('Phase (rad)')
legend('Input', 'Output')

% Filtered Phase Error
phi_error_ss_expected = 0;
phi_error_filtered_ss_expected = 2*pi * F_offset;
figure
plot(phi_error)
hold on
plot(V, 'g')
hold on
plot(phi_error_ss_expected * ones(nIterations,1), '-.k')
hold on
plot(phi_error_filtered_ss_expected *  ones(nIterations,1), '--r')
title('Loop Filter Output')
legend('Phase Error', ...
    'Filtered Phase Error', ...
    'Expected Phase Error Steady-State', ...
    'Expected Filtered Steady-State')
xlabel('Sample')
ylabel('Error (rad)')

% Input vs. Output Sinusoids
figure
subplot(121)
plot(real(s_in))
hold on
plot(real(s_loop), 'r')
ylabel('Amplitude')
xlabel('Sample')
title('Real part')
legend('Input', 'Output')

subplot(122)
plot(imag(s_in))
hold on
plot(imag(s_loop), 'r')
ylabel('Amplitude')
xlabel('Sample')
title('Imaginary part')
legend('Input', 'Output')
