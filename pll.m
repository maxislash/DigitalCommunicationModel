clear all;
close all;

Fs = 10000;
Fc_tx = 400;
N = 100;

NCO(1)  = sqrt(pi)*randn;
Vi = 0;

%% To fine tune
Kp = 0.6;
Ki = 0.01;

f_offset = 0.02;     % Absolute freq error (in Hz)
delta_phi_in = 2*pi*(Fc_tx + f_offset);
phase_tx = 0;
pn_var      = 1e-9;
phase_noise  = sqrt(pn_var)*randn(N, 1);
phase = phase_tx + phase_noise;


Tx_I_BB = zeros(1, N);
dt = 1/Fs; 
t = 0 : dt : (length(Tx_I_BB) - 1) * dt;
in       = t.' * delta_phi_in + phase;
s_in         = exp(1j * in);


Tx_I = 1 .* (sqrt(2)*cos(2*pi*Fc_tx*t+phase));
Tx_Q = 1 .* (sqrt(2)*sin(2*pi*Fc_tx*t+phase));
Tx = Tx_I - Tx_Q;
s_in = Tx;

for k = 1 : N
  
  % Phase error detector
##  e = sin(in(k) - NCO(k)); % phase error between input and nco
##  err(k) = e;
  s_loop(k) = exp(1j*NCO(k));
  dds_mult(k) = s_in(k) * conj(s_loop(k));
  err(k) = imag(dds_mult(k));
  
  % Loop filter
  Vp = Kp * err(k);
  Vi = Vi + Ki * err(k);
  V(k) = Vp + Vi;
  
  % Numerically Controlled Oscillator
  NCO(k+1) = NCO(k) + 2*pi*Fc_tx/Fs + V(k);
  
endfor

figure(1)
plot(err);

figure(2)
plot(sin(in));
hold on
plot(sin(NCO));
legend("input", "output");
