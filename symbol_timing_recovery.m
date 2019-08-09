newSymbol = 0;
CNT_next = 1;
mu_next = 0;
Vi = 0;
receivedSymbols_I = zeros(1, numberOfSymbols + 1);
receivedSymbols_Q = zeros(1, numberOfSymbols + 1);
n = 2;
%K     ->  TEQ Gain
% K0     ->  Counter (interpolator controller) gain
% eta    ->  Damping factor
% Bn_Ts  ->  PLL Noise Bandwith normalized by the symbol rate (multiplied by the symbol period Ts)
test = 0;
Kp = 0;
Ki = 0;

%% PID constant
K = (g(Ns/2 + 1) - g(Ns/2 - 1))/(2/Ns);
% Convert (Bn*Ts), i.e. multiplied by symbol period, to (Bn*T), i.e.
% multiplied by the sampling period.
Bn_times_T = Bn_Ts / Ns;

% Theta_n (a definition in terms of T and Bn)
Theta_n = (Bn_times_T)/(eta + (1/(4*eta)));

% Constants obtained by analogy to the continuous time transfer function:
K_K0_K1 = (4 * eta * Theta_n) / (1 + 2*eta*Theta_n + Theta_n^2);
K_K0_K2 = (4 * Theta_n^2) / (1 + 2*eta*Theta_n + Theta_n^2);

% K1 and K2 (PI constants):
Kp = Kp_K0_K1/(K*K0);
Ki = Kp_K0_K2/(K*K0);

for i = Ns/2 : length(Rx_I_BB)
  
  CNT = CNT_next;
  mu = mu_next;
  
  if newSymbol == 1
    
      k = i - 1;
      
      %% Interpolator
      receivedSymbols_I(n) = mu * Rx_I_BB(k + 1) + (1 - mu) * Rx_I_BB(k);
      receivedSymbols_Q(n) = mu * Rx_Q_BB(k + 1) + (1 - mu) * Rx_Q_BB(k);
      
      %% Gardner TED
      e = (receivedSymbols_I(n) - receivedSymbols_I(n-1)) *  test;
      e = e + (receivedSymbols_Q(n) - receivedSymbols_Q(n-1)) *  test;
      
      err(n) = e; %for plotting
      
      n = n + 1;
  
  else
      e = 0;
      
  endif
  
  %% Loop filter
  Vp = Kp * e; %Proportional
  Vi = Vi + Ki * e; %Integral
  V = Vp + Vi;

  %% Modulo-1 Counter
  W = 1/Ns + V;
  CNT_next = mod(CNT - W, 1);
 
  if CNT < W
    newSymbol = 1;
    mu_next = CNT/W;
  else
    newSymbol = 0;
    mu_next = mu;
  endif
endfor

receivedSymbols = complex(receivedSymbols_I, receivedSymbols_Q);