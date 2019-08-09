step = 1;
e = zeros(1, numberOfSymbols);
k(1) = 2 * Ns * G;
k(2) = k(1) + Ns;

for i = 2 : numberOfSymbols
  
  receivedSymbols_I(i-1) = Rx_I_BB(k(i));
  receivedSymbols_Q(i-1) = Rx_Q_BB(k(i));
  
  e(i) = (Rx_I_BB(k(i)) - Rx_I_BB (k(i-1))) * Rx_I_BB(k(i) - Ns/2);
  e(i) = e(i) + (Rx_Q_BB(k(i)) - Rx_Q_BB (k(i-1))) * Rx_Q_BB(k(i) - Ns/2);
  
  if mean(e) > 0
    delta(i) = -step;
  elseif mean(e) < 0
    delta(i) = step;
  else
    delta(i) = 0;
  endif
  
  k(i+1) = k(i) + Ns + delta(i);
  k_plot(i+1) = mod(k(i+1) +Ns/2, Ns);
  
endfor

receivedSymbolsT = complex(receivedSymbols_I, receivedSymbols_Q);
##figure(21);
##stem(2*Ns*G :length(Rx_I_BB), Rx_I_BB(2*Ns*G: end));
##hold on
##stem( k( 1 : end-1), Rx_I_BB( k( 1: end-1) ), 'filled');
##
##figure(22);
##stem(Rx_I_BB(k(1:end-1)), 'filled');

figure(23);
stem(k_plot);

figure(24)
scatter(real(receivedSymbolsT), imag(receivedSymbolsT) , "filled");
axis([-2 2 -2 2])
title("Received Symbols T")