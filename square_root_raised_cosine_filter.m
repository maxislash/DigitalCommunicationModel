clear all;
close all;
%----------------------------------------------------
%---------------- Initial Parameters ----------------
%----------------------------------------------------

n = 24;                              %number of bits
M = 4;                              %4-QAM
k = log2(M);                        %number of bits in one symbol
m = n/k;                            %number of symbols
Fs = 1000;                          %sampling frequency
dt = 1/Fs;                          %simulation time differential
Rs = 250;                            %Symbol rate
Ns = Fs/Rs;                         %Number of samples per symbol

dataIn = randi(2,1,n) - 1;  % Generate vector of binary data

%----------------------------------------------------
%----------------- QAM Mapper -----------------------
%----------------------------------------------------

numberOfSymbols = length(dataIn)/k; %number of symbols in the FIFO

% Define mapping table applying Gray mapping

mappingTable(1) =  1 + 1*i;
mappingTable(2) = -1 + 1*i;
mappingTable(3) =  1 - 1*i;
mappingTable(4) = -1 - 1*i;

mappedSymbols = zeros(1, numberOfSymbols);

% Map bits to symbols
for i = 1:k:length(dataIn)   
  symbolBits = dataIn(i:i+k-1);
  symbolIndex = 2^1 * symbolBits(1) + 2^0 * symbolBits(2);
  
   % Mapping
   mappedSymbols((i - 1)/k + 1) = mappingTable( symbolIndex + 1);   
endfor

%----------------------------------------------------
%----------------- Pulse shaping --------------------
%----------------------------------------------------

for i = 1:length(mappedSymbols)
    mappedSymbols_upsampled( (i-1)*Ns+1 : i*Ns) = [mappedSymbols(i) zeros(1, Ns-1)];
end

A = 1;

% Sqrt-Raised-Cosine Impulse Response Filter
alpha = 0.35;     % excess bandwidth factor = 0.35 for IS 136
L = Ns; % samples per symbol = Fs/Rs
G = 4;

n = -L*G : L*G;
for i = 1 : length(n)
    if(n(i) == 0)
      B(i) = (1/L) * ( 1 - alpha + 4 * (alpha/pi));
      
    elseif (abs(n(i)) == L/(4*alpha))
      B(i) = (1 + (2/pi)) * sin(pi / (4*alpha)) + (1 - (2/pi)) * cos(pi / (4*alpha));
      B(i) = (1/L) * (alpha/sqrt(2)) * B(i);
      
    else
      B(i) = sin(pi * (1-alpha) * n(i)/L) + 4*alpha*n(i)/L * cos(pi * (1+alpha) * n(i)/L);
      B(i) = B(i) / (pi * n(i)/L * (1 - (4*alpha*n(i)/L)^2) );
      B(i) = (1/L) * B(i);
    end 
end

% Normalize by the energy
B = B/sqrt(B*B');

mappedSymbols_upsampled_and_padded = [real(mappedSymbols_upsampled) zeros(1, L*G)];
Tx_I_BB = filter(B, A, real(mappedSymbols_upsampled_and_padded)); 
Rx_I_BB = filter(B, A, Tx_I_BB); 
%----------------------------------------------------
%----------------- Figures --------------------------
%----------------------------------------------------
figure(1);
stem(real(mappedSymbols_upsampled));
figure(2);
plot(Tx_I_BB);
hold on
stem(Tx_I_BB);
hold on
stem(1:4:length(Tx_I_BB), Tx_I_BB(1:4:end), 'filled');

figure(3);
plot(Rx_I_BB);
hold on
stem(Rx_I_BB);
hold on
stem(1:4:length(Rx_I_BB), Rx_I_BB(1:4:end), 'filled');