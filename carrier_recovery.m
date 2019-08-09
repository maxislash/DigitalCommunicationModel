clear all

fs   = 180000;   % Sample at 180kHz
%fs   = 48000;   % Sample at 180kHz
fc1  = 49000;    % Carrier at 50kHz (with some error)
fc2  = 51000;    % Carrier at 50kHz (with some error)
fsym = 6000;     % Symbol rate = 6000 symbols/sec
fLPF = 15000;    % 12kHz signal bandwidth
T    = fs/fsym;  % Samples per symbol
N    = 500;      % Number of symbols to simulate
An   = 0.1;

fc   = [fc1*ones(T*N/2, 1); fc2*ones(T*N/2, 1)];
t    = [0:T*N-1]';
time = t ./ fs;
bi   = sign(rand(N,1) - 0.5);
bq   = sign(rand(N,1) - 0.5);
bio  = reshape((bi(1:N)*ones(1,T))', N*T, 1);
bqo  = reshape((bq(1:N)*ones(1,T))', N*T, 1);
y    = bio .* cos(2.*pi.*(fc./fs) .* t) + ...
       bqo .* sin(2.*pi.*(fc./fs) .* t) + ...
       An * randn(T*N, 1);

[B,A]   = butter(2, fLPF/(fs/2));
[B2,A2] = butter(2, 30/(fs/2));

Ts    = 1/fs;
Ts = 1/20000000;
%
% Loop filter
%
Kv  = 50000;
z   = 0.5;
wn  = 3000*2*pi;
t1  = Kv/wn^2;
t2  = 2*z/wn - 1/Kv;
Ba(1) = (Ts+2*t2)/(Ts+2*t1);
Ba(2) = (Ts-2*t2)/(Ts+2*t1);
Aa(1) = 1.0;
Aa(2) = (Ts-2*t1)/(Ts+2*t1);

ri    = zeros(T*N,1);
rq    = zeros(T*N,1);
Ri    = zeros(T*N,1);
Rq    = zeros(T*N,1);
Fo    = zeros(T*N,1);
mult  = zeros(T*N,1);
vctrl = zeros(T*N,1);
theta = zeros(T*N,1);

pha = 2*pi*fc(2)/fs;
theta(2)=pha;
for i=3:T*N
   ri(i) = y(i) * cos(pha);
   rq(i) = y(i) * sin(pha);
   theta(i) = pha;

   % LPF I & Q arms
   Ri(i) = B(1)*ri(i) + B(2)*ri(i-1) + B(3)*ri(i-2) ...
                      - A(2)*Ri(i-1) - A(3)*Ri(i-2);
   Rq(i) = B(1)*rq(i) + B(2)*rq(i-1) + B(3)*rq(i-2) ...
                      - A(2)*Rq(i-1) - A(3)*Rq(i-2);

   % VCO
   mult(i) = sign(Rq(i)) * Ri(i) - sign(Ri(i)) * Rq(i);
   vctrl(i) = Ba(1)*mult(i) + Ba(2)*mult(i-1) - Aa(2)*vctrl(i-1);
   fo = 2*pi*50000 + Kv*vctrl(i);
   pha = rem(pha + fo/fs, 2*pi);
   Fo(i) = fo/2/pi;
end

%
% Demodulate data
%
ci = conv(Ri, ones(T, 1));
cq = conv(Rq, ones(T, 1));
di = sign(ci(T:T:length(ci)));
dq = sign(cq(T:T:length(cq)));

i=sqrt(-1);
a = diff(angle(bi+i*bq));
b = diff(angle(di+i*dq));
l = length(a);
bita = [];
bitb = [];

for j=1:l
   if b(j) == 0
      bitb = [bitb; 0; 0];
   elseif (b(j) == pi/2) | (b(j) == -3*pi/2)
      bitb = [bitb; 0; 1];
   elseif (b(j) == pi) | (b(j) == -pi)
      bitb = [bitb; 1; 1];
   elseif (b(j) == 3*pi/2) | (b(j) == -pi/2)
      bitb = [bitb; 1; 0];
   end
   if a(j) == 0
      bita = [bita; 0; 0];
   elseif (a(j) == pi/2) | (a(j) == -3*pi/2)
      bita = [bita; 0; 1];
   elseif (a(j) == pi) | (a(j) == -pi)
      bita = [bita; 1; 1];
   elseif (a(j) == 3*pi/2) | (a(j) == -pi/2)
      bita = [bita; 1; 0];
   end
end

%
% Demodulate with perfect information
%
xi = y .* cos(2.*pi.*(fc./fs) .* t);
xq = y .* sin(2.*pi.*(fc./fs) .* t);
Xi = filter(B, A, xi);
Xq = filter(B, A, xq);
Ci = conv(Xi, ones(T, 1));
Cq = conv(Xq, ones(T, 1));
Di = sign(Ci(T:T:length(Ci)));
Dq = sign(Cq(T:T:length(Cq)));

plot(Fo)