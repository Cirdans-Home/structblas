%% Test for circulant operation

n = 10;
c = zeros(n,1);
c(1) = 1;
c(n) = -1;
C = gallery("circul",c);
x = rand(n,1);

y = C'*x;

y2 = ifft(conj(fft(C(:,1))).*fft(x));

norm(y-y2)