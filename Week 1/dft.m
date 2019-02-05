function dft = dft(x,N)
% Calulates the N-Point DFT of a vector x
dft = dftmtx(N)*x;