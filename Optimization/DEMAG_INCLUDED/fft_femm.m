function [ah,bh]=fft_femm(x)

x2=x(:,2);

% remember to drop the last point of x2 if equal
% to the first one before passing it to fft

N = length(x2); % number of samples
if x2(N)==x2(1)
    x2=x2(1:N-1);
end
X = 2/N*fft(x2); % make the complex FFT
X = X(1:floor(N/2)); % drop the last half
X(1) = X(1)/2; % mean value correction
ah = real(X); % cos harmonic coefficients
bh = -imag(X); % sin harmonic coefficients
end