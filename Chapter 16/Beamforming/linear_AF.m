%************************************************************************
%	Z = Linear_AF(N,d,w,angles_0,E_pattern)
%************************************************************************
%	It computes the array factor and returns the level of the nulls (in dB)
%
%	Input Parameters Description
%	----------------------------
%	- N          number of elements
%	- d          inter-element spacing (in wavelength)
%                default value is 0.5 (Nyquist rate)
%	- w          weight coefficients (row or column vector)
%   - angle_0    direction of the nulls
%	- E_pattern  samples of element pattern (column vector)
%                if interested in array factor only, enter 1
%                default value is 1
%
%	Output Parameters Description
%	-----------------------------
%   - Z          returns the level of the nulls (in dB)
%*************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%-----------------------------------------------------------------------------

function Z = linear_AF(N,d,w,angles_0,E_pattern)

k0 = 2*pi;

m     = linspace(0,N-1,N);
theta = linspace(-pi/2,pi/2,181);
    
AF = E_pattern .* sum(diag(w)*exp(i*(k0*d*m'*sin(theta))));
AF = 20*log10(abs(AF)./max(abs(AF)));
    
Z = AF(91+angles_0(2:length(angles_0)));
