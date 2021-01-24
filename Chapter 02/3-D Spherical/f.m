function y = f(THETA,PHI)

% Please edit the field pattern here.
% Make sure to use dot operations(".*","./" and ".^") in the
% expression instead of regular operations("*","/" and "^")!

% Example 1
% y=abs(sin(THETA).*sin(PHI).^2);

% Example 2
lambda=1;
k=2*pi/lambda;
d=0.25*lambda;
N=10;
beta=-(k*d+pi/N);
kesi=k*d*cos(THETA)+beta; 
% y=abs(sin(N*kesi/2)./sin(kesi/2)/N);

% Example3
% lambda=1;
% k=2*pi/lambda;
% l=0.5*lambda;
% d=0.25*lambda;
% N=4;
% y=abs((cos(k*l/2*cos(THETA))-cos(k*l/2))./sin(THETA));

% Example 4
y=abs(sin(THETA).^1.5);





