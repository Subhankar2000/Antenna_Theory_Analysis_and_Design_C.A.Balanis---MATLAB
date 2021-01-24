function y = f(theta)

% Example 1
y = abs(sin(theta).^1.5);
%Example 2
a1=4.9886;
a2=3.0851;
a3=1;
d=0.7124;
u=pi*d*cos(theta);
AF=a1*cos(u)+a2*cos(3*u)+a3*cos(5*u);
E=cos(pi*cos(theta)/2)./sin(theta);
y=abs(E.*AF);
y=abs(AF);
