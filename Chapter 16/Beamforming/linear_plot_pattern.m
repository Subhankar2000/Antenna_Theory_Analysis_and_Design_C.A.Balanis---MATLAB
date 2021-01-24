%************************************************************************
%	Linear_plot_pattern(w,N,d,E_pattern,AF_limit,TYPE,rticks,line_style)
%************************************************************************
%	PLOT_PATTERN is a MATLAB function that plots 2-D patterns in
%	polar and rectangular coordinates where:
%
%	Input Parameters Description
%	----------------------------
%	- w          weight coefficients (row or column vector)
%	- N          number of elements
%	- d          inter-element spacing (in wavelength)
%                default value is 0.5 (Nyquist rate)
%	- E_pattern  samples of element pattern (column vector)
%                if interested in array factor only, enter 1
%                default value is 1
%	- AF_limit   lower radial tick limit for polar plot (in dB)
%                default value is -40
%	- TYPE       options are:
%			     'half' plots from -90 to 90 degrees
%			     'full' plots from -180 to 180 degrees
%                default value is 'half'   
%   - rticks     is the # of radial ticks (or circles) desired.
%	             default value is 4
%   - line_style is solid (e.g., '-') or dashed (e.g., '--')
%                default value is '-'                 
%************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%-----------------------------------------------------------------------------

function linear_plot_pattern(sig,w,N,d,E_pattern,AF_limit,TYPE,rticks,line_style)

k0 = 2*pi;

%%%%%%%%%%%%%%% Parameters initialization %%%%%%%%%%%%%%%
%---- default values ----%
def_d = .5;
def_E_pattern = 1;
def_AF_limit = -40;
def_TYPE = 'half';
def_rticks = 4;
def_line_style = '-';
H_FIG = 1;
%------------------------%

switch nargin
case 2
    d = def_d; E_pattern = def_E_pattern; AF_limit = def_AF_limit;
    TYPE = def_TYPE; rticks = def_rticks; line_style = def_line_style;
case 3
    E_pattern = def_E_pattern; AF_limit = def_AF_limit;
    TYPE = def_TYPE; rticks = def_rticks; line_style = def_line_style;
case 4
    AF_limit = def_AF_limit;
    TYPE = def_TYPE; rticks = def_rticks; line_style = def_line_style;
case 5
    TYPE = def_TYPE; rticks = def_rticks; line_style = def_line_style;
case 6
    rticks = def_rticks; line_style = def_line_style;
case 7
    line_style = def_line_style;
end

m = linspace(0,N-1,N);
switch TYPE
case 'half'
    theta_1 = -pi/2; theta_2 = pi/2;
    theta = linspace(-pi/2,pi/2,181);
case 'full'
    theta_1 = -pi; theta_2 = pi;
    theta = linspace(-pi,pi,361);
otherwise
    disp('ERROR => Only two options are allowed for TYPE: half or full');
    return
end 
%-------------------------------------------------------%
    
%%%%%%%%%%%%%% Total Pattern [Element * AF] %%%%%%%%%%%%%
AF = E_pattern .* sum(diag(w)*exp(i*(k0*d*m'*sin(theta))));
AF = 20*log10(abs(AF)./max(abs(AF)));
%-------------------------------------------------------%
    
%%%%%%%%%%%%%%%%%%%%%%%% Polar Plot %%%%%%%%%%%%%%%%%%%%%
figure(H_FIG);
switch TYPE
case 'half'
    linear_semipolar_dB(theta*180/pi,AF,AF_limit,0,rticks,line_style);
case 'full'
    linear_polar_dB(theta*180/pi,AF,AF_limit,0,rticks,line_style);
end

title(['Linear Array Beamforming Pattern ( N=', num2str(N),', d=', num2str(d), '\lambda, SOI \theta=', num2str(sig(1,2)), '\circ)']);

%-------------------------------------------------------%
    
%%%%%%%%%%%%%%%%% Rectangularr Plot %%%%%%%%%%%%%%%%%%%%%
H_FIG = H_FIG + 1; figure(H_FIG);
plot(theta*180/pi,AF,line_style);
axis([theta_1*180/pi theta_2*180/pi AF_limit 0]); grid on;
xlabel('\theta (degrees)','fontsize',14);
ylabel('Magnitude (dB)','fontsize',14);
a=sprintf('\theta (degrees)','fontsize',14);
fs  = get(gca, 'FontSize'); set(gca, 'FontSize', 14);
title(['Linear Array Beamforming Pattern ( N=', num2str(N),', d=', num2str(d), '\lambda, SOI \theta=', num2str(sig(1,2)), '\circ)']);

%-------------------------------------------------------%

