%************************************************************************
% planar_plot_pattern(N,d,w,E_pattern,AF_limit,METHOD,PHI_SCAN)
%************************************************************************
%	PLOT_PATTERN is a MATLAB function that plots 3-D patterns in
%	rectangular coordinates x axis represents phi and y reresents theta
%************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%-----------------------------------------------------------------------------
function AF=planar_plot_pattern(N,d,w,E_pattern,AF_limit,METHOD,PHI_SCAN)

H_FIG = 1;

k0 = 2*pi;

switch METHOD
case 'linear'
    %%%%%%%%%%%%%%% Parameters initialization %%%%%%%%%%%%%%%
    m     = linspace(0,N-1,N);
    theta = linspace(-pi/2,pi/2,181);
    %-------------------------------------------------------%
    
    %%%%%%%%%%%%%% Total Pattern [Element * AF] %%%%%%%%%%%%%
    AF = E_pattern .* sum(diag(w)*exp(i*(k0*d*m'*sin(theta))));
    AF = 20*log10(abs(AF)./max(abs(AF)));
    %-------------------------------------------------------%
    
    figure(H_FIG);
    semipolar(theta,AF,AF_limit,0,4,'-'); grid on;
    
    AF(find(AF <= AF_limit)) = AF_limit;

    H_FIG = H_FIG + 1; figure(H_FIG);
    plot(theta*180/pi,AF); axis([-90 90 AF_limit 0]); grid on;
    xlabel('\theta (degrees)','fontsize',14);
    ylabel('Magnitude (dB)','fontsize',14);
    fs  = get(gca, 'FontSize');
    set(gca, 'FontSize', 14);
    
case 'planar'
    %%%%%%%%%%%%%%% Parameters initialization %%%%%%%%%%%%%%%
    theta = linspace(0,pi/2,361);
    phi   = linspace(0,2*pi,361);
    AF    = zeros(length(theta),length(phi));
    k     = 0;
    %-------------------------------------------------------%
    
    [THETA,PHI] = meshgrid(theta,phi);
    
    PSI_x = k0*d(1)*sin(THETA).*cos(PHI);
    PSI_y = k0*d(2)*sin(THETA).*sin(PHI);
    
    for m = 1:N(1),
        for n = 1:N(2),
            k = k + 1;
            PSI = (m-1)*PSI_x + (n-1)*PSI_y;
            AF  = AF + w(k) * exp(i*PSI);
        end;
    end;
    
    %%%%%%%%%%%%%% Total Pattern [Element * AF] %%%%%%%%%%%%%
    AF = E_pattern .* AF;
    
    AF = 20*log10(abs(AF)./max(max(abs(AF)))); AF = AF';
    %-------------------------------------------------------%

    H_FIG = H_FIG + 1; figure(H_FIG);                % H_FIG = figure(2)
    contourf(180/pi*phi,180/pi*theta,AF,50); hold on; shading('faceted');
    ylabel('Elevation');
    xlabel('Azimuth');
    caxis([AF_limit 0])
    colorbar;
    grid on;
    fs  = get(gca, 'FontSize');
    set(gca, 'FontSize', 14);
 
    if nargin == 7
        H_FIG = H_FIG + 1;
        plot_cut(AF,AF_limit,PHI_SCAN,H_FIG);
    end;

end;

