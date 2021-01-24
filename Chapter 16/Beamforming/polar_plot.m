%**************************************************************************
% polar_plot.m
%*************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%-----------------------------------------------------------------------------

function polar_plot(p,d,w,AF_limit,E_plane)
% polar_plot(p,d,w)

k = 0;
k0 = 2*pi;

theta = linspace(0,pi/2,361);
phi   = linspace(0,2*pi,361);

[THETA,PHI] = meshgrid(theta,phi);

PSI_x = -k0*d(1)*sin(THETA).*cos(PHI);
PSI_y = k0*d(2)*sin(THETA).*sin(PHI);

AF = zeros(length(theta),length(phi));
for m = 1:p(1),
    for n = 1:p(2),
	k = k + 1;
        PSI = (m-1)*PSI_x + (n-1)*PSI_y;
        AF = AF + w(k) * exp(i*PSI);
    end;
end;

AF = E_plane .* AF;

AF = 20*log10(abs(AF)./max(max(abs(AF))));

AF(find(AF<=AF_limit)) = AF_limit;

figure(1);
[Y,X]=pol2cart(PHI,THETA);
surf(X,Y,AF(1:1:361,1:1:361));
shading('interp'); axis off; %colorbar;
view(-45,75);

% Draw x, y, and z axes
OX = X(181,181); OY = Y(271,271); s = 0.5;
set(line([OX;OX],[OY;min(min(Y))-s],[min(min(AF));min(min(AF))]),'Color','k'); %PHI=0
set(line([OX;max(max(X))+s],[OY;OY],[min(min(AF));min(min(AF))]),'Color','k'); %PHI=90
set(line([OX;OX],[OY;max(max(Y))+s],[min(min(AF));min(min(AF))]),'Color','k'); %PHI=180
set(line([OX;min(min(X))-s],[OY;OY],[min(min(AF));min(min(AF))]),'Color','k'); %PHI=270

% Label x, y, and z axes
text(min(min(OX)),min(min(Y))-s,min(min(AF)),'\phi=0\circ'  ,'FontSize',14,'Color','k'); %PHI=0
text(max(max(X))+s,min(min(OY)),min(min(AF)),'\phi=90\circ' ,'FontSize',14,'Color','k'); %PHI=90 
text(min(min(OX)),max(max(Y))+s,min(min(AF)),'\phi=180\circ','FontSize',14,'Color','k'); %PHI=180
text(min(min(X))-s,min(min(OY)),min(min(AF)),'\phi=270\circ','FontSize',14,'Color','k'); %PHI=270 

axis normal