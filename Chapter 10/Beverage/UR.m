function y = UR(TH,PHI)

global h l K

y =  ( sin( 2*pi*h.*sin(TH).*cos(PHI) ).^2 ).* ...
     (sin(TH).^2).*( cos(pi*l.*cos(TH)) - cos(pi*l*K )  ).^2./(K^2-cos(TH).^2).^2 ;
