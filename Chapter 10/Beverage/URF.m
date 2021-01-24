function y = URF(TH)

global l K

y = (sin(TH).^2).*( cos(pi*l.*cos(TH)) - cos(pi*l*K )  ).^2./(K^2-cos(TH).^2).^2 ;
