% this program synthesizes a line source using the Woodward-Lawson method

% Written by Marios Gkatzianas	
% Arizona State University, September 2002


% ATTENTION: fill in str and d

% [theta_even,SF_rec_even,theta_samp1,b1,z1,I1]=wood(len,sf_mode,filename,pstring,2*M,Ntheta,theta,SF);
% [theta_odd,SF_rec_odd,theta_samp2,b2,z2,I2]=wood(len,sf_mode,filename,pstring,2*M+1,Ntheta,theta,SF);

AF=SF;   % fix this

[theta_even,SF_rec_even,theta_samp1,b1,z1,I1]=wood('s',0.5,len,sf_mode,filename,pstring,2*M,Ntheta,theta,SF);
[theta_odd,SF_rec_odd,theta_samp2,b2,z2,I2]=wood('s',0.5,len,sf_mode,filename,pstring,2*M+1,Ntheta,theta,SF);


% Figure 1
% ********
plot(theta,abs(SF)); hold on; grid; plot(theta_even,abs(SF_rec_even),'r'); 
plot(theta_odd,abs(SF_rec_odd),'g');

legend('Desired',[num2str(2*M),' samples'],[num2str(2*M+1),' samples']);
h1=plot(theta_samp1,abs(b1),'kd');  h2=plot(theta_samp2,abs(b2),'ks');
set([h1 h2],{'markerfacecolor','markersize'},{'k' 7;'k' 6});
set(gca,'fontsize',12,'xlim',[0 180]);
xlabel('\theta (in \circ)','fontsize',12); 
ylabel('| SF |','fontsize',12);
title('SF synthesis using the Woodward-Lawson method','fontsize',12);


% Figure 2
% ********
figure(2);
plot(z1/len,abs(I1),'r'); grid; hold on; plot(z2/len,abs(I2),'g');
set(gca,'fontsize',12);
xlabel('z/length','fontsize',12); ylabel('| I(z) |','fontsize',12);
title('Amplitude distribution for the current line source','fontsize',12);
legend([num2str(2*M),' samples'],[num2str(2*M+1),' samples']);

