function [stnice]=cosm(x);

stnice=num2str(x,'%8.4f');
stnice(stnice=='i')='j';
stnice=cellstr(stnice);                      
           
for ik=1:length(stnice),
   pp=findstr(stnice{ik},'+');  
   if ~isempty(pp),
      stnice{ik}=[stnice{ik}(1:pp-1) ' + ' stnice{ik}(pp+1:end)];
   end;
   pm=findstr(stnice{ik},'-');
   pm(pm==1)=[];     % ignore minus on the real part 
   if ~isempty(pm),
      stnice{ik}=[stnice{ik}(1:pm-1) ' - ' stnice{ik}(pm+1:end)];
   end;   
end;   

