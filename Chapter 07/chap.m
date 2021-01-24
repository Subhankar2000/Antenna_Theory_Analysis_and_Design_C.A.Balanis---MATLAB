function [st]=chap(SFAF,str_mod);
st=0;
if norm(SFAF-fliplr(SFAF))>=1e-3,
    disp(' ');
    disp(['Warning: Specified ' str_mod ' is not symmetrical with respect to 90 deg.' ...
          ' Phase information is missing and']); 
    disp(['cannnot be uniquely reconstructed. Program will continue but another synthesis method ' ...
         'is strongly suggested.']);
    st=1;    
end;   

if abs(SFAF(1))>=1e-3|abs(SFAF(end))>=1e-3,
   disp(' ');
   disp(['Warning: Specified ' str_mod ' does not have compact support. Invisible region cannot be ' ...
         'uniquely defined.']);
   disp('Program will continue but another synthesis method is strongly suggested.');
   st=1;
end;   
