function value = roe_bucket(roe,roe30,roe70)
%BM_BUCKET Summary of this function goes here
%   row is a table class 

if isnan(roe)
    value=blanks(1); 
elseif roe<=roe30
    value='W';
elseif roe<=roe70
    value='M';
elseif roe>roe70
    value='R';
else
    value=blanks(1);
    
end

