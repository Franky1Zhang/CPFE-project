function value = tagr_bucket(tagr,tagr30,tagr70)
%BM_BUCKET Summary of this function goes here
%   row is a table class 

if isnan(tagr)
    value=blanks(1); 
elseif tagr<=tagr30
    value='A';
elseif tagr<=tagr70
    value='M';
elseif tagr>tagr70
    value='C';
else
    value=blanks(1);
    
end

