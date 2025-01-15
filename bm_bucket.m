function value = bm_bucket(beme,bm30,bm70)
%BM_BUCKET Summary of this function goes here
%   row is a table class 
if isnan(beme)
    value=blanks(1); 
elseif beme<=bm30
    value='L';
elseif beme<=bm70
    value='M';
elseif beme>bm70
    value='H';
else
    value=blanks(1);
    
end

