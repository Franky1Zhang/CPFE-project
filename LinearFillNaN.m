function result=LinearFillNaN(input)
%Do Linear interpolation when meet NaN.Input and return are arrays.
for i=1:length(input)
    if isnan(input(i))
        n=0;
        value_left=input(i-1);
        while isnan(input(i+n))
            n=n+1;
        end
        value_right=input(i+n);
        delta=(value_right-value_left)/(n+1);
        for j=1:n
            input(i+j-1)=value_left+j*delta;
        end
    end
end

result=input;
end