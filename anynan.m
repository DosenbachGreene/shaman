function tf = anynan(x)
    arguments 
        x {isNumeric}
    end
    % Returns true if any element in matrix x is NaN, false otherwise.
    
    tf = any(isnan(data.fmri(:)));
end