function tf = anynan(x)
    arguments 
        x {mustBeNumeric}
    end
    % Returns true if any element in matrix x is NaN, false otherwise.
    
    tf = any(isnan(x(:)));
end