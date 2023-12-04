function mustBeUnique(x)
    assert(isvector(x));
    assert(length(x) == length(unique(x)));
end