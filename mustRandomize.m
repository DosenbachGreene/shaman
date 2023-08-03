function mustRandomize(x)
    arguments
        x RandomizationMethod
    end
    assert(x ~= RandomizationMethod.None);
end