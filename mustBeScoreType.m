function mustBeScoreType(x)
    if iscell(x)
        assert(all(cellfun(@(x) isa(x,"ScoreType"), x)));
    else
        assert(all(isa(x,"ScoreType")));
    end
end