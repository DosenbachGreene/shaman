classdef ScoreType
    enumeration
        FalsePositive, FalseNegative, TwoSided
    end
    methods (Static)
        function obj = getDefaultValue()
           obj = ScoreType.TwoSided;
        end
    end
end