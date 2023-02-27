classdef ScoreType
    % Type of motion impact score, e.g. false positive.
    enumeration
        FalsePositive, FalseNegative, TwoSided
    end
    methods
        function str = to_string(this)
            switch(this)
                case "FalsePositive"
                    str = "False Positive";
                case "FalseNegative"
                    str = "False Negative";
                case "TwoSided"
                    str = "Two Sided";
                otherwise
                    error("Unknown ScoreType")
            end
        end
    end
    methods (Static)
        function obj = getDefaultValue()
           obj = ScoreType.TwoSided;
        end
    end
end