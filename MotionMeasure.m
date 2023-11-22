classdef MotionMeasure
    % Measure used to quantify motion.
    enumeration
        DVAR
        FD
    end
    methods (Static)
        function obj = getDefaultValue()
            % Default method for non-parametric combining is Stouffer's Z score.
            obj = MotionMeasure.FD;
        end
    end
end