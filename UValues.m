classdef UValues
    properties
        score_type ScoreType {mustBeScalar,mustBeNonempty} = ScoreType.getDefaultValue()
        u {mustBeNumeric} = []
        t_thresh {mustBeScalar,mustBeNumeric,mustBeNonnegative} = 0;
    end
    methods
        function this = UValues(Args)
            arguments
                Args.score_type ScoreType {mustBeScalar,mustBeNonempty} = ScoreType.getDefaultValue()
                Args.u {mustBeNumeric} = []
                Args.t_thresh {mustBeScalar,mustBeNumeric,mustBeNonnegative} = 0;
            end
            this.score_type = Args.score_type;
            this.u = Args.u;
            this.t_thresh = Args.t_thresh;
        end
    end
end