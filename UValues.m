classdef UValues
    properties
        score_type ScoreType {mustBeScalar,mustBeNonempty} = ScoreType.getDefaultValue()
        u {mustBeNumeric} = []
        t_thresh {mustBeScalar,mustBeNumeric,mustBeNonnegative} = 0;
        x_names (1,:) string = []
    end
    methods
        function this = UValues(Args)
            arguments
                Args.score_type ScoreType {mustBeScalar,mustBeNonempty} = ScoreType.getDefaultValue()
                Args.u {mustBeNumeric} = []
                Args.t_thresh {mustBeScalar,mustBeNumeric,mustBeNonnegative} = 0;
                Args.x_names = []
            end
            this.score_type = Args.score_type;
            this.u = Args.u;
            this.t_thresh = Args.t_thresh;
            this.x_names = Args.x_names;
        end
    end
end