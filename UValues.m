classdef UValues
    % Matrix of u-values and associated metadata.
    %
    % U-values are used in non-parametric combining.  See:
    % Winkler AM, Webster MA, Brooks JC, Tracey I, Smith SM, Nichols TE.
    % Non-parametric combination and related permutation tests for neuroimaging.
    % Hum Brain Mapp. 2016 Apr;37(4):1486-511. doi: 10.1002/hbm.23115
    properties
        score_type ScoreType {mustBeScalar,mustBeNonempty} = ScoreType.getDefaultValue() % ScoreType, Default: TwoSided
        u {mustBeNumeric} = [] % Matrix of u-values. By convention the dimensions are permutations (or 1 row if not permuted) x edges (columns) x predictor variables. Default: []
        t_thresh {mustBeScalar,mustBeNumeric,mustBeNonnegative} = 0; % Zero or a positive number used as a t-value threshold for FalsePositive and FalseNegative scores. Default: 0
        x_names string {mustBeVectorOrEmpty} = [] % String vector with names of predictor variables in the third dimension of UValues.u. Default: []
    end
    methods
        function this = UValues(Args)
            arguments
                Args.score_type ScoreType {mustBeScalar,mustBeNonempty} = ScoreType.getDefaultValue()
                Args.u {mustBeNumeric} = []
                Args.t_thresh {mustBeScalar,mustBeNumeric,mustBeNonnegative} = 0
                Args.x_names string {mustBeVectorOrEmpty} = []
            end
            % Construct a new UValues object. If called with no arguments an
            % empty UValues object is created. Optional arguments are score_type
            % u, t_thresh, and x_names, which have the same meaning as they
            % do as properties of Uvalues.

            this.score_type = Args.score_type;
            this.u = Args.u;
            this.t_thresh = Args.t_thresh;
            this.x_names = Args.x_names;
        end
    end
end