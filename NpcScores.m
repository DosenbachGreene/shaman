classdef NpcScores < handle
    properties
        randomization_method RandomizationMethod {mustBeScalar,mustBeNonempty,mustRandomize} = RandomizationMethod.getDefaultValue() % How were permutations randomized
        npc_method NpcMethod {mustBeScalar,mustBeNonempty} = NpcMethod.getDefaultValue() % NpcMethod, Default: Stouffer
        score_type ScoreType {mustBeScalar,mustBeNonempty} = ScoreType.getDefaultValue() % ScoreType, Default: TwoSided
        t_thresh {mustBeNumeric,mustBeNonnegative,mustBeScalar} = 0 % t-value threshold for FalsePositive or FalseNegative scores, Default: NaN
        scores {mustBeNumeric} = [] % Matrix of non-parametric combining scores. By convention the dimensions are permutations (or bootstraps or 1 row if neither) x edges or nodes (or 1 column if combining across an entire connectivity matrix) x predictor variables x subsamples. Default: []
        p_values {mustBeNumeric,mustBeNonnegative} = [] % Matrix of p-values associated with scores. By convention the dimensions are permutations (or bootstraps or 1 row if neither) x edges or nodes (or 1 column if combining across an entire connectivity matrix) x predictor variables x subsamples. Default: []
        x_names string {mustBeVectorOrEmpty} = [] % String vector with names of predictor variables in the third dimension of NpcScores.scores. Default: []
        subsamples uint32 {mustBeUnique,mustBeInteger,mustBeVector,mustBeNonnegative,mustBeNonempty} = [0] % vector of subsample sizes; a size of 0 means no subsampling
    end
    methods
        function this = NpcScores(Args)
            arguments
                Args.npc_method NpcMethod {mustBeScalar,mustBeNonempty} = NpcMethod.getDefaultValue()
                Args.score_type ScoreType {mustBeScalar,mustBeNonempty} = ScoreType.getDefaultValue()
                Args.t_thresh {mustBeNumeric,mustBeNonnegative,mustBeScalar} = 0
                Args.scores {mustBeNumeric} = []
                Args.p_values {mustBeNumeric,mustBeNonnegative} = []
                Args.x_names string {mustBeVectorOrEmpty} = []
                Args.randomization_method RandomizationMethod {mustBeScalar,mustBeNonempty,mustRandomize} = RandomizationMethod.None
                Args.subsamples {mustBeUnique,mustBeInteger,mustBeVector,mustBeNonnegative,mustBeNonempty} = [0]
            end
            % Construct a new NpcScores object. If called with no arguments an
            % empty NpcScores object is created. Optional arguments are
            % npc_method, score_type, t_thresh, scores, p_values, and x_names,
            % which have the same meaning as they do as properties of NpcScores.

            this.npc_method = Args.npc_method;
            this.score_type = Args.score_type;
            this.t_thresh = Args.t_thresh;
            this.scores = Args.scores;
            this.p_values = Args.p_values;
            this.x_names = Args.x_names;
            if ~isempty(Args.randomization_method)
                this.randomization_method = Args.randomization_method;
            end
        end
        function compute_p_values(this, other)
            arguments
                this NpcScores
                other NpcScores
            end
            % Compute p-values for these scores using the permutations in other.
            %
            % other: An NpcScores object with more than one row of scores.
            this.p_values = zeros(size(this.scores));
            for i=1:size(this.scores,1)
                this.p_values(i,:,:,:) = mean(this.scores(i,:,:,:) < other.scores);
            end
        end
        function tbl = to_table(this)
            % Render motion impact scores and p-values as a table.

            % Just render the full sample as a table, no bootstrapping or
            % permutation.
            assert(this.subsamples(1) == 0);
            assert(size(this.scores,2) == 1);
            assert(size(this.scores,3) == length(this.x_names));
            
            tbl = table;
            tbl = addprop(tbl, ["score_type", "npc_method", "t_thresh", "randomization_method"], ["table", "table", "table", "table"]);
            scores = squeeze(this.scores(1,:,:,1));
            tbl.(this.score_type.to_string() + " Score") = scores(:);
            tbl.Properties.RowNames = this.x_names;
            tbl.Properties.CustomProperties.score_type = this.score_type;
            tbl.Properties.CustomProperties.npc_method = this.npc_method;
            tbl.Properties.CustomProperties.t_thresh = this.t_thresh;
            tbl.Properties.CustomProperties.randomization_method = this.randomization_method;
            tbl.Properties.Description = ...
                this.score_type.to_string() + ...
                " Motion Impact Score, " + ...
                this.npc_method.to_string() + ...
                " Non-parametric Combining Method, t-thresdhold = " + ...
                sprintf('%0.2f', this.t_thresh);

            if ~isempty(this.p_values)
                p_values = squeeze(this.p_values(1,:,:,1));
                tbl.(this.score_type.to_string() + " p Value") = p_values(:);
            end
        end
    end
end