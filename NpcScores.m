classdef NpcScores < handle
    properties
        npc_method NpcMethod {mustBeScalar,mustBeNonempty} = NpcMethod.getDefaultValue()
        score_type ScoreType {mustBeScalar,mustBeNonempty} = ScoreType.getDefaultValue()
        t_thresh {mustBeNumeric,mustBeNonnegative,mustBeScalar} = 0
        scores {mustBeNumeric} = []
        p_values {mustBeNumeric,mustBeNonnegative} = []
        x_names (1,:) string = []
    end
    methods
        function this = NpcScores(Args)
            arguments
                Args.npc_method NpcMethod {mustBeScalar,mustBeNonempty} = NpcMethod.getDefaultValue()
                Args.score_type ScoreType {mustBeScalar,mustBeNonempty} = ScoreType.getDefaultValue()
                Args.t_thresh {mustBeNumeric,mustBeNonnegative,mustBeScalar} = 0
                Args.scores {mustBeNumeric} = []
                Args.p_values {mustBeNumeric,mustBeNonnegative} = []
                Args.x_names (1,:) string = []
            end
            this.npc_method = Args.npc_method;
            this.score_type = Args.score_type;
            this.t_thresh = Args.t_thresh;
            this.scores = Args.scores;
            this.p_values = Args.p_values;
            this.x_names = Args.x_names;
        end
        function compute_p_values(this, other)
            arguments
                this NpcScores
                other NpcScores
            end
            assert(size(this.scores,1) == 1);
            assert(size(this.scores,3) == size(other.scores,3));
            this.p_values = zeros(1,1,size(this.scores,3));
            for i=1:size(this.scores,3)
                this.p_values(i) = mean(this.scores(1,1,i) < other.scores(:,1,i));
            end
        end
        function tbl = to_table(this)
            % Render motion impact scores and p-values as a table.
            assert(size(this.scores, 1) == 1, "Cannot render motion impact scores from multiple permutations as table.");
            assert(size(this.scores,2) == 1);
            assert(size(this.scores,3) == length(this.x_names));
            
            tbl = table;
            tbl = addprop(tbl, ["score_type", "npc_method", "t_thresh"], ["table", "table", "table"]);
            tbl.(this.score_type.to_string() + " Score") = this.scores(:);
            tbl.Properties.RowNames = this.x_names;
            tbl.Properties.CustomProperties.score_type = this.score_type;
            tbl.Properties.CustomProperties.npc_method = this.npc_method;
            tbl.Properties.CustomProperties.t_thresh = this.t_thresh;
            tbl.Properties.Description = ...
                this.score_type.to_string() + ...
                " Motion Impact Score, " + ...
                this.npc_method.to_string() + ...
                " Non-parametric Combining Method, t-thresdhold = " + ...
                sprintf('%0.2f', this.t_thresh);

            if length(this.p_values) == length(this.x_names)
                tbl.(this.score_type.to_string() + " p Value") = this.p_values;
            end
        end
    end
end