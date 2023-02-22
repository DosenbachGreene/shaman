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
    end
end