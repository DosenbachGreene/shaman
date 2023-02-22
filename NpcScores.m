classdef NpcScores < handle
    properties
        npc_method NpcMethod {mustBeScalar,mustBeNonempty} = NpcMethod.getDefaultValue()
        score_type ScoreType {mustBeScalar,mustBeNonempty} = ScoreType.getDefaultValue()
        t_thresh {mustBeNumeric,mustBeNonnegative,mustBeScalar} = 0
        scores {mustBeNumeric} = []
        p_values {mustBeNumeric,mustBeNonnegative} = []
    end
    methods
        function this = NpcScores(Args)
            arguments
                Args.npc_method NpcMethod {mustBeScalar,mustBeNonempty} = NpcMethod.getDefaultValue()
                Args.score_type ScoreType {mustBeScalar,mustBeNonempty} = ScoreType.getDefaultValue()
                Args.t_thresh {mustBeNumeric,mustBeNonnegative,mustBeScalar} = 0
                Args.scores {mustBeNumeric} = []
                Args.p_values {mustBeNumeric,mustBeNonnegative} = []
            end
            this.npc_method = Args.npc_method;
            this.score_type = Args.score_type;
            this.t_thresh = Args.t_thresh;
            this.scores = Args.scores;
            this.p_values = Args.p_values;
        end
    end
end