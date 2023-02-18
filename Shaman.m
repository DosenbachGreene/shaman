classdef Shaman < handle
    properties (SetAccess=protected, GetAccess=public)
        data_provider % where data is sourced from
        x % names of variables used as main predictors
        split_model_fit % struct array of observed regression coefficients and t-values for the observed (not permuted) data, one element for each variable in x
        full_model_fit
        permutations % n permutations (rows) x edges (columns) x variables in x (3rd dimension) array of t-values
        intercept % whether the model fit includes an interecept column
        covariates % cell array of variables included as covariates
        motion_covariate % whether motion was included as a covariate
        n_nodes
        n_edges
    end
    properties
        show_progress % whether to show a progress indicator
    end
    methods
        function this = Shaman(data_provider, x, OptionalArgs)
            arguments
                data_provider DataProvider
                x (1,:) cell {mustBeText,mustBeNonempty} % names of variables to use as main predictors
                OptionalArgs.nperm uint16 = 0
                OptionalArgs.intercept logical = true % include an intercept term
                OptionalArgs.motion_covariate logical = true % include motion as a covariate
                OptionalArgs.covariates (1,:) cell = {} % cell array of variable names to include as covarriates
                OptionalArgs.show_progress logical = true % display progress
            end

            % Store arguments.
            this.data_provider = data_provider;
            this.x = x;
            this.intercept = OptionalArgs.intercept;
            this.motion_covariate = OptionalArgs.motion_covariate;
            this.covariates = OptionalArgs.covariates;
            this.show_progress = OptionalArgs.show_progress;

            % Fit the full model.
            if this.show_progress
                fprintf('Loading data for full model: ');
            end
            model = FullModel(this.data_provider, "Progress", this.show_progress);
            if this.show_progress
                fprintf('Fitting variables to full model: ');
            end
            this.full_model_fit = FullModelFit(model, this.x, "intercept", this.intercept, "motion_covariate", this.motion_covariate, "covariates", this.covariates, "show_progress", this.show_progress);

            % Fit the observed (not permuted) split model.
            if this.show_progress
                fprintf('Loading data for split model: ');
            end
            model = SplitModel(this.data_provider, "show_progress", this.show_progress);
            if this.show_progress
                fprintf('Fitting variables to split model: ');
            end
            this.split_model_fit = SplitModelFit(model, this.x, "intercept", this.intercept, "motion_covariate", this.motion_covariate, "covariates", this.covariates, "show_progress", this.show_progress);
            clear model;

            % Initialize permutation test.
            this.permutations = Permutations(this.data_provider, this.x, 'covariates', this.covariates, 'intercept', this.intercept, 'motion_covariate', this.motion_covariate, 'show_progress', this.show_progress);

            % Perform permutations.
            if OptionalArgs.nperm > 0
                this.permutations.nperm = OptionalArgs.nperm;
            end

            % Store number of nodes and edges.
            this.n_edges = size(this.full_model_fit.t,2);
            this.n_nodes = (1 + sqrt(1 + 8*this.n_edges)) / 2; % infer fro # of edges
        end
        function [u0, u] = get_u_values(this, OptionalArgs)
            arguments
                this Shaman
                OptionalArgs.x = []
                OptionalArgs.score_type ScoreType = ScoreType.getDefaultValue()
                OptionalArgs.t_thresh {mustBeNumeric,mustBeScalar} = 2
            end

            % Validate x argument and convert to indices in this.x.
            xidx = this.xtoidx(OptionalArgs.x);

            % Preallocate memory for u-values.
            u0 = zeros(length(xidx), size(this.split_model_fit.t,2));
            if nargout == 2
                u = zeros(this.permutations.nperm, size(this.split_model_fit.t,2), length(xidx));
            end

            % Iterate over each variable in x and compute u-values.
            for i=1:length(xidx)
                u0(i,:) = Shaman.get_u_values_(this.split_model_fit.t(xidx(i),:), this.permutations.null_model_t(:,:,xidx(i)), "full_model_t", this.full_model_fit.t(xidx(i),:), "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh).u;
                if nargout == 2
                    for j=1:this.permutations.nperm
                        u(j,:,i) = Shaman.get_u_values_(this.permutations.null_model_t(j,:,xidx(i)), this.permutations.null_model_t(:,:,xidx(i)), "full_model_t", this.full_model_fit.t(xidx(i),:), "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh).u;
                    end
                end
            end

            % Package result as a UValues object.
            u0 = UValues("u", u0, "score_type", OptionalArgs.score_type);
            if nargout == 2
                u = UValues("u", u, "score_type", OptionalArgs.score_type);
            end
        end
        function [npc0, npc] = get_npc_scores(this, OptionalArgs)
            arguments
                this Shaman
                OptionalArgs.x = []
                OptionalArgs.score_type ScoreType = ScoreType.getDefaultValue()
                OptionalArgs.npc_method NpcMethod = NpcMethod.getDefaultValue()
                OptionalArgs.t_thresh {mustBeNumeric,mustBeScalar} = 2
                OptionalArgs.nodes {mustBeVector,mustBeInteger,mustBePositive} = []
                OptionalArgs.edges {mustBeVector,mustBeInteger,mustBePositive} = []
            end

            % Validate x argument and convert to indices in this.x.
            xidx = this.xtoidx(OptionalArgs.x);

            % Validate nodes/edges arguments and convert to edges.
            edges = this.to_edges(OptionalArgs.nodes, OptionalArgs.edges);

            % Preallocate memory for u-values.
            npc0 = zeros(length(xidx), 1);
            if nargout == 2
                npc = zeros(this.permutations.nperm, size(this.split_model_fit.t,2), length(xidx));
            end

            % Iterate over each variable in x and compute u-values.
            for i=1:length(xidx)
                u0(i,:) = Shaman.get_u_values_(this.split_model_fit.t(xidx(i),:), this.permutations.null_model_t(:,:,xidx(i)), "full_model_t", this.full_model_fit.t(xidx(i),:), "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh).u;
                if nargout == 2
                    for j=1:this.permutations.nperm
                        u(j,:,i) = Shaman.get_u_values_(this.permutations.null_model_t(j,:,xidx(i)), this.permutations.null_model_t(:,:,xidx(i)), "full_model_t", this.full_model_fit.t(xidx(i),:), "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh).u;
                    end
                end
            end

            % Compute u-values.
            if nargout == 1
                u0 = this.get_u_values("score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh);
            else
                [u0, u] = this.get_u_values("score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh);
            end

            % Perform non-parametric combining on a subset of nodes.
            npc0 = NpcMethod.npc(u0(edges), npc_method);
            p_value = [];
            if nargout > 1
                npc = zeros(1,this.permutations.nperm);
                for i=1:this.permutations.nperm
                    npc(i) = NpcMethod.npc(u(i,edges), npc_method);
                end
                p_value = sum(npc > npc0) / this.permutations.nperm;
            end
            
            % Package into NpcScores objects.
            npc0 = NpcScores("scores", npc0, "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh, "npc_method", OptionalArgs.npc_method, "p_values", p_value);
            npc = NpcScores("scores", npc, "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh, "npc_method", OptionalArgs.npc_method);
        end
        function get_motion_impact_score(this, OptionalArgs)
            arguments
                this Shaman
                OptionalArgs.x = []
                OptionalArgs.score_type (1,:) cell {mustBeScoreType}
                OptionalArgs.t_thresh {mustBeNumeric,mustBeScalar} = 2
                OptionalArgs.nodes {mustBeVector,mustBeInteger,mustBePositive} = []
                OptionalArgs.edges {mustBeVector,mustBeInteger,mustBePositive} = []
            end

            % Validate x argument and convert to indices in this.x.
            xidx = this.xtoidx(OptionalArgs.x);

            % We can work on a subset of nodes or edges, but not both.
            assert(and(isempty(OptionalArgs.nodes), isempty(OptionalArgs.nodes)) || xor(~isempty(OptionalArgs.nodes), ~isempty(OptionalArgs.nodes)), "Cannot specify a subset of nodes and a subset of edges at the same time.");

            error('Not implemented.');
        end
        function get_motion_impact_score_by_node(this)
            error("Not Implemented.");
        end
    end
    methods (Static)
        function u = get_u_values_(t0, tperm, OptionalArgs)
            arguments
                t0 {mustBeNumeric,mustBeVector}
                tperm {mustBeNumeric,ismatrix}
                OptionalArgs.full_model_t {mustBeNumeric,mustBeVectorOrEmpty} = []
                OptionalArgs.score_type ScoreType = ScoreType.getDefaultValue()
                OptionalArgs.t_thresh {mustBeNumeric,mustBeScalar} = 2
            end

            % Validate arguments.
            if OptionalArgs.score_type ~= ScoreType.TwoSided
                assert(~isempty(OptionalArgs.full_model_t), 'Need t-values from full model to compute false positive or false negative motion impact score.');
            end

            % If we are computing false negative motion impact score,
            % simply flip the signs on the full model's t-values.
            if OptionalArgs.score_type == ScoreType.FalseNegative
                OptionalArgs.full_model_t = -OptionalArgs.full_model_t;
            end

            % Begin by computing two-sided u-values.
            u = sum(abs(tperm) > abs(t0)) ./ size(tperm,1);

            % Then compute one-sided u-values for edges above the t-value
            % threshold.
            ft = OptionalArgs.full_model_t;
            thresh = OptionalArgs.t_thresh;
            if OptionalArgs.score_type ~= ScoreType.TwoSided
                u(ft < -thresh) = sum(tperm(:, ft < -thresh) < t0(ft < -thresh)) ./ size(tperm,1);
                u(ft > thresh) = sum(tperm(:, ft > thresh) > t0(ft > thresh)) ./ size(tperm,1);
            end

            % Package result as a UValues object.
            u = UValues("u", u, "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh);
        end
    end
    methods (Access = private)
        function xidx = xtoidx(this, x)
            if iscell(x)
                % Find indicdes in this.x that match variable names in x.
                xidx = find(cellfun(@(a) any(cellfun(@(b) a == b || strcmp(a, b), x)), this.x));
                assert(length(xidx) == length(x), "Couldn't find all of x in Shaman.x");
            elseif isvector(x)
                if isstring(foo)
                    % Find indicdes in this.x that match variable names in x.
                    xidx = arrayfun(@(x) find(this.x == x), x);
                    assert(length(xidx) == length(x), "Couldn't find all of x in Shaman.x");
                else
                    % Make sure x are valid indices in this.x.
                    mustBePositive(x);
                    mustBeInteger(x);
                    assert(max(x) <= length(this.x));
                    xidx = x;
                end
            elseif isempty(x)
                % If no x is provided, default to doing every variable in
                % this.x
                xidx = 1:length(this.x);
            else
                error("x does not index the values of Shaman.x");
            end
        end
        function edges = toedges(this, nodes, edges)
            arguments
                this SHaman
                nodes {mustBeInteger,mustBePositive,mustBeVector}
                edges {mustBeInteger,mustBePositive,mustBeVector}
            end
            if isempty(nodes) && isempty(edges)
                edges = 1:this.n_edges;
            elseif isempty(nodes) && ~isempty(edges)
                assert(all(edges < this.n_edges));
            elseif ~isempty(OptionalArgs.nodes) && isempty(OptionalArgs.edges)
                edges = nodes_to_edges("total_nodes", this.n_nodes, "nodes", nodes);
            else
                error("Cannot specify a subset of nodes and a subset of edges at the same time.");
            end
        end
    end
end