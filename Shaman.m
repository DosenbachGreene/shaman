classdef Shaman < handle
    properties (SetAccess=protected, GetAccess=public)
        data_provider % where data is sourced from
        x_names (1,:) string % names of variables used as main predictors
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
        function this = Shaman(data_provider, x_names, OptionalArgs)
            arguments
                data_provider DataProvider
                x_names (1,:) string {mustBeNonempty} % names of variables to use as main predictors
                OptionalArgs.nperm uint16 = 0
                OptionalArgs.intercept logical = true % include an intercept term
                OptionalArgs.motion_covariate logical = true % include motion as a covariate
                OptionalArgs.covariates (1,:) string = [] % cell array of variable names to include as covarriates
                OptionalArgs.show_progress logical = true % display progress
            end

            % Store arguments.
            this.data_provider = data_provider;
            this.x_names = x_names;
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
            this.full_model_fit = FullModelFit(model, this.x_names, "intercept", this.intercept, "motion_covariate", this.motion_covariate, "covariates", this.covariates, "show_progress", this.show_progress);

            % Fit the observed (not permuted) split model.
            if this.show_progress
                fprintf('Loading data for split model: ');
            end
            model = SplitModel(this.data_provider, "show_progress", this.show_progress);
            if this.show_progress
                fprintf('Fitting variables to split model: ');
            end
            this.split_model_fit = SplitModelFit(model, this.x_names, "intercept", this.intercept, "motion_covariate", this.motion_covariate, "covariates", this.covariates, "show_progress", this.show_progress);
            clear model;

            % Initialize permutation test.
            this.permutations = Permutations(this.data_provider, this.x_names, 'covariates', this.covariates, 'intercept', this.intercept, 'motion_covariate', this.motion_covariate, 'show_progress', this.show_progress);

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
                OptionalArgs.t_thresh {mustBeNumeric,mustBeScalar,mustBeNonnegative} = 2
                OptionalArgs.nodes {mustBeVectorOrEmpty,mustBeInteger,mustBePositive} = []
                OptionalArgs.edges {mustBeVectorOrEmpty,mustBeInteger,mustBePositive} = []
                OptionalArgs.show_progress logical = this.show_progress % display progress
            end

            % Validate x argument and convert to indices in this.x.
            xidx = this.xtoidx(OptionalArgs.x);

            % Validate nodes/edges arguments and convert to edges.
            edges = this.to_edges(OptionalArgs.nodes, OptionalArgs.edges);

            % Preallocate memory for u-values.
            u0 = zeros(1, length(edges), length(xidx));
            if nargout == 2
                u = zeros(this.permutations.nperm, length(edges), length(xidx));
            end

            % Display progress indicator.
            if OptionalArgs.show_progress
                line_length1 = fprintf("Computing %s u-values for variable ", OptionalArgs.score_type.to_string());
                line_length2 = fprintf("%d of %d", 1, length(xidx));
            end

            % Iterate over each variable in x and compute u-values.
            for i=1:length(xidx)
                if OptionalArgs.show_progress
                    fprintf(repmat('\b', 1, line_length2));
                    line_length2 = fprintf("%d of %d", i, length(xidx));
                end
                u0(1,:,i) = Shaman.compute_u_values(this.split_model_fit.t(xidx(i),edges), this.permutations.null_model_t(:,edges,xidx(i)), "full_model_t", this.full_model_fit.t(xidx(i),:), "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh).u;
                if nargout == 2
                    for j=1:this.permutations.nperm
                        u(j,:,i) = Shaman.compute_u_values(this.permutations.null_model_t(j,edges,xidx(i)), this.permutations.null_model_t(:,edges,xidx(i)), "full_model_t", this.full_model_fit.t(xidx(i),:), "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh).u;
                    end
                end
            end

            % Package result as a UValues object.
            u0 = UValues("u", u0, "score_type", OptionalArgs.score_type, "x_names", this.x_names(xidx), "t_thresh", OptionalArgs.t_thresh);
            if nargout == 2
                u = UValues("u", u, "score_type", OptionalArgs.score_type, "x_names", this.x_names(xidx), "t_thresh", OptionalArgs.t_thresh);
            end

            if OptionalArgs.show_progress
                fprintf(repmat('\b', 1, line_length1 + line_length2));
                fprintf("Computed %s u-values for %d variables.\n", u0.score_type.to_string(), length(xidx));
            end
        end
        function [npc0, npc] = get_npc_scores(this, OptionalArgs)
            arguments
                this Shaman
                OptionalArgs.x = []
                OptionalArgs.score_type ScoreType = ScoreType.getDefaultValue()
                OptionalArgs.npc_method NpcMethod = NpcMethod.getDefaultValue()
                OptionalArgs.t_thresh {mustBeNumeric,mustBeScalar,mustBeNonnegative} = 2
                OptionalArgs.nodes {mustBeVectorOrEmpty,mustBeInteger,mustBePositive} = []
                OptionalArgs.edges {mustBeVectorOrEmpty,mustBeInteger,mustBePositive} = []
                OptionalArgs.show_progress = this.show_progress;
            end

            % Validate x argument and convert to indices in this.x_names.
            xidx = this.xtoidx(OptionalArgs.x);

            % Validate nodes/edges arguments and convert to edges.
            edges = this.to_edges(OptionalArgs.nodes, OptionalArgs.edges);

            % Get u-values.
            if nargout == 1
                u0 = this.get_u_values("x", xidx, "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh, "edges", edges, "show_progress", OptionalArgs.show_progress);
            else
                [u0, u] = this.get_u_values("x", xidx, "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh, "edges", edges, "show_progress", OptionalArgs.show_progress);
            end

            % Convert u-values to npc scores.
            if nargout == 1
                npc0 = Shaman.compute_npc_scores(u0, "npc_method", OptionalArgs.npc_method, "show_progress", OptionalArgs.show_progress);
            else
                npc0 = Shaman.compute_npc_scores(u0, "npc_method", OptionalArgs.npc_method, "show_progress", false);
                npc = Shaman.compute_npc_scores(u, "npc_method", OptionalArgs.npc_method, "show_progress", OptionalArgs.show_progress);
            
                % Compute p-values.
                npc0.compute_p_values(npc);
            end
        end
        function tbl = get_scores_as_table(this, OptionalArgs)
            % Convenience method, takes same arguments as get_npc_scores()
            arguments
                this Shaman
                OptionalArgs.x = []
                OptionalArgs.score_type ScoreType = ScoreType.getDefaultValue()
                OptionalArgs.npc_method NpcMethod = NpcMethod.getDefaultValue()
                OptionalArgs.t_thresh {mustBeNumeric,mustBeScalar,mustBeNonnegative} = 2
                OptionalArgs.nodes {mustBeVectorOrEmpty,mustBeInteger,mustBePositive} = []
                OptionalArgs.edges {mustBeVectorOrEmpty,mustBeInteger,mustBePositive} = []
                OptionalArgs.show_progress logical = this.show_progress
                OptionalArgs.compute_p_values logical = true
            end
            compute_p_values = OptionalArgs.compute_p_values;
            OptionalArgs = rmfield(OptionalArgs, "compute_p_values");
            OptionalArgs = namedargs2cell(OptionalArgs);
            if compute_p_values
                [npc0, ~] = this.get_npc_scores(OptionalArgs{:});
            else
                npc0 = this.get_npc_scores(OptionalArgs{:});
            end
            tbl = npc0.to_table();
        end
        function npc0 = get_scores_by_node(this, OptionalArgs)
            arguments
                this Shaman
                OptionalArgs.x = []
                OptionalArgs.score_type ScoreType = ScoreType.getDefaultValue()
                OptionalArgs.npc_method NpcMethod = NpcMethod.getDefaultValue()
                OptionalArgs.t_thresh {mustBeNumeric,mustBeScalar,mustBeNonnegative} = 2
                OptionalArgs.show_progress logical = this.show_progress
                OptionalArgs.compute_p_values logical = true
            end

            xidx = this.xtoidx(OptionalArgs.x);

            if OptionalArgs.show_progress
                line_length1 = fprintf("Computing %s motion impact score\nusing %s non-parametric combining\non node ", OptionalArgs.score_type.to_string(), OptionalArgs.npc_method.to_string());
                line_length2 = fprintf("%d of %d", 1, this.n_nodes);
            end

            % Preallocate memory.
            npc0 = NpcScores("score_type", OptionalArgs.score_type, "npc_method", OptionalArgs.npc_method, "t_thresh", OptionalArgs.t_thresh, "x_names", this.x_names(xidx));
            npc0.scores = zeros(1, this.n_nodes, length(xidx));
            if OptionalArgs.compute_p_values
                npc0.p_values = zeros(1, this.n_nodes, length(xidx));
            end

            % Compute scores for each node.
            for i=1:this.n_nodes
                if OptionalArgs.show_progress
                    fprintf(repmat('\b', 1, line_length2));
                    line_length2 = fprintf("%d of %d", i, this.n_nodes);
                end
                if OptionalArgs.compute_p_values
                    [npc0i, ~] = this.get_npc_scores("nodes", [i], "x", xidx, "score_type", OptionalArgs.score_type, "npc_method", OptionalArgs.npc_method, "t_thresh", OptionalArgs.t_thresh, "show_progress", false);
                else
                    npc0i = this.get_npc_scores("nodes", [i], "x", xidx, "score_type", OptionalArgs.score_type, "npc_method", OptionalArgs.npc_method, "t_thresh", OptionalArgs.t_thresh, "show_progress", false);
                end
                npc0.scores(1,i,:) = npc0i.scores;
                if OptionalArgs.compute_p_values
                    npc0.p_values(1,i,:) = npc0i.p_values;
                end
            end

            if OptionalArgs.show_progress
                fprintf(repmat('\b', 1, line_length1 + line_length2));
                fprintf("Computed %s motion impact score\nusing %s non-parametric combining\non %d nodes.\n", OptionalArgs.score_type.to_string(), OptionalArgs.npc_method.to_string(), this.n_nodes);
            end
        end
    end
    methods (Static)
        function u = compute_u_values(t0, tperm, OptionalArgs)
            arguments
                t0 {mustBeNumeric,mustBeVector}
                tperm {mustBeNumeric,ismatrix}
                OptionalArgs.full_model_t {mustBeNumeric,mustBeVectorOrEmpty} = []
                OptionalArgs.score_type ScoreType = ScoreType.getDefaultValue()
                OptionalArgs.t_thresh {mustBeNumeric,mustBeScalar,mustBeNonnegative} = 2
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
            if OptionalArgs.score_type ~= ScoreType.TwoSided
                ft = OptionalArgs.full_model_t;
                thresh = OptionalArgs.t_thresh;
                u(ft < -thresh) = sum(tperm(:, ft < -thresh) < t0(ft < -thresh)) ./ size(tperm,1);
                u(ft > thresh) = sum(tperm(:, ft > thresh) > t0(ft > thresh)) ./ size(tperm,1);
            end

            % Package result as a UValues object.
            u = UValues("u", u, "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh);
        end
        function npc = compute_npc_scores(u, OptionalArgs)
            arguments
                u UValues
                OptionalArgs.npc_method NpcMethod = NpcMethod.getDefaultValue()
                OptionalArgs.show_progress logical = true
            end
            
            if OptionalArgs.show_progress
                line_length1 = fprintf("Performing %s non-parametric combining on variable ", OptionalArgs.npc_method.to_string());
                line_length2 = fprintf("%d of %d", 1, size(u.u,3));
            end

            % Perform non-parametric combining of u-values on the specified
            % edges using the specified method.
            npc = zeros(size(u.u,1), 1, size(u.u,3));
            for i=1:size(u.u,3)
                if OptionalArgs.show_progress
                    fprintf(repmat('\b', 1, line_length2));
                    line_length2 = fprintf("%d of %d", i, size(u.u,3));
                end
                for j=1:size(u.u,1)
                    npc(j,1,i) = NpcMethod.npc(u.u(j,:,i), OptionalArgs.npc_method);
                end
            end
            
            % Package into NpcScores objects.
            npc = NpcScores("scores", npc, "score_type", u.score_type, "t_thresh", u.t_thresh, "npc_method", OptionalArgs.npc_method, "x_names", u.x_names);

            if OptionalArgs.show_progress
                fprintf(repmat('\b', 1, line_length1 + line_length2));
                fprintf("Performed %s non-parametric combining on %d variables.\n", OptionalArgs.npc_method.to_string(), size(u.u,3));
            end
        end
    end
    methods (Access = private)
        function xidx = xtoidx(this, x)
            if iscell(x)
                % Find indicdes in this.x_names that match variable names in x.
                xidx = find(cellfun(@(a) any(cellfun(@(b) a == b || strcmp(a, b), x)), this.x_names));
                assert(length(xidx) == length(x), "Couldn't find all of x in Shaman.x");
            elseif isvector(x)
                if isstring(x)
                    % Find indicdes in this.x_names that match variable names in x.
                    xidx = arrayfun(@(x) find(this.x_names == x), x);
                    assert(length(xidx) == length(x), "Couldn't find all of x in Shaman.x");
                else
                    % Make sure x are valid indices in this.x_names.
                    mustBePositive(x);
                    mustBeInteger(x);
                    assert(max(x) <= length(this.x_names));
                    xidx = x;
                end
            elseif isempty(x)
                % If no x is provided, default to doing every variable in
                % this.x_names
                xidx = 1:length(this.x_names);
            else
                error("x does not index the values of Shaman.x");
            end
        end
        function edges = to_edges(this, nodes, edges)
            arguments
                this Shaman
                nodes {mustBeInteger,mustBePositive,mustBeVectorOrEmpty}
                edges {mustBeInteger,mustBePositive,mustBeVectorOrEmpty}
            end
            if isempty(nodes) && isempty(edges)
                edges = 1:this.n_edges;
            elseif isempty(nodes) && ~isempty(edges)
                assert(all(edges <= this.n_edges));
            elseif ~isempty(nodes) && isempty(edges)
                edges = nodes_to_edges("total_nodes", this.n_nodes, "nodes", nodes);
            else
                error("Cannot specify a subset of nodes and a subset of edges at the same time.");
            end
        end
    end
end