classdef Shaman < handle
    % The Shaman class brings all the functionality you need to compute
    % motion impact score into one convenient place. Thus is the class you
    % will interact with most.
    %
    % To instantiate a Shaman object you need a DataProvider and the names
    % of the variables you want to compute motion impact score for.
    %
    %     data_provider = MyDataProvider("where_to_find_my_data");
    %     shaman = Shaman(data_provider, ["trait1", "trait2", "trait3"]);
    %
    % Once you have a shaman object, you will probably want to run some
    % permutations so that you can perform statistics on the motion impact
    % scores.
    %
    %     shaman.permutations.nperm = 1000; % do a thousand permutations
    %
    % Note that you can do more permutations without losing your progress.
    %
    %     shaman.permutations.nperm = 2000; % do another 1000 permutations
    %
    % To quicky get a table of motion impact scores and p-values:
    %
    %     tbl = shaman.get_scores_as_table("score_type", ScoreType.FalsePositive);

    properties (SetAccess=protected, GetAccess=public)
        data_provider DataProvider = NullDataProvider() % A DataProvider, the source of imaging and non-imaging (e.g. behavioral, biophysical, etc) data.
        x_names string {mustBeVector,mustBeNonempty} = [0] % Names of non-imaging variables provided by the DataProvider upon which to compute motion impact score. A separate score is computed independently for each variable. Example: ["trait1", "trait2"]
        full_model_fit FullModelFit % ModelFit object for the "full" (not split-half) connectivity matrix for each variable.
        split_model_fit SplitModelFit % ModelFit object for the split-half (by motion) connectivity matrix for each variable. This is the raw motion impact before permutation and non-parametric combining to get a motion impact score.
        subsamples {mustBeUnique,mustBeInteger,mustBeVector,mustBeNonnegative,mustBeNonempty} = [0] % vector of subsample sizes; a size of 0 means no subsampling
        permutations Permutations % Permutations object.  Set shaman.permutations.nperm = 1000 to run a thousand (or whatever number you want) permutations.
        intercept logical = true % Whether to model an intercept term. Default: true
        motion_covariate logical = true % Whether to model motion as a covariate. Default: true
        covariates string {mustBeVectorOrEmpty} % Names of additional non-imaging variables, besides motion and an intercept term, to include as covariates. The default is not to include any extra covariates. Example: ["cov1", "cov2"]
        n_nodes {mustBeInteger,mustBePositive} % Number of nodes (i.e. voxels, vertices, regions, parcels) in the model.
        n_edges {mustBeInteger,mustBePositive} % Number of edges (i.e. pairwise connections) in the model.
        tbl table % Table of non-imaging data.
        motion {mustBeVector,mustBeNumeric} = [0] % Motion ccovariate for full model.
    end
    properties
        max_par_workers uint32 = maxNumCompThreads() % maximum number of parallel workers to use. Default: maxNumCompThreads()
        show_progress logical = true % Whether to show a progress indicator for operations that could take a long time. Default: true
    end
    methods
        function this = Shaman(data_provider, x_names, OptionalArgs)
            arguments
                data_provider DataProvider
                x_names string {mustBeVector,mustBeNonempty}
                OptionalArgs.subsamples {mustBeUnique,mustBeInteger,mustBeVector,mustBeNonnegative,mustBeNonempty} = [0]
                OptionalArgs.nboot {mustBeNonnegative,mustBeInteger,mustBeScalar} = 1 
                OptionalArgs.nperm {mustBeInteger,mustBeNonnegative} = 0
                OptionalArgs.intercept logical = true
                OptionalArgs.motion_covariate logical = true
                OptionalArgs.covariates string {mustBeVectorOrEmpty} = []
                OptionalArgs.randomization_method RandomizationMethod {mustRandomize} = RandomizationMethod.getDefaultValue();
                OptionalArgs.max_par_workers {mustBeNumeric,mustBeNonnegative,mustBeScalar} = maxNumCompThreads()
                OptionalArgs.show_progress logical = true
            end
            % Construct a new Shaman object.
            %
            % You must supply a DataProvider and a list of non-imaging
            % variables upon which you would like to compute a motion
            % impact score.  For example:
            %
            %     shaman = Shaman(my_data_provider, ["var1", "var2"]);
            %
            % You can specify optional arguments using the syntax:
            %
            %     shaman = Shaman(data_provider, x_names, "name", value);
            %
            % Optional arguments:
            %
            % nperm:
            %     Number of permutations to perform. Default: 0
            %     You can also set the property
            %     Shaman.permutations.nperm to perform additional permutations
            %     after the Shaman object has been created.
            % intercept: Whether to model an intercept term. Default: true
            % motion_covariate:
            %     Whether to model motion as a covariate. Default: true
            % covariates:
            %     Names of additional non-imaging variables to model as
            %     covariates. Default: []
            %     Example: ["cov1", "cov2"]
            % randomization_method:
            %     Randomization method to use for permutation testing.
            %     See RandomizationMethod.
            % max_par_workers:
            %     Use at most this many parallel workers.
            %     default: maxNumCompThreads()
            % show_progress: Whether to show a progress indicator. Default: true

            % Store arguments in self.
            this.data_provider = data_provider;
            this.x_names = x_names;
            this.subsamples = unique(OptionalArgs.subsamples); % sort and make sure subsamples are unique
            this.intercept = OptionalArgs.intercept;
            this.motion_covariate = OptionalArgs.motion_covariate;
            this.covariates = OptionalArgs.covariates;
            this.max_par_workers = OptionalArgs.max_par_workers;
            this.show_progress = OptionalArgs.show_progress;
            nboot = OptionalArgs.nboot;

            % Load data for the full model.
            if this.show_progress
                fprintf('Loading data for full model: ');
            end
            model = FullModel(this.data_provider, "show_progress", this.show_progress);
            % Save non-imaging data from the model.
            this.tbl = model.tbl;
            this.motion = model.motion;
            % Fit the full model.
            if this.show_progress
                fprintf('Fitting variables to full model: ');
            end
            this.full_model_fit = FullModelFit(model, this.x_names, "intercept", this.intercept, "motion_covariate", this.motion_covariate, "covariates", this.covariates, "show_progress", this.show_progress, "subsamples", this.subsamples, "nboot", nboot);

            % Fit the observed (not permuted) split model.
            if this.show_progress
                fprintf('Loading data for split model: ');
            end
            model = SplitModel(this.data_provider, "show_progress", this.show_progress);
            if this.show_progress
                fprintf('Fitting variables to split model: ');
            end
            this.split_model_fit = SplitModelFit(model, this.x_names, "intercept", this.intercept, "motion_covariate", this.motion_covariate, "covariates", this.covariates, "subsamples", this.subsamples, "nboot", nboot, "show_progress", this.show_progress);
            clear model;

            % Initialize permutation test.
            this.permutations = Permutations(this.data_provider, this.x_names, "covariates", this.covariates, "intercept", this.intercept, "motion_covariate", this.motion_covariate, "randomization_method", OptionalArgs.randomization_method, "max_par_workers", this.max_par_workers, "subsamples", this.subsamples, "show_progress", this.show_progress);

            % Perform permutations.
            if OptionalArgs.nperm > 0
                this.permutations.nperm = OptionalArgs.nperm;
            end

            % Store number of nodes and edges.
            this.n_edges = size(this.full_model_fit.t,2);
            this.n_nodes = (1 + sqrt(1 + 8*this.n_edges)) / 2;
        end

        function set.max_par_workers(this, val)
            arguments
                this Shaman
                val {mustBeNumeric,mustBeNonnegative,mustBeScalar}
            end
            this.max_par_workers = val;
            if ~isempty(this.permutations)
            	this.permutations.max_par_workers = val;
            end
        end
        function [u0, u] = get_u_values(this, OptionalArgs)
            arguments
                this Shaman
                OptionalArgs.x = []
                OptionalArgs.score_type ScoreType = ScoreType.getDefaultValue()
                OptionalArgs.t_thresh {mustBeNumeric,mustBeScalar,mustBeNonnegative} = 2
                OptionalArgs.nodes {mustBeVectorOrEmpty,mustBeInteger,mustBePositive} = []
                OptionalArgs.edges {mustBeVectorOrEmpty,mustBeInteger,mustBePositive} = []
                OptionalArgs.show_progress logical = this.show_progress
            end
            % Compute u-values, the first step in non-parametric combining.
            %
            % Compares the t-value at each edge in the split half model to
            % the permuted t-values for each corresponding edge in the null
            % model to obtain a u-value. The u-value is, conceptually, an
            % (not corrected for multiple comparisons) p-value. Anderson
            % Winkler recommends the terminology "u-value" to draw a
            % distinction from the omnibus p-value obtained at the end of
            % non-parametric combining.
            %
            % u0 = Shaman.get_u_values() returns just the u-values for the
            % not-permuted split half model in a 1 x edges x variables
            % matrix wrapped in a UValues object.
            %
            % [u0, u] = Shaman.get_u_values() returns the u-values for the
            % split half model as u0 and the u-values for each permutation
            % in an nperm x edges x variables matrix u, each wrapped in a
            % UValues object.
            %
            % Optional arguments:
            %
            % x:
            %     Vector of names of variables, e.g. ["var1", "var2"], or
            %     indices of variables in Shaman.x_names, e.g. [1, 2], to
            %     compute u-values for.  Defaults to [], for which all the
            %     variables in Shaman.x_names are used.
            % score_type: A ScoreType. Default: ScoreType.TwoSided
            % t_thresh:
            %     A t-value threshold. The threshold is used for FalsePositive
            %     or FalseNegative tests. Default: 2
            % nodes:
            %     Vector specifying a subset of nodes to analyze, e.g.
            %     [1,2,3]. Defaults to [], for which all nodes are analyzed.
            %     Cannot be mixed with edges.
            % edges:
            %     Vector specifying a subset of edges to analyze, e.g.
            %     [1, 2, 3]. Default to []], for which all edges are analyzed.
            %     Cannot be mixed with nodes.
            % show_progress:
            %     Whether to show a progress indicator. Default to the
            %     Shaman.show_progress property.

            % Validate x argument and convert to indices in this.x_names.
            xidx = this.xtoidx(OptionalArgs.x);

            % Validate nodes/edges arguments and convert to edges.
            edges = this.to_edges(OptionalArgs.nodes, OptionalArgs.edges);

            % TODO check if xidx and edges arguments are respected.
            u0 = Shaman.compute_u_values(this.split_model_fit.t, this.permutations.null_model_t, "randomization_method", this.permutations.randomization_method, "full_model_t", this.full_model_fit.t, "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh, "x_names", this.x_names(xidx), "subsamples", this.subsamples);
            if nargout > 1
                u = Shaman.compute_u_values(this.permutations.null_model_t, this.permutations.null_model_t, "randomization_method", this.permutations.randomization_method, "full_model_t", this.full_model_fit.t, "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh, "x_names", this.x_names(xidx), "subsamples", this.subsamples);
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
            % Compute scores, the final step in non-parametric combining.
            %
            % This function automatically computes u-values, the first step
            % in non-parametric combining. See Shaman.get_u_values() for
            % details. Then this function perfoms non-parametric combining
            % (npc) across edges to compute an omnibus motion impact score.
            %
            % npc0 = Shaman.get_npc_scores() returns just the motion impact
            % score for not-permuted split half model in a
            % 1 x 1 x variables matrix wrapped in an NpcScores object.
            %
            % [npc0, npc] = Shaman.get_npc_scores() returns the motion
            % impact score for the split half model as npc0 and the scores
            % for each permutation in an nperm x 1 x variables matrix npc,
            % each wrapped in an NpcScores object. It also automatically
            % computes the p-values in npc0.
            %
            % Optional arguments:
            %
            % x:
            %     Vector of names of variables, e.g. ["var1", "var2"], or
            %     indices of variables in Shaman.x_names, e.g. [1, 2], to
            %     compute u-values for. Defaults to [], for which all the
            %     variables in Shaman.x_names are used.
            % score_type: A ScoreType. Default: ScoreType.TwoSided
            % npc_method: An NpcMethod, Default: Stouffer
            % t_thresh:
            %     A t-value threshold. The threshold is used for FalsePositive
            %     or FalseNegative tests. Default: 2
            % nodes:
            %     Vector specifying a subset of nodes to analyze, e.g.
            %     [1, 2, 3]. Defaults to [], for which all nodes are analyzed.
            %     Cannot be mixed with edges.
            % edges:
            %     Vector specifying a subset of edges to analyze, e.g.
            %     [1, 2, 3]. Defaults to [], for which all edges are analyzed.
            %     Cannot be mixed with nodes.
            % show_progress:
            %     Whether to show a progress indicator. Default to the
            %     Shaman.show_progress property.

            % Validate x argument and convert to indices in this.x_names.
            xidx = this.xtoidx(OptionalArgs.x);

            % Validate nodes/edges arguments and convert to edges.
            edges = this.to_edges(OptionalArgs.nodes, OptionalArgs.edges);

            % TODO check if xidx and edges arguments are respected.
            [u0, u] = this.get_u_values("x", xidx, "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh, "nodes", OptionalArgs.nodes, "edges", OptionalArgs.edges, "show_progress", OptionalArgs.show_progress); % TODO args
            npc0 = Shaman.compute_npc_scores(u0, "npc_method", OptionalArgs.npc_method, "show_progress", OptionalArgs.show_progress);
            if nargout > 1
                npc = Shaman.compute_npc_scores(u, "npc_method", OptionalArgs.npc_method, "show_progress", OptionalArgs.show_progress);
                npc0.compute_p_values(npc);
            end
            
            
        end
        function tbl = get_scores_as_table(this, OptionalArgs)
            % Convenience method to get motion impact scores in table form.
            %
            % Takes the same optional arguments as Shaman.get_npc_scores().
            % This function Generates a human-readable table of scores
            % intead of an NpcScores object.
            %
            % Call as Shaman.get_scores_as_table("compute_p_values", false)
            % to suppress computation of p-values.
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
            % Convenience method to compute a score for each node.
            %
            % Takes the same optional arguments as Shaman.get_npc_scores().
            % Returns a 1 x nodes x variables matrix of motion impact
            % scores in npc0, wrapped inside an NpcScores object. This is
            % useful for generating figures showing motion impact score
            % projected onto the brain.
            %
            % Call as Shaman.get_scores_by_node("compute_p_values", false)
            % to suppress computation of p-values.
            arguments
                this Shaman
                OptionalArgs.x = []
                OptionalArgs.score_type ScoreType = ScoreType.getDefaultValue()
                OptionalArgs.npc_method NpcMethod = NpcMethod.getDefaultValue()
                OptionalArgs.t_thresh {mustBeNumeric,mustBeScalar,mustBeNonnegative} = 2
                OptionalArgs.show_progress logical = this.show_progress
                OptionalArgs.compute_p_values logical = true
            end

            % Validate x argument and convert to indices in this.x_names.
            xidx = this.xtoidx(OptionalArgs.x);

            % Preallocate memory.
            npc0 = NpcScores("score_type", OptionalArgs.score_type, "npc_method", OptionalArgs.npc_method, "t_thresh", OptionalArgs.t_thresh, "x_names", this.x_names(xidx), "randomization_method", this.permutations.randomization_method);
            npc0.scores = zeros(1, this.n_nodes, length(xidx));
            if OptionalArgs.compute_p_values
                npc0.p_values = zeros(1, this.n_nodes, length(xidx));
            end

            % Compute scores for each node.
            if OptionalArgs.show_progress
                line_length1 = fprintf("Computing %s motion impact score\nusing %s non-parametric combining\non node ", OptionalArgs.score_type.to_string(), OptionalArgs.npc_method.to_string());
                line_length2 = fprintf("%d of %d", 1, this.n_nodes);
            end
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
        function pa = power_analysis(this)
            pa = PowerAnalysis;
            pa.nperm = this.permutations.nperm;
            pa.nboot = size(this.full_model_fit.t,1);
            pa.subsamples = this.subsamples;
            pa.x_names = this.x_names;
            pa.randomization_method = this.permutations.randomization_method;
            pa.npc_method = NpcMethod.Stouffer;
            [npc0, ~] = this.get_npc_scores("score_type", ScoreType.FalsePositive);
            pa.p_values_fp = npc0.p_values;
            [npc0, ~] = this.get_npc_scores("score_type", ScoreType.FalseNegative);
            pa.p_values_fn = npc0.p_values;
        end
    end
    methods (Static)
        function u = compute_u_values(t0, tperm, OptionalArgs)
            arguments
                t0 {mustBeNumeric}
                tperm {mustBeNumeric}
                OptionalArgs.randomization_method RandomizationMethod {mustBeScalar,mustBeNonempty,mustRandomize} = RandomizationMethod.None
                OptionalArgs.full_model_t {mustBeNumeric} = []
                OptionalArgs.score_type ScoreType = ScoreType.getDefaultValue()
                OptionalArgs.t_thresh {mustBeNumeric,mustBeScalar,mustBeNonnegative} = 2
                OptionalArgs.x_names string {mustBeVector} = []
                OptionalArgs.subsamples {mustBeUnique,mustBeInteger,mustBeVector,mustBeNonnegative,mustBeNonempty} = [0]
            end
            % Use Shaman.get_u_values() instead.
            %
            % Static method for computing u-values.  See the documentation
            % for Shaman.get_u_values(). This method is called by
            % Shaman.get_u_values() to compute a single row of u-values.

            % Validate arguments.
            if OptionalArgs.score_type ~= ScoreType.TwoSided
                assert(~isempty(OptionalArgs.full_model_t), 'Need t-values from full model to compute false positive or false negative motion impact score.');
            end
            assert(size(tperm,1) > 1, "Cannot compute u values without running any permutations.");

            % If we are computing false negative motion impact score,
            % simply flip the signs on the full model's t-values.
            if OptionalArgs.score_type == ScoreType.FalseNegative
                OptionalArgs.full_model_t = -OptionalArgs.full_model_t;
            end

            % Prellocate memory for u values.
            u = zeros(size(t0));

            % Iterate over rows in t0.
            for i=1:size(t0,1)
                % Begin by computing two-sided u-values.
                u(i,:,:,:) = sum(abs(tperm) > abs(t0(i,:,:,:))) ./ size(tperm,1);
    
                % Then compute one-sided u-values for edges above the t-value
                % threshold.
                if OptionalArgs.score_type ~= ScoreType.TwoSided
                    if size(t0,1) == size(OptionalArgs.full_model_t,1)
                        % There is one bootstrap for each row in t0.
                        ft = OptionalArgs.full_model_t(i,:,:,:);
                    else
                        % Randomly select a bootstrap for comparison.
                        j = randi(size(OptionalArgs.full_model_t,1),1,1);
                        ft = OptionalArgs.full_model_t(j,:,:,:);
                    end
                    thresh = OptionalArgs.t_thresh;
                    i_thresh = ft < -thresh;
                    u(i, i_thresh) = sum(tperm(:, i_thresh) < t0(i, i_thresh)) ./ size(tperm,1);
                    i_thresh = ft > thresh;
                    u(i, i_thresh) = sum(tperm(:, i_thresh) > t0(i, i_thresh)) ./ size(tperm,1);
                end
            end

            % Package result as a UValues object.
            u = UValues("u", u, "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh, "randomization_method", OptionalArgs.randomization_method, "x_names", OptionalArgs.x_names, "subsamples", OptionalArgs.subsamples);
        end
        function npc = compute_npc_scores(u, OptionalArgs)
            arguments
                u UValues
                OptionalArgs.npc_method NpcMethod = NpcMethod.getDefaultValue()
                OptionalArgs.show_progress logical = true
            end
            % Use Shaman.get_npc_scores() instead.
            %
            % Static method for computing non-parametric combining scores.
            % See the documentation for Shaman.get_npc_scores(). This
            % method is called by Shaman.get_npc_scores() to compute npc
            % scores from u-values.

            % Preallocate memory.
            npc = zeros(size(u.u,1), 1, size(u.u,3), size(u.u,4));   

            for i=1:size(u.u,1)
                npc(i,:,:,:) = NpcMethod.npc(squeeze(u.u(i,:,:,:)), OptionalArgs.npc_method);
            end

            % Package into an NpcScores object.
            npc = NpcScores("scores", npc, "score_type", u.score_type, "t_thresh", u.t_thresh, "npc_method", OptionalArgs.npc_method, "x_names", u.x_names, "randomization_method", u.randomization_method);
        end
    end
    methods (Access = private)
        function xidx = xtoidx(this, x)
            % x can be a cell array, string vector, or vector of indices.
            % Convert x into indices in this.x_names.
            % Make sure all the entries of x are valid.
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
            % Take a vector of nodes and a vector of edges.
            % Make sure at least one of them is empty.
            % Make sure the entries are valid.
            % If we have a vector of nodes, convert it to a vector of edges.

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
