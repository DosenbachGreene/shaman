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
        permutations Permutations % Permutations object.  Set shaman.permutations.nperm = 1000 to run a thousand (or whatever number you want) permutations.
        intercept logical = true % Whether to model an intercept term. Default: true
        motion_covariate logical = true % Whether to model motion as a covariate. Default: true
        covariates string {mustBeVectorOrEmpty} % Names of additional non-imaging variables, besides motion and an intercept term, to include as covariates. The default is not to include any extra covariates. Example: ["cov1", "cov2"]
        n_nodes {mustBeInteger,mustBePositive} % Number of nodes (i.e. voxels, vertices, regions, parcels) in the model.
        n_edges {mustBeInteger,mustBePositive} % Number of edges (i.e. pairwise connections) in the model.
    end
    properties
        show_progress logical = true % Whether to show a progress indicator for operations that could take a long time. Default: true
    end
    methods
        function this = Shaman(data_provider, x_names, OptionalArgs)
            arguments
                data_provider DataProvider
                x_names string {mustBeVector,mustBeNonempty}
                OptionalArgs.nperm {mustBeInteger,mustBeNonnegative} = 0
                OptionalArgs.intercept logical = true
                OptionalArgs.motion_covariate logical = true
                OptionalArgs.covariates string {mustBeVectorOrEmpty} = []
                OptionalArgs.randomization_method RandomizationMethod {mustRandomize} = RandomizationMethod.getDefaultValue();
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
            % show_progress: Whether to show a progress indicator. Default: true

            % Store arguments in self.
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
            model = FullModel(this.data_provider, "show_progress", this.show_progress);
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
            this.permutations = Permutations(this.data_provider, this.x_names, "covariates", this.covariates, "intercept", this.intercept, "motion_covariate", this.motion_covariate, "randomization_method", OptionalArgs.randomization_method, "show_progress", this.show_progress);

            % Perform permutations.
            if OptionalArgs.nperm > 0
                this.permutations.nperm = OptionalArgs.nperm;
            end

            % Store number of nodes and edges.
            this.n_edges = size(this.full_model_fit.t,2);
            this.n_nodes = (1 + sqrt(1 + 8*this.n_edges)) / 2;
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

            % Preallocate memory for u-values.
            u0 = zeros(1, length(edges), length(xidx));
            if nargout > 1
                u = zeros(this.permutations.nperm, length(edges), length(xidx));
            end

            % Cache some variables on the function stack to make parfor
            % happier.
            ntraits = length(xidx);
            nperm = this.permutations.nperm;
            nedges = length(edges);
            score_type = OptionalArgs.score_type;
            t_thresh = OptionalArgs.t_thresh;
            full_model_t = this.full_model_fit.t(:, edges, :);
            split_model_t = this.split_model_fit.t(:, edges, :);
            null_model_t = this.permutations.null_model_t(:, edges, :);
            randomization_method = this.permutations.randomization_method;
            show_progress = OptionalArgs.show_progress;
            nout = nargout;

            % Start parallel pool.
            p = gcp('nocreate');
            if isempty(p)
                parpool;
            end

            % Show progress indicator
            running = 0;
            finished = 0;
            line_length = 0;
            function update_progress(status)
                if strcmp(status, 'start')
                    running = running + 1;
                elseif strcmp(status, 'finish')
                    finished = finished + 1;
                    running = running - 1;
                end
                fprintf(repmat('\b', 1, line_length));
                line_length = fprintf('Computing %s u-values on %d workers.\nFinished %d of %d traits.', score_type.to_string(), running, finished, ntraits);
            end
            data_queue = parallel.pool.DataQueue; % For unclear reasons, Matlab requires data_queue to exist even if show_progress is false.
            if show_progress
                afterEach(data_queue, @update_progress);
            end

            % Iterate over traits and compute u-values.
            parfor i=1:ntraits
                % Update progress indicator.
                if show_progress
                    send(data_queue,'start');
                end

                % Compute u-values for the not-permuted model.
                u0(1,:,i) = Shaman.compute_u_values(split_model_t(1, :, xidx(i)), null_model_t(:,:,xidx(i)), randomization_method, "full_model_t", full_model_t(1,:,xidx(i)), "score_type", score_type, "t_thresh", t_thresh).u;
                if nout > 1
                    % Compute u-values for each permutation.
                    u_ = zeros(nperm, nedges); % assign to thread-local variable to make parfor happy
                    for j=1:nperm
                        u_(j,:) = Shaman.compute_u_values(null_model_t(j,:,xidx(i)), null_model_t(:,:,xidx(i)), randomization_method, "full_model_t", full_model_t(1,:,xidx(i)), "score_type", score_type, "t_thresh", t_thresh).u;
                    end
                    u(:,:,i) = u_;
                end

                % Update progress indicator.
                if show_progress
                    send(data_queue,'finish');
                end
            end

            % Package results as UValues objects.
            u0 = UValues("u", u0, "score_type", OptionalArgs.score_type, "x_names", this.x_names(xidx), "t_thresh", OptionalArgs.t_thresh, "randomization_method", this.permutations.randomization_method);
            if nargout > 1
                u = UValues("u", u, "score_type", OptionalArgs.score_type, "x_names", this.x_names(xidx), "t_thresh", OptionalArgs.t_thresh, "randomization_method", this.permutations.randomization_method);
            end
            if OptionalArgs.show_progress
                fprintf(repmat('\b', 1, line_length));
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

            % Preallocate memory.
            npc0 = zeros(1, 1, length(xidx));
            if nargout > 1
                p_values = zeros(1,length(xidx));
                npc = zeros(this.permutations.nperm, 1, length(xidx));
            end

            % Cache some variables on the function stack to make parfor
            % happier.
            ntraits = length(xidx);
            nperm = this.permutations.nperm;
            nedges = length(edges);
            score_type = OptionalArgs.score_type;
            npc_method = OptionalArgs.npc_method;
            t_thresh = OptionalArgs.t_thresh;
            full_model_t = this.full_model_fit.t(:, edges, :);
            split_model_t = this.split_model_fit.t(:, edges, :);
            null_model_t = this.permutations.null_model_t(:, edges, :);
            randomization_method = this.permutations.randomization_method;
            x_names = this.x_names;
            show_progress = OptionalArgs.show_progress;
            nout = nargout;

            % Start parallel pool.
            p = gcp('nocreate');
            if isempty(p)
                parpool;
            end

            % Show progress indicator
            running = 0;
            finished = 0;
            line_length = 0;
            function update_progress(status)
                if strcmp(status, 'start')
                    running = running + 1;
                elseif strcmp(status, 'finish')
                    finished = finished + 1;
                    running = running - 1;
                end
                fprintf(repmat('\b', 1, line_length));
                line_length = fprintf('Performing %s non-parametric combining on %d workers.\nFinished %d of %d traits.', score_type.to_string(), running, finished, ntraits);
            end
            data_queue = parallel.pool.DataQueue; % For unclear reasons, Matlab requires data_queue to exist even if show_progress is false.
            if show_progress
                afterEach(data_queue, @update_progress);
            end

            % Iterate over xidx and compute npc values.
            % This is more memory-efficient than delegating to
            % get_u_values() because we do not need to instantiate the
            % entire u matrix for all elements in xidx in memory at once.
            parfor i_xidx = 1:ntraits
                % Update progress indicator.
                if show_progress
                    send(data_queue,'start');
                end

                % Compute u-values for the not-permuted model.
                u0_i = Shaman.compute_u_values(split_model_t(1, :, xidx(i_xidx)), null_model_t(:,:,xidx(i_xidx)), randomization_method, "full_model_t", full_model_t(1,:,xidx(i_xidx)), "score_type", score_type, "t_thresh", t_thresh);
                u_i = 0; % suppress MATLAB:mir_warning_maybe_uninitialized_temporary for u_i in if statement below
                if nout > 1
                    % Compute u-values for each permutation.
                    u_i = zeros(nperm, nedges); % assign to thread-local variable to make parfor happy
                    for j_perm=1:nperm
                        u_i(j_perm,:) = Shaman.compute_u_values(null_model_t(j_perm,:,xidx(i_xidx)), null_model_t(:,:,xidx(i_xidx)), randomization_method, "full_model_t", full_model_t(1,:,xidx(i_xidx)), "score_type", score_type, "t_thresh", t_thresh).u;
                    end
                    u_i = UValues("u", u_i, "score_type", score_type, "x_names", x_names(xidx(i_xidx)), "t_thresh", t_thresh, "randomization_method", randomization_method);
                end
                
                % Convert u-values to npc scores.
                if nout < 2
                    npc0(:,:,i_xidx) = Shaman.compute_npc_scores(u0_i, "npc_method", npc_method, "show_progress",false).scores;
                else
                    npc0_i = Shaman.compute_npc_scores(u0_i, "npc_method", npc_method, "show_progress", false);
                    npc_i = Shaman.compute_npc_scores(u_i, "npc_method", npc_method, "show_progress", false);
                    npc0_i.compute_p_values(npc_i);
                    npc0(:,:,i_xidx) = npc0_i.scores;
                    p_values(:,i_xidx) = npc0_i.p_values;
                    npc(:,:,i_xidx) = npc_i.scores;
                end

                % Update progress indicator.
                if show_progress
                    send(data_queue,'finish');
                end
            end

            % Package up results into NpcScores object(s).
            npc0 = NpcScores("scores", npc0, "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh, "npc_method", OptionalArgs.npc_method, "x_names", this.x_names(xidx), "randomization_method", this.permutations.randomization_method);
            if nargout > 1
                npc0.p_values = p_values;
                npc = NpcScores("scores", npc, "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh, "npc_method", OptionalArgs.npc_method, "x_names", this.x_names(xidx), "randomization_method", this.permutations.randomization_method);
            end

            % Final progress message.
            if OptionalArgs.show_progress
                fprintf(repmat('\b', 1, line_length));
                fprintf("Performed %s non-parametric combining on %d variables.\n", OptionalArgs.npc_method.to_string(), length(xidx));
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
    end
    methods (Static)
        function u = compute_u_values(t0, tperm, randomization_method, OptionalArgs)
            arguments
                t0 {mustBeNumeric,mustBeVector}
                tperm {mustBeNumeric,ismatrix}
                randomization_method RandomizationMethod {mustBeScalar,mustBeNonempty,mustRandomize}
                OptionalArgs.full_model_t {mustBeNumeric,mustBeVectorOrEmpty} = []
                OptionalArgs.score_type ScoreType = ScoreType.getDefaultValue()
                OptionalArgs.t_thresh {mustBeNumeric,mustBeScalar,mustBeNonnegative} = 2
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
            u = UValues("u", u, "score_type", OptionalArgs.score_type, "t_thresh", OptionalArgs.t_thresh, "randomization_method", randomization_method);
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
            npc = zeros(size(u.u,1), 1, size(u.u,3));

            % Perform non-parametric on each row of u (across it's second
            % dimension) separately for each variable in u (3rd dimension).
            if OptionalArgs.show_progress
                line_length1 = fprintf("Performing %s non-parametric combining on variable ", OptionalArgs.npc_method.to_string());
                line_length2 = fprintf("%d of %d", 1, size(u.u,3));
            end
            for i=1:size(u.u,3) % iterate over variables
                if OptionalArgs.show_progress
                    fprintf(repmat('\b', 1, line_length2));
                    line_length2 = fprintf("%d of %d", i, size(u.u,3));
                end
                for j=1:size(u.u,1) % iterate over rows
                    % Here is where we actually perform the non-parametric
                    % combining.
                    npc(j,1,i) = NpcMethod.npc(u.u(j,:,i), OptionalArgs.npc_method);
                end
            end
            if OptionalArgs.show_progress
                fprintf(repmat('\b', 1, line_length1 + line_length2));
                fprintf("Performed %s non-parametric combining on %d variables.\n", OptionalArgs.npc_method.to_string(), size(u.u,3));
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