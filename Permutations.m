classdef Permutations < handle
    properties (SetAccess=protected, GetAccess=public)
        data_provider DataProvider = NullDataProvider() % A DataProvider, the source of imaging and non-imaging (e.g. behavioral, biophysical, etc) data.
        x_names string {mustBeVector,mustBeNonempty} = [0] % Names of non-imaging variables provided by the DataProvider upon which to compute motion impact score. A separate score is computed independently for each variable. Example: ["trait1", "trait2"]
        null_model_t {mustBeNumeric} % n permutations (rows) x edges (columns) x variables in x (3rd dimension) matrix of t-values
        intercept logical = true % Whether the model fit includes an interecept column. Default: true
        covariates string {mustBeVectorOrEmpty} = [] % Names of additional non-imaging variables, besides motion and an intercept term, to include as covariates. The default is not to include any extra covariates. Example: ["cov1", "cov2"]
        motion_covariate logical = true % Whether to include motion as a covaraite. Default: true
        randomization_method RandomizationMethod {mustRandomize} = RandomizationMethod.getDefaultValue(); % How permutations are randomized.
    end
    properties
        nperm {mustBeInteger,mustBeNonnegative} = 0 % Number of permutations. Defaults to 0. Setting this property to a number X times larger than its initial value Y will keep the first Y permutations and run X-Y additional permutations. If X < Y it will generate an error.
        show_progress logical = true % Whether to show a progress indicator. Default: true
    end
    methods
        function this = Permutations(data_provider, x_names, OptionalArgs)
            arguments
                data_provider DataProvider
                x_names string {mustBeNonempty,mustBeVector}
                OptionalArgs.nperm {mustBeNonnegative, mustBeInteger} = 0
                OptionalArgs.intercept logical = true
                OptionalArgs.motion_covariate logical = true
                OptionalArgs.covariates string {mustBeVectorOrEmpty} = []
                OptionalArgs.show_progress logical = true
                OptionalArgs.randomization_method RandomizationMethod {mustRandomize} = RandomizationMethod.getDefaultValue();
            end
            % Set up a permutation test.
            %
            % See the documentation for SplitModelFit for an explanation of
            % arguments. The nperm argument defaults to zero, in which case
            % the actual permutations will be run later by setting the value of
            % the Permutations.nperm property. If you wish, you may set nperm in
            % the constructor to begin running permutations immediately.
            %
            % Note that permutations are run in a parallel for loop.

            % Store arguments.
            this.data_provider = data_provider;
            this.x_names = x_names;
            this.intercept = OptionalArgs.intercept;
            this.motion_covariate = OptionalArgs.motion_covariate;
            this.covariates = OptionalArgs.covariates;
            this.show_progress = OptionalArgs.show_progress;
            this.randomization_method = OptionalArgs.randomization_method;
            
            % Peek at the data to figure out how many edges there are.
            nedges = this.data_provider.size_voxels(); % actually # of voxels
            if nedges == 0
                data_provider.reset();
                [data, ~] = data_provider.nextData();
                assert(~isempty(data));
                nedges = size(data.fmri, 2); % actually # of voxels
                data_provider.reset();
            end
            % Infer number of edges from number of nodes/voxels.
            nedges = (nedges * nedges - nedges) / 2;

            % Initialize the matrix of permuted t-values.
            this.null_model_t = zeros(0, nedges, length(this.x_names));

            % Perform permutations.
            if OptionalArgs.nperm > 0
                this.nperm = OptionalArgs.nperm;
            end
        end
        function set.nperm(this, nperm)
            % Increase the number of permutations.
            arguments
                this
                nperm {mustBePositive, mustBeInteger}
            end
            assert(nperm >= 8, 'Must perform at least 8 permutations.')
            assert(nperm >= this.nperm, 'Cannot decrease the number of permutations.')
            
            % Allocate matrix of t-values.
            % Note that we cannot reallocate the existing matrix due to the
            % way Matlab handles memory in parfor loops.
            t = zeros(nperm-this.nperm, size(this.null_model_t, 2), length(this.x_names));

            % Cache some variables on the function stack to make parfor
            % happier.
            data_provider = this.data_provider;
            x = this.x_names;
            intercept = this.intercept;
            motion_covariate = this.motion_covariate;
            covariates = this.covariates;
            randomization_method = this.randomization_method;
            show_progress = this.show_progress;

            % Start parallel pool.
            p = gcp('nocreate');
            if isempty(p)
                parpool;
            end

            % Show progress indicator
            running = 0;
            finished = this.nperm;
            line_length = 0;
            function update_progress(status)
                if strcmp(status, 'start')
                    running = running + 1;
                elseif strcmp(status, 'finish')
                    finished = finished + 1;
                    running = running - 1;
                end
                fprintf(repmat('\b', 1, line_length));
                line_length = fprintf('Running permutations on %d workers.\nFinished %d of %d permutations.', running, finished, nperm);
            end
            if show_progress
                data_queue = parallel.pool.DataQueue;
                afterEach(data_queue, @update_progress);
            end

            % Perform the permutations.
            parfor i=1:(nperm-this.nperm)
                % Update progress indicator.
                if show_progress
                    send(data_queue,'start');
                end

                % Make a copy of the DataProvider to prevent a race
                % condition that could arise from multiple workers
                % operating on the same handle object.
                worker_local_data_provider = copy(data_provider);
                
                % Generate a permuted split model.
                model = SplitModel(worker_local_data_provider, "randomization_method", randomization_method, "show_progress", false);

                % Fit model to the predictors and get t-values.
                fit = SplitModelFit(model, x, "Intercept", intercept, "motion_covariate", motion_covariate, "covariates", covariates, "show_progress", false);
                t(i, :, :) = fit.t;

                % Update progress indicator.
                if show_progress
                    send(data_queue, 'finish');
                end
            end

            if show_progress
                fprintf(repmat('\b', 1, line_length));
                fprintf('Finished %d permutations.\n', nperm);
            end

            % Update the matrix of t-values.
            this.null_model_t = [this.null_model_t; t];

            % Finally, change this.nperm.
            this.nperm = nperm;
        end
    end
end