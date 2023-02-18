classdef Permutations < handle
    properties (SetAccess=protected, GetAccess=public)
        data_provider % where data is sourced from
        x % names of variables used as main predictors
        null_model_t % n permutations (rows) x edges (columns) x variables in x (3rd dimension) array of t-values
        intercept % whether the model fit includes an interecept column
        covariates % cell array of variables included as covariates
        motion_covariate % whether motion was included as a covariate
    end
    properties
        nperm = 0 % number of permutations
        show_progress % whether to show a progress indicator
    end
    methods
        function this = Permutations(data_provider, x, OptionalArgs)
            arguments
                data_provider DataProvider
                x (1,:) cell {mustBeText,mustBeNonempty} % names of variables to use as main predictors
                OptionalArgs.nperm {mustBeNonnegative, mustBeInteger} = 0
                OptionalArgs.intercept logical = true % include an intercept term
                OptionalArgs.motion_covariate logical = true % include motion as a covariate
                OptionalArgs.covariates (1,:) cell {mustBeText} = {} % cell array of variable names to include as covarriates
                OptionalArgs.show_progress logical = true % display progress
            end

            % Store arguments.
            this.data_provider = data_provider;
            this.x = x;
            this.intercept = OptionalArgs.intercept;
            this.motion_covariate = OptionalArgs.motion_covariate;
            this.covariates = OptionalArgs.covariates;
            this.show_progress = OptionalArgs.show_progress;

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
            this.null_model_t = zeros(0, nedges, length(this.x));

            % Perform permutations.
            if OptionalArgs.nperm > 0
                this.nperm = OptionalArgs.N;
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
            t = zeros(nperm-this.nperm, size(this.null_model_t, 2), length(this.x));

            % Cache some variables on the function stack to make parfor
            % happier.
            data_provider = this.data_provider;
            x = this.x;
            intercept = this.intercept;
            motion_covariate = this.motion_covariate;
            covariates = this.covariates;
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
                
                % Generat a permuted split model.
                model = SplitModel(data_provider, "Randomize", true, "show_progress", false);
                % Fit model to the predictors and get t-values.

                fit = SplitModelFit(model, x, "Intercept", intercept, "motion_covariate", motion_covariate, "covariates", covariates, "show_progress", false);
                t(i, :, :) = fit.t';

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