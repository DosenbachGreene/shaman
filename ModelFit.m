classdef ModelFit
    % Regress a model to fit its beta- and t-values.
    %
    % A model with one main predictor and zero or more covariates is fit to
    % obtain the beta-values (regression coefficients) and t-values for the main
    % predictor.
    %
    % A ModelFit can fit any class that inherits from Model, i.e. a FullModel or
    % a SplitModel. If you need a type of ModelFit that fits a specific subclass
    % of model (so that you don't fit the wrong subclass by accident) then see
    % FullModelFit and SplitModelFit.
    properties (SetAccess=protected, GetAccess=public)
        x_names string {mustBeVector,mustBeNonempty} = [0] % names of predictors, in the same order as the 3rd dimension of b and t
        b {mustBeNumeric} % 1 (row) x edges (columns) x predictors matrix of beta values
        t {mustBeNumeric} % bootstrap x edges (columns) x predictors x subsamples matrix t-values
        covariates string {mustBeVectorOrEmpty} = [] % names of covariates
        intercept logical = true % whether the model includes an intercept column
        motion_covariate logical = true % whether the model includes motion as a covariate
        n_participants uint32 % number of participants (rows) in the fitted model (without subsampling)
        subsamples uint32 {mustBeUnique,mustBeInteger,mustBeVector,mustBeNonnegative,mustBeNonempty} = [0] % vector of subsample sizes; a size of 0 means no subsampling
    end
    methods
        function this = ModelFit(model, x_names, OptionalArgs)
            % Fit parameters to model using x as the main predictor.
            %
            % model: An object of the Model class.
            % x_names:
            %     A string vector of the names of one or more predictors in
            %     model.tbl
            %
            % Optional arguments:
            %
            % intercept: Whether to include an intercept column. Default: true
            % motion_covariate:
            %     Whether to include model.motion as a covariate. Default: true
            % covariates:
            %     A string vector of zero or more columns in model.tbl to
            %     include as covariates. Default: [] (no extra covariates).
            % show_progress: Show a progress indicator. Default: true
            arguments
                model Model
                x_names string {mustBeVector,mustBeNonempty}
                OptionalArgs.intercept logical = true
                OptionalArgs.motion_covariate logical = true
                OptionalArgs.covariates string {mustBeVectorOrEmpty} = []
                OptionalArgs.nboot {mustBeNonnegative,mustBeInteger,mustBeScalar} = 1
                OptionalArgs.subsamples {mustBeUnique,mustBeInteger,mustBeVector,mustBeNonnegative,mustBeNonempty} = [0]
                OptionalArgs.show_progress logical = true
            end

            % Store arguments.
            this.x_names = x_names;
            this.covariates = OptionalArgs.covariates;
            this.intercept = OptionalArgs.intercept;
            this.motion_covariate = OptionalArgs.motion_covariate;
            this.n_participants = size(model.con,1);
            this.subsamples = unique(OptionalArgs.subsamples); % sort
            assert(all(this.subsamples < this.n_participants));
            nboot = OptionalArgs.nboot; % can be inferred from ModelFit.t
            
            % Closure to make the design matrices.
            function X = make_design_matrix(model, xi)
                X = model.tbl.(xi); % predictor of interest in 1st column
                X = [X, table2array(model.tbl(:, OptionalArgs.covariates))]; % covariates
                if OptionalArgs.motion_covariate
                    X = [X, model.motion]; % motion as a covariate
                end
                if OptionalArgs.intercept
                    X = [X, ones(size(X,1),1)]; % intercept column
                end
            end

            % Allocate memory.
            this.b = zeros(1, size(model.con,2), length(x_names));
            this.t = zeros(nboot, size(model.con,2), length(x_names), length(this.subsamples));

            % Fit model for each subsample.
            % Progress indicator.
            if OptionalArgs.show_progress
                line_length_1 = fprintf('Fitting subsample ');
                line_length_2 = 0;
                line_length_3 = 0;
                line_length_4 = 0;
                line_length_5 = 0;
                line_length_6 = 0;
            end
            for i_subsample = 1:length(this.subsamples)
                subsample_size = this.subsamples(i_subsample);

                % Progress indicator.
                if OptionalArgs.show_progress
                    fprintf(repmat('\b', 1, line_length_2 + line_length_3 + line_length_4 + line_length_5 + line_length_6));
                    line_length_2 = fprintf('%d of %d', i_subsample, length(this.subsamples));
                    line_length_3 = fprintf(", bootstrap ");
                    line_length_4 = 0;
                    line_length_5 = 0;
                    line_length_6 = 0;
                end
                for i_boot = 1:nboot
                    if OptionalArgs.show_progress
                        fprintf(repmat('\b', 1, line_length_4 + line_length_5 + line_length_6));
                        line_length_4 = fprintf('%d of %d', i_boot, nboot);
                        line_length_5 = fprintf(', variable ');
                        line_length_6 = 0;
                    end

                    if subsample_size ~= 0
                        % Randomly subsample the model with replacement.
                        submodel = model.resample(subsample_size);
                    end

                    % Fit the model for each predictor in x.
                    for i_x=1:length(this.x_names)
                        % Progress indicator.
                        if OptionalArgs.show_progress
                            fprintf(repmat('\b', 1, line_length_6));
                            line_length_6 = fprintf('%d of %d', i_x, length(this.x_names));
                        end
        
                        % Perform regression.
                        if subsample_size == 0
                            [this.b(1,:,i_x), this.t(1,:,i_x,i_subsample)] = ModelFit.regress(make_design_matrix(model, this.x_names(i_x)), model.con);
                        else
                            [~, this.t(i_boot,:,i_x,i_subsample)] = ModelFit.regress(make_design_matrix(submodel, this.x_names(i_x)), submodel.con);
                        end
                    end

                    if subsample_size == 0
                        % Don't bootstrap for special size of zero = full
                        % sample.
                        this.t(:,:,:,i_subsample) = repmat(this.t(1,:,:,i_subsample), size(this.t,1), 1);
                        break;
                    end
                end
            end
            if OptionalArgs.show_progress
                fprintf(repmat('\b', 1, line_length_1 + line_length_2 + line_length_3 + line_length_4 + line_length_5 + line_length_6));
                fprintf('Fitted %d subsamples, %d bootstraps, and %d variables.\n', length(this.subsamples), nboot, length(this.x_names));
            end
        end

    end
    methods(Static, Access=private)
        function [b, t] = regress(x, y)
            % Perform massive univariate regression with observed variable y
            % and design matrix x. If y has multiple columns then each column
            % of y is regressed independently. Returns beta- and t-values for
            % the first column in the design matrix x.
            
            % Strip design matrix rows with NaN out of x and y.
            [nan_rows, ~] = find(isnan(x));
            nan_rows = unique(nan_rows);
            x(nan_rows,:) = [];
            y(nan_rows,:) = [];
            
            % Do ordinary least squares regression.
            xpinv = pinv(x); % Moore-Penrose pseudoinverse
            b = xpinv * y; % beta coefficients
            mse = y - x*b; % residuals
            mse = sum(mse.*mse) ./ (size(x,1) - size(x,2) - 1); % mean squared error
            b = b(1,:); % we only need to beta values for the first column in b predictor
            se = sum(xpinv .* xpinv, 2); % equivalent to se = diag(xpinv*xpinv'), gets variance associated with predictors
            se = sqrt(se(1) .* mse); % standard error for the first predictor/column in x
            clear mse xpinv
            t = b ./ se; % t-values for the main predictor
        end
    end
end