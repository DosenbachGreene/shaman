classdef ModelFit
    % Perform regression to fit b-values (regression coefficients) and
    % t-values to a model.
    properties (SetAccess=protected, GetAccess=public)
        x_names string {mustBeVector,mustBeNonempty} = [0] % names of predictors, in the same order as the 3rd dimension of b and t
        b {mustBeNumeric} % 1 (row) x edges (columns) x predictors matrix of beta values
        t {mustBeNumeric} % 1 (row) x edges (columns) x predictors matrix t-values
        covariates string {mustBeVectorOrEmpty} = [] % names of covariates
        intercept logical = true % whether the model includes an intercept column
        motion_covariate logical = true % whether the model includes motion as a covariate
    end
    methods
        function this = ModelFit(model, x_names, OptionalArgs)
            % Fit parameters to model using x as the main predictor.
            %
            % Example: Fit the model using the column 'foo' from model.tbl
            % as the main predictor.  Use an intercept, motion, and the
            % column 'bar' from model.tbl as covariates.
            %   fit_model = ModelFit(model, 'foo', {'bar'});
            arguments
                model Model
                x_names (1,:) string {mustBeNonempty} % array of names of columns in model.tbl to use as predictor
                OptionalArgs.intercept logical = true % include an intercept term
                OptionalArgs.motion_covariate logical = true % include motion as a covariate
                OptionalArgs.covariates string {mustBeVectorOrEmpty} = []
                OptionalArgs.show_progress logical = true
            end

            % Store list of predictors.
            this.x_names = x_names;
            this.covariates = OptionalArgs.covariates;
            this.intercept = OptionalArgs.intercept;
            this.motion_covariate = OptionalArgs.motion_covariate;
            
            % Closure to make the design matrices.
            function X = make_design_matrix(xi)
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
            this.t = zeros(1, size(model.con,2), length(x_names));

            % Progress indicator.
            if OptionalArgs.show_progress
                fprintf('Fitting variable ');
                line_length = fprintf('1 of %d', length(this.x_names));
            end

            % Fit the model for each predictor in x.
            for i=1:length(this.x_names)
                % Progress indicator.
                if OptionalArgs.show_progress
                    fprintf(repmat('\b', 1, line_length));
                    line_length = fprintf('%d of %d', i, length(this.x_names));
                end

                % Perform regression.
                [this.b(1,:,i), this.t(1,:,i)] = ModelFit.regress(make_design_matrix(this.x_names(i)), model.con);
            end

            % Progress indicator.
            if OptionalArgs.show_progress
                fprintf(repmat('\b', 1, line_length + 17))
                fprintf('Fitted %d variables.\n', length(this.x_names));
            end
        end
    end
    methods(Static, Access=private)
        function [b, t] = regress(x, y)
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