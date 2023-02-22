classdef ModelFit
    % Perform regression to fit b-values (regression coefficients) and
    % t-values to a model.
    properties (SetAccess=protected, GetAccess=public)
        x % list of predictors
        b % predictor x edges (columns) vector of beta values (regression coefficients) for x
        t % predictor x edges (columns) vector of t-values for x
    end
    methods
        function this = ModelFit(model, x, OptionalArgs)
            % Fit parameters to model using x as the main predictor.
            %
            % Example: Fit the model using the column 'foo' from model.tbl
            % as the main predictor.  Use an intercept, motion, and the
            % column 'bar' from model.tbl as covariates.
            %   fit_model = ModelFit(model, 'foo', {'bar'});
            arguments
                model Model
                x (1,:) cell {mustBeText,mustBeNonempty} % cell array of names of columns in model.tbl to use as predictor
                OptionalArgs.intercept logical = true % include an intercept term
                OptionalArgs.motion_covariate logical = true % include motion as a covariate
                OptionalArgs.covariates (1,:) cell {mustBeText} = {} % cell array of column names in model.tbl to include as covarriates
                OptionalArgs.show_progress logical = true % show progres indicator
            end

            % Store list of predictors.
            this.x = x;
            
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
            this.b = zeros(1, size(model.con,2), length(x));
            this.t = zeros(1, size(model.con,2), length(x));

            % Progress indicator.
            if OptionalArgs.show_progress
                fprintf('Fitting variable ');
                line_length = fprintf('1 of %d', length(this.x));
            end

            % Fit the model for each predictor in x.
            for i=1:length(x)
                % Progress indicator.
                if OptionalArgs.show_progress
                    fprintf(repmat('\b', 1, line_length));
                    line_length = fprintf('%d of %d', i, length(this.x));
                end

                % Perform regression.
                [this.b(1,:,i), this.t(1,:,i)] = ModelFit.regress(make_design_matrix(this.x{i}), model.con);
            end

            % Progress indicator.
            if OptionalArgs.show_progress
                fprintf(repmat('\b', 1, line_length + 17))
                fprintf('Fitted %d variables.\n', length(this.x));
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