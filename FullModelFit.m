classdef FullModelFit < ModelFit
    % Subclass of ModelFit that specificall fits a FullModel.
    % See also: ModelFit, SplitModelFit
    methods
        function this = FullModelFit(model, x_names, OptionalArgs)
            % Fit parameters to model using x as the main predictor.
            %
            % Example: Fit the model using the column 'foo' from model.tbl
            % as the main predictor.  Use an intercept, motion, and the
            % column 'bar' from model.tbl as covariates.
            %   fit_model = ModelFit(model, 'foo', {'bar'});
            arguments
                model FullModel
                x_names (1,:) string {mustBeNonempty} % cell array of names of columns in model.tbl to use as predictor
                OptionalArgs.intercept logical = true % include an intercept term
                OptionalArgs.motion_covariate logical = true % include motion as a covariate
                OptionalArgs.covariates string {mustBeVectorOrEmpty} = []
                OptionalArgs.show_progress logical = true
            end
            
            % Delegate to superclass constructor.
            OptionalArgs = namedargs2cell(OptionalArgs);
            this = this@ModelFit(model, x_names, OptionalArgs{:});
       end
    end
end