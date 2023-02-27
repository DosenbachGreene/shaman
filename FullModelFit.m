classdef FullModelFit < ModelFit
    % Subclass of ModelFit that only fits a FullModel.
    methods
        function this = FullModelFit(model, x_names, OptionalArgs)
            % See documentation for the constructor of ModelFit.
            arguments
                model FullModel
                x_names string {mustBeVector,mustBeNonempty}
                OptionalArgs.intercept logical = true
                OptionalArgs.motion_covariate logical = true
                OptionalArgs.covariates string {mustBeVectorOrEmpty} = []
                OptionalArgs.show_progress logical = true
            end
            
            % Delegate to superclass constructor.
            OptionalArgs = namedargs2cell(OptionalArgs);
            this = this@ModelFit(model, x_names, OptionalArgs{:});
       end
    end
end