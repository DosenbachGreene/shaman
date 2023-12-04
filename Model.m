classdef (Abstract) Model
    % Abstract parent class for SplitModel and FullModel.
    properties (Abstract)
        con {mustBeNumeric} % Participant (rows) x node (column) connectivity matrices. Each row is a connectivity matrix that has been vectorized using corrmat_vectorize(). Unvectorize using corrmat_unvectorie().
        motion {mustBeNumeric,mustBeVectorOrEmpty} % Participant (rows) x 1 vector of motion (e.g. mean framewise displacement).
        tbl table % Participant (rows) x covariate (column) table of variables. By convention, the participant id is in column named "id".
    end
    methods
        function newmodel = resample(this, v)
            arguments
                this Model
                v {mustBeInteger,mustBeNumeric,mustBePositive,mustBeNonempty}
            end
            % Generate a new model sampling the observations indexed by v.
            % If v is a scaler, generate a new model with v observations
            % sampled randomly from the original model, with replacement.

            if isscalar(v)
                assert(v > 0)
                % Make a random sampling vector with replacement.
                v = randi(size(this.con,1), v, 1);
            end
            assert(isvector(v));
            
            newmodel = this;
            newmodel.con = newmodel.con(v,:);
            newmodel.motion = newmodel.motion(v);
            newmodel.tbl = newmodel.tbl(v,:);
        end
    end
end