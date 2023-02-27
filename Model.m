classdef (Abstract) Model
    % Abstract parent class for SplitModel and FullModel.
    properties (Abstract, SetAccess=protected, GetAccess=public)
        con {mustBeNumeric} % Participant (rows) x node (column) connectivity matrices. Each row is a connectivity matrix that has been vectorized using corrmat_vectorize(). Unvectorize using corrmat_unvectorie().
        motion {mustBeNumeric,mustBeVectorOrEmpty} % Participant (rows) x 1 vector of motion (e.g. mean framewise displacement).
        tbl table % Participant (rows) x covariate (column) table of variables. By convention, the participant id is in column named "id".
    end
end