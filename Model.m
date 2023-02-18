classdef (Abstract) Model
    % Abstract parent class for SplitModel and FullModel.
    properties (Abstract, SetAccess=protected, GetAccess=public)
        con % participant (rows) x node (column) connectivity matrices; unvectorize using corrmat_unvectorie()
        motion % participant (rows) x 1 table of motion
        tbl % participant (rows) x covariate (column) table of variables
    end
end