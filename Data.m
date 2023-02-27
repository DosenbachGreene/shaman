classdef Data
    % The fMRI and non-imaging (e.g. behavioral) data for a single participant.
    %
    % See also: DataProvider
    properties
        fmri {mustBeNumeric} % time (rows) x voxels/vertices (columns) matrix of fMRI data
        motion {mustBeNumeric,mustBeVectorOrEmpty} % time x 1 vector of motion, e.g. framewise displacement
        tbl table % table of covariates, by convention "id" column is study-assigned participant id as a string or number
    end
end