classdef Data
    % Simple structure provided by a DataProvider.
    % Contains the fMRI and behavioral data for a single study participant.
    properties
        fmri % time (rows) x voxels/vertices (columns) matrix fMRI data
        motion % time x 1 vector of motion, e.g. framewise displacement
        tbl % table of covariates, by convention "id" column is study-assigned participant id as a string or number
    end
end