classdef (Abstract) DataProvider < handle
    % Matlab class to provide neuroimaging data and associated metadata
    % such as motion/framewise displacement and non-imaging variables (e.g.
    % behavioral data). You must implement this class to feed your data
    % into Shaman. See SimulatedDataProvider for a reference
    % implementation.
    %
    % The DataProvider functions much like an iterator over study
    % participants. Participant's data is produced one-at-a-time, which
    % helps to economize memory use. Data can be obtained in an arbitrary
    % order using nextData(), which allows the underlying implementation to
    % use asynchronous IO operations. Data can also be retrieved for a
    % specific participant with dataAt().
    methods (Abstract)
        % Return the data for the participant identified by index.
        % Index can be 1 through n_participants.
        [data] = dataAt(this, index)
        % Return the next participant's data along with the index of that
        % participant. Each sequential call returns data from a different
        % participant. Sequential calls can return data in any order (not
        % necessarily 1, 2, 3...). Once data from every participant has
        % been returned, subsequent calls to nextData() will return an
        % empty data matrix [] and an index of 0, and isMoreData() will
        % return false.
        [data, index] = nextData(this)
        % Returns true if a call to nextData() will return more data, false
        % if there are no more participants to get data from.
        moreData = isMoreData(this)
        % Resets the data counter so that isMoreData() will return true and
        % the next call to nextData() will return data from the beginning
        % again.
        reset(this)
    end
    methods
        function n_participants = size_hint_participants(this)
            % Size hint.  The actual returned sizes may be zero if unknown,
            % or it may be inaccurate (actual size may be < or > the hint).
            % n_participants = number of participants
            n_participants = 0;
        end
        function n_voxels = size_voxels(this)
            % Returns number of voxels, parcels, nodes, etc; size(data,2).
            % Unlike size_hint_participants(), the number of voxels
            % returned must be accurate (cannot be an educated guess), or
            % else return 0 if the number of voxels is unknown.
            n_voxels = 0;
        end
    end
end