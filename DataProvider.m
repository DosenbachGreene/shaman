classdef (Abstract) DataProvider < handle
    % Matlab class to provide neuroimaging data and associated metadata
    % such as motion/framewise displacement and non-imaging variables (e.g.
    % behavioral data). You must implement this class to feed your data
    % into Shaman. Specifically, you must implement nextData(),
    % isMoreData(), and reset(). See SimulatedDataProvider for a reference
    % implementation.
    %
    % The DataProvider functions much like an iterator over study
    % participants. Participant's data is produced one-at-a-time, which
    % helps to economize memory use. Data can be obtained in an arbitrary
    % order using nextData(), which allows the underlying implementation to
    % use asynchronous IO operations.
    methods (Abstract)
        % Return the next participant's data. Participants may be returned
        % in any order, and the order may be different with different
        % instantiations of the DataProvider or after a call to reset().
        % For example, one instance may return 3 participants in the order
        % 1, 2, 3. After calling reset it may return 1, 3, 2. Another
        % instance may return 3, 2, 1. A duplicate subject may not be
        % returned twice before a call to reset(), e.g. 1, 1, 2, 3 is not
        % allowed. Once data from all participants have been returned,
        % subsequent calls to nextData() should return an empty matrix []
        % until reset() is called.
        data = nextData(this)
        % Returns true if a call to nextData() will return more data, false
        % if there are no more participants to get data from.
        moreData = isMoreData(this)
        % Reset the state of the DataProvider. Equivalent to instantiating
        % a new data provider. After calling reset(), subsequent calls to
        % nextData() will return all of the participants again.
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