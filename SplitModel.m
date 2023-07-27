classdef SplitModel < Model
    % A model of split-half (as opposed to full) connectivity.
    %
    % 1. A connectivity matrix is generated from the high- and low-motion
    % halves of each participants MRI data.
    % 2. Between-participant motion is regressed out of the high- and
    % low-motion connectivity matrices.
    % 3. The difference between the residuals of the high- and low-motion
    % connectivity matrices is stored in SplitModel.con.
    %
    % Rather than splitting fMRI timeseires data into high- and low-motion
    % halves, the split can be randomized with respect to motion for the
    % purpose of permutation testing.
    properties (SetAccess=protected, GetAccess=public)
        con % See Model.con.
        motion % See Model.motion.
        tbl % See model.tbl.
        randomized logical = false % Whether split was randomized with respect to motion.
    end
    methods
        function this = SplitModel(data_provider, OptionalArgs)
            arguments
                data_provider DataProvider
                OptionalArgs.show_progress logical = true
                OptionalArgs.randomize logical = false
            end
            % Generate the model using data from the supplied DataProvider.
            %
            % Optional arguments:
            %
            % show_progress: Whether to show a progress indicator. Default: true
            % randomize:
            %     Whether to randomize the order of the split.
            %     Default: false, order the split from low to high motion.

            % Store arguments.
            this.randomized = OptionalArgs.randomize;
            
            % Rewind the data provider to its beginning.
            data_provider.reset();
            
            % Display progress.
            if OptionalArgs.show_progress
                fprintf('Processing participant 1');
                line_length = 1;
            end
            
            % Get the first participant's data.
            data = data_provider.nextData();
            
            % Compute connectivity from the data.
            [con_low, con_high, motion_low, motion_high] = this.corr(data, this.randomized);

            % Average motion.
            this.motion = mean(data.motion);

            % Preallocate memory for the remaining participants.
            size_hint = data_provider.size_hint_participants() - 1;
            if size_hint == 0
                size_hint = 1;
            end
            con_low = [con_low; zeros(size_hint - 1, length(con_low))];
            con_high = [con_high; zeros(size_hint - 1, length(con_high))];
            motion_low = [motion_low; zeros(size_hint - 1, 1)];
            motion_high = [motion_high; zeros(size_hint - 1, 1)];
            this.motion = [this.motion; zeros(size_hint - 1, 1)];
            this.tbl = [data.tbl; table('Size', [size_hint - 1, size(data.tbl,2)], 'VariableTypes', varfun(@class,data.tbl,'OutputFormat','cell'), 'VariableNames', data.tbl.Properties.VariableNames)];
            
            % Load the remaining participants.
            i = 1;
            while data_provider.isMoreData()
                % Increment participant number.
                i = i + 1;

                % Display progress.
                if OptionalArgs.show_progress
                    fprintf(repmat('\b',1,line_length));
                    line_length = fprintf('%d of %d', i, data_provider.size_hint_participants());
                end
                
                % Load data.
                data = data_provider.nextData();

                % Compute connectivity.
                [con_low(i,:), con_high(i,:), motion_low(i), motion_high(i)] = this.corr(data,this.randomized);

                % Average motion.
                this.motion(i) = mean(data.motion);

                % Append non-imaging data to the table.
                if size(this.tbl,1) >= i
                    this.tbl(i,:) = data.tbl;
                else
                    % Resize table if needed.
                    this.tbl = [this.tbl; data.tbl];
                end
            end

            % Trim data structures if size_hint was larger than the actual
            % number of participants.
            if i < size_hint
                con_low = con_low(1:i, :);
                con_high = con_high(1:i, :);
                motion_low = motion_low(1:i);
                motion_high = motion_high(1:i);
                this.motion = this.motion(1:i);
                this.tbl = this.tbl(1:i, :);
            end

            % Display progress.
            if OptionalArgs.show_progress
                fprintf(repmat('\b',1,line_length + 23));
                line_length = fprintf('Computing split half connectivity...\n');
            end

            % Regress motion out of each half of the data.
            con_low = this.residualize(con_low, motion_low);
            con_high = this.residualize(con_high, motion_high);

            % Compute the difference between high- and low-motion
            % connectivity matrices.
            this.con = con_high - con_low;

            % Display progress.
            if OptionalArgs.show_progress
                fprintf(repmat('\b',1,line_length));
                fprintf('Processed %d participants.\n', i);
            end
        end
    end
    methods (Access=private, Static)
        function [con_low, con_high, motion_low, motion_high] = corr(data, randomize)
            arguments
                data Data
                randomize logical
            end
            % Split data in half by motion or at random.
            % Compute connectivity from the low- and high-motion halves.
            % Average motion from the low- and high-motion halves.

            if randomize
                % Sort fMRI data randomly.
                i = permute_blocks(data);
            else
                % Sort fMRI data from low motion to high motion.
                [~,i] = sort(data.motion);
            end
            i_low = i(1:floor(length(i)/2));
            i_high = i(floor(length(i)/2)+1:length(i));

            % Generate connectivity matrices from low- and high-data.
            con_low = corrmat_vectorize(atanh(corr(data.fmri(i_low,:))));
            con_high = corrmat_vectorize(atanh(corr(data.fmri(i_high,:))));

            % Average motion from the low- and high-data.
            motion_low = mean(data.motion(i_low));
            motion_high = mean(data.motion(i_high));
        end

        function blocks = motion_based_blocking(data)
            arguments
                data Data
            end
        
            % Based on FD from midpoint of sorted motion array
            % Loops through timeseries and organizes blocks based on periods of high/low motion
            % Shuffles blocks to reconstruct timeseries 
            data.motion
            % Initialize variables
            fd = median(data.motion);          % Get average motion and use it as framewise displacement
            fd
            motion_category = data.motion(1) >= fd;  % Starting category based on first index
            num_timepoints = 0;             % Number of time points in a block
            block_start_idx = 1;            % Starting index of the current block
        
            % Initialize the data structure to store block information
            blocks = struct('start_idx', [], 'num_timepoints', []);
        
            % Loop through the motion data
            for i = 1:length(data.motion)
                if data.motion(i) >= fd
                    new_category = true;  % Check if it's in the same category as the current block
                else
                    new_category = false;
                end
        
                if new_category == motion_category
                    % Increment num_timepoints if still in the same category
                    num_timepoints = num_timepoints + 1;
                else
                    % We have reached the end of a block
                    % Save the current block information
                    blocks(end+1).start_idx = block_start_idx;
                    blocks(end).num_timepoints = num_timepoints;
        
                    % Reset num_timepoints and update motion_category
                    num_timepoints = 1;  % Set to 1 since we have a new block
                    motion_category = new_category;
                    block_start_idx = i;
                end
            end
        end

        function permuted_indices = permute_blocks(data)
            arguments
                data Data
            end

            % Call the motion_based_blocking function to obtain the original blocks
            blocks = motion_based_blocking(data);
        
            % Shuffle the order of the blocks to create the permuted timeseries
            num_blocks = numel(blocks);
            permuted_order = randperm(num_blocks);
            shuffled_blocks = blocks(permuted_order);
        
            % Concatenate the permuted blocks to form the permuted timeseries array
            permuted_indices = [];
            for i = 1:num_blocks
                block_idx = shuffled_blocks(i).start_idx;
                num_timepoints = shuffled_blocks(i).num_timepoints;
                permuted_indices = [permuted_indices, block_idx:block_idx + num_timepoints - 1];
            end
        end

        function resid = residualize(fmri, motion)
            arguments
                fmri {mustBeNumeric,mustBeNonempty}
                motion (:,1) {mustBeNumeric,mustBeNonempty}
            end
            % Fit the model fmri ~ 1 + motion
            % and return the intercept + the residuals.

            assert(size(fmri,1) == length(motion));
            x = [ones(length(motion),1), motion];
            b = x\fmri;
            resid = fmri - motion*b(2,:);
        end
    end
end