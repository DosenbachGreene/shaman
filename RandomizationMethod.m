classdef RandomizationMethod
    % Non-parametric combining method, e.g. Stouffer's Z.
    enumeration
        None
        Unconstrained
        MotionBlocks
    end
    methods
        function i = randomize(this, motion)
            arguments
                this RandomizationMethod
                motion {mustBeNumeric,mustBeVectorOrEmpty} % time x 1 vector of motion, e.g. framewise displacement
            end
            switch(this)
                case RandomizationMethod.Unconstrained
                    i = RandomizationMethod.randomize_unconstrained(motion);
                case RandomizationMethod.MotionBlocks
                    i = RandomizationMethod.randomize_motion_blocks(motion);
            end
        end
    end
    methods (Static)
        function obj = getDefaultValue()
            % Default method for non-parametric combining is Stouffer's Z score.
            obj = RandomizationMethod.MotionBlocks;
        end
        function i = randomize_none(motion)
            arguments
                motion {mustBeNumeric,mustBeVectorOrEmpty} % time x 1 vector of motion, e.g. framewise displacement
            end
            % Don't perform randomization.  Simply return the elements of
            % the timeseries vector in order.

            i = 1:length(motion);
        end
        function i = randomize_unconstrained(motion)
            arguments
                motion {mustBeNumeric,mustBeVectorOrEmpty} % time x 1 vector of motion, e.g. framewise displacement
            end
            % Randomize individual timepoints without any constraints.
            %
            % Input: time x 1 vector of motion, e.g. framewise displaement
            % Output: time x 1 vector of randomized indices

            i = randperm(length(motion));
        end
        function i = randomize_motion_blocks(motion)
            arguments
                motion {mustBeNumeric,mustBeVectorOrEmpty} % time x 1 vector of motion, e.g. framewise displacement
            end
            % Divide the data into blocks of consecutive timepoints whose
            % motion is either above or below the median, mimicking the
            % sizes of consecutive blocks in the unrandomized split half
            % model. Then randomize the order of the blocks while
            % preserving the number, sizes, and ordering of the data within
            % each block.
            %
            % Input: time x 1 vector of motion, e.g. framewise displaement
            % Output: time x 1 vector of randomized indices

            % Initialize state.
            median_motion = median(motion); % The median motion.
            motion_was_high = motion(1) >= median_motion; % True if a time point is above median motion.
            num_timepoints = 0; % Number of time points a block.
            start_idx = 1; % Starting index of a block.
            blocks = []; % Data structure to store block information.
        
            % Iterate over timepoints from start to end.
            for j = 1:length(motion)
                motion_is_high = motion(j) >= median_motion;
        
                if motion_is_high == motion_was_high
                    % Increment number of timepoints in this block.
                    num_timepoints = num_timepoints + 1;
                end
                if motion_is_high ~= motion_was_high || j == length(motion)
                    % We have reached the end of a block.
                    % Append current block information to blocks array.
                    blocks = [blocks, struct('start_idx', start_idx, 'num_timepoints', num_timepoints)];
        
                    % Reset state for next block.
                    num_timepoints = 1;
                    start_idx = j;
                end
                if motion_is_high ~= motion_was_high && j == length(motion)
                    % Add final block if it has just one timepoint.
                    blocks = [blocks, struct('start_idx', start_idx, 'num_timepoints', num_timepoints)];
                end

                % Store motion state for next iteration.
                motion_was_high = motion_is_high;
            end

            % Shuffle the order of the blocks.
            blocks = blocks(randperm(numel(blocks)));
        
            % Generate randomized timeseries indices from the blocks.
            i = zeros(1, length(motion)); % prellocate memory
            j_idx = 1; % starting position
            for j_block = 1:numel(blocks)
                start_idx = blocks(j_block).start_idx; % starting index for this block
                num_timepoints = blocks(j_block).num_timepoints; % number of timepoints in this block
                i(j_idx : j_idx+num_timepoints-1) = start_idx : (start_idx + num_timepoints - 1); % fill in indices for this block
                j_idx = j_idx + num_timepoints; % advance position in i for next iteration
            end
        end
    end
end