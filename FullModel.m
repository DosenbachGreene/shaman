classdef FullModel < Model
    % Load in data from a DataProvider to generate full (as opposed to
    % split-half) connectivity matrices and a table of covarariates.
    properties (SetAccess=protected, GetAccess=public)
        con % participant (rows) x node (column) connectivity matrices; unvectorize using corrmat_unvectorie()
        motion % participant (rows) x 1 table of motion
        tbl % participant (rows) x covariate (column) table of variables
    end
    methods
        function this = FullModel(data_provider, OptionalArgs)
            % Use the supplied DataProvider to populate the connectivity
            % matrices and table of regressors for this model.

            arguments
                data_provider DataProvider
                OptionalArgs.Progress logical = true
            end
            
            % Rewind the data provider to its beginning.
            data_provider.reset();
            
            % Display progress.
            if OptionalArgs.Progress
                fprintf('Processing participant 1');
                line_length = 1;
            end
            
            % Get the first participant's data.
            data = data_provider.nextData();
            
            % Compute connectivity from the data.
            this.con = corrmat_vectorize(atanh(corr(data.fmri)));

            % Average motion.
            this.motion = mean(data.motion);

            % Preallocate memory for the remaining participants.
            size_hint = data_provider.size_hint_participants() - 1;
            if size_hint == 0
                size_hint = 1;
            end
            this.con = [this.con; zeros(size_hint, length(this.con))];
            this.motion = [this.motion; zeros(size_hint, 1)];
            this.tbl = [data.tbl; table('Size', [size_hint, size(data.tbl,2)], 'VariableTypes', varfun(@class,data.tbl,'OutputFormat','cell'), 'VariableNames', data.tbl.Properties.VariableNames)];
            
            % Load the remaining participants.
            i = 1;
            while data_provider.isMoreData()
                % Increment participant count.
                i = i + 1;

                % Display progress.
                if OptionalArgs.Progress
                    fprintf(repmat('\b',1,line_length));
                    line_length = fprintf('%d of %d', i, data_provider.size_hint_participants());
                end
                
                % Load data.
                data = data_provider.nextData();

                % Compute connectivity.
                this.con(i,:) = corrmat_vectorize(atanh(corr(data.fmri)));

                % Average motion.
                this.motion(i) = mean(data.motion);
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
                this.con = con(1:i, :);
                this.motion = this.motion(1:i);
                this.tbl = this.tbl(1:i, :);
            end
            
            % Display progress.
            if OptionalArgs.Progress
                fprintf(repmat('\b',1,line_length + 23));
                fprintf('Processed %d participants.\n', i);
            end
        end
    end
end