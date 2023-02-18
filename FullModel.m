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
            this.con = [this.con; zeros(data_provider.size_hint_participants()-1, length(this.con))];
            this.motion = [this.motion; zeros(data_provider.size_hint_participants()-1, 1)];
            this.tbl = [data.tbl; table('Size', [data_provider.size_hint_participants()-1, size(data.tbl,2)], 'VariableTypes', varfun(@class,data.tbl,'OutputFormat','cell'), 'VariableNames', data.tbl.Properties.VariableNames)];
            
            % Load the remaining participants.
            i = 2;
            while data_provider.isMoreData()
                % Display progress.
                if OptionalArgs.Progress
                    fprintf(repmat('\b',1,line_length));
                    line_length = fprintf('%d of %d', i, data_provider.size_hint_participants());
                end
                
                data = data_provider.nextData();
                this.con(i,:) = corrmat_vectorize(atanh(corr(data.fmri)));
                this.motion(i) = mean(data.motion);
                if size(this.tbl,1) >= i
                    this.tbl(i,:) = data.tbl;
                else
                    this.tbl = [this.tbl; data.tbl];
                end
                i = i + 1;
            end
            
            % Display progress.
            if OptionalArgs.Progress
                fprintf(repmat('\b',1,line_length + 23));
                fprintf('Processed %d participants.\n', i-1);
            end
        end
    end
end