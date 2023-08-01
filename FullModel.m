classdef FullModel < Model
    % A model of full (as opposed to split-half) connectivity.
    properties (SetAccess=protected, GetAccess=public)
        con % See Model.con.
        motion % See Model.motion.
        tbl % See model.tbl.
    end
    methods
        function this = FullModel(data_provider, OptionalArgs)
            arguments
                data_provider DataProvider
                OptionalArgs.show_progress logical = true
            end
            % Generate the model using data from the supplied DataProvider.
            %
            % Optional arguments:
            %
            % show_progress: Whether to show a progress indicator. Default: true
            
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
            % Workaround Matlab case # 06060506 by defering assignment to
            % this until after con is computed.
            con = atanh(corrmat_vectorize(corr(data.fmri)));

            % Average motion.
            this.motion = mean(data.motion);

            % Preallocate memory for the remaining participants.
            size_hint = data_provider.size_hint_participants() - 1;
            if size_hint == 0
                size_hint = 1;
            end
            con = [con; zeros(size_hint, length(con))];
            this.motion = [this.motion; zeros(size_hint, 1)];
            this.tbl = [data.tbl; table('Size', [size_hint, size(data.tbl,2)], 'VariableTypes', varfun(@class,data.tbl,'OutputFormat','cell'), 'VariableNames', data.tbl.Properties.VariableNames)];
            
            % Load the remaining participants.
            i = 1;
            while data_provider.isMoreData()
                % Increment participant count.
                i = i + 1;

                % Display progress.
                if OptionalArgs.show_progress
                    fprintf(repmat('\b',1,line_length));
                    line_length = fprintf('%d of %d', i, data_provider.size_hint_participants());
                end
                
                % Load data.
                data = data_provider.nextData();

                % Compute connectivity.
                con(i,:) = corrmat_vectorize(atanh(corr(data.fmri)));

                % Average motion.
                this.motion(i) = mean(data.motion);

                % Skip participants with NaN values.
                if anynan(con(i,:)) || anynan(this.motion(i))
                    i = i - 1;
                    continue;
                end

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
                con = con(1:i, :);
                this.motion = this.motion(1:i);
                this.tbl = this.tbl(1:i, :);
            end

            % Workaround Matlab case # 06060506 by defering assignment to
            % this until after con is computed.
            this.con = con;

            % Display progress.
            if OptionalArgs.show_progress
                fprintf(repmat('\b',1,line_length + 23));
                fprintf('Processed %d participants.\n', i);
            end
        end
    end
end