classdef SplitModel < Model
    % Load in data from a DataProvider to generate a split-half difference
    % connectivity matrices and a table of covariates.
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
        con % participant (rows) x node (column) connectivity matrices; unvectorize using corrmat_unvectorie()
        motion % participant (rows) x 1 table of motion
        tbl % participant (rows) x covariate (column) table of variables
        randomized % whether split was randomized with respect to motion
    end
    methods
        function this = SplitModel(data_provider, OptionalArgs)
            % Use the supplied DataProvider to populate the connectivity
            % matrices and table of regressors for this model.

            arguments
                data_provider DataProvider
                OptionalArgs.show_progress logical = true % dispay progress indicator
                OptionalArgs.randomize = false % randomize split with respect to motion
            end
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
            con_low = [con_low; zeros(data_provider.size_hint_participants()-1, length(con_low))];
            con_high = [con_high; zeros(data_provider.size_hint_participants()-1, length(con_high))];
            motion_low = [motion_low; zeros(data_provider.size_hint_participants()-1, 1)];
            motion_high = [motion_high; zeros(data_provider.size_hint_participants()-1, 1)];
            this.motion = [this.motion; zeros(data_provider.size_hint_participants()-1, 1)];
            this.tbl = [data.tbl; table('Size', [data_provider.size_hint_participants()-1, size(data.tbl,2)], 'VariableTypes', varfun(@class,data.tbl,'OutputFormat','cell'), 'VariableNames', data.tbl.Properties.VariableNames)];
            
            % Load the remaining participants.
            i = 2;
            while data_provider.isMoreData()
                % Display progress.
                if OptionalArgs.show_progress
                    fprintf(repmat('\b',1,line_length));
                    line_length = fprintf('%d of %d', i, data_provider.size_hint_participants());
                end
                
                data = data_provider.nextData();
                [con_low(i,:), con_high(i,:), motion_low(i), motion_high(i)] = this.corr(data,this.randomized);
                this.motion(i) = mean(data.motion);
                if size(this.tbl,1) >= i
                    this.tbl(i,:) = data.tbl;
                else
                    this.tbl = [this.tbl; data.tbl];
                end
                i = i + 1;
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
                fprintf('Processed %d participants.\n', i-1);
            end
        end
    end
    methods (Access=private, Static)
        function [con_low, con_high, motion_low, motion_high] = corr(data, randomize)
            arguments
                data Data
                randomize logical
            end

            if randomize
                % Sort fMRI data randomly.
                i = randperm(length(data.motion));
            else
                % Sort fMRI data from low motion to high motion.
                [~,i] = sort(data.motion);
            end

            % Generate connectivity matrices from low- and high-data.
            i_low = i(1:floor(length(i)/2));
            i_high = i(floor(length(i)/2)+1:length(i));
            motion_low = mean(data.motion(i_low));
            motion_high = mean(data.motion(i_high));
            con_low = corrmat_vectorize(atanh(corr(data.fmri(i_low,:))));
            con_high = corrmat_vectorize(atanh(corr(data.fmri(i_high,:))));
        end
        function resid = residualize(fmri, motion)
            % Fit the model fmri ~ 1 + motion
            % and return the residuals.
            arguments
                fmri double
                motion (:,1) double
            end
            assert(size(fmri,1) == length(motion));
            x = [ones(length(motion),1), motion];
            b = x\fmri;
            resid = fmri - x*b;
        end
    end
end