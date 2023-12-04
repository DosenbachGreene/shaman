classdef PowerAnalysis < handle
    % The results of a power analysis, see Shaman.power_analysis().

    properties
        nperm double = 0 % Number of permutations at each sample size.
        nboot double = 0 % Number of un-permuted bootstraps at each sample size.
        x_names string {mustBeVector,mustBeNonempty} = [0] % Names of non-imaging variables provided by the DataProvider upon which to compute motion impact score. A separate score is computed independently for each variable. Example: ["trait1", "trait2"]
        randomization_method RandomizationMethod {mustBeScalar,mustBeNonempty,mustRandomize} = RandomizationMethod.getDefaultValue() % How were permutations randomized
        npc_method NpcMethod {mustBeScalar,mustBeNonempty} = NpcMethod.getDefaultValue() % NpcMethod, Default: Stouffer
        p_values_fp double % Matrix of omnibus p-values for false positive scores. By convention the dimensions are bootstraps (or 1 row if not permuted) x 1 x predictor variables x sample size. Default: []
        p_values_fn double % Matrix of omnibus p-values for false negative scores. By convention the dimensions are bootstraps (or 1 row if not permuted) x 1 x predictor variables x sample size. Default: []
        subsamples uint32 {mustBeUnique,mustBeInteger,mustBeVector,mustBeNonnegative,mustBeNonempty} = [0] % vector of subsample sizes; a size of 0 means no subsampling
    end

    methods
        function x = sensitivity_fp(this, p_thresh)
            % Sensitivity to false positive scores, ranging from 0 to 1.
            % p_thresh is the alpha threshold.

            % First subsample should be full, un-bootstrapped sample.
            assert(this.subsamples(1) == 0);
            % Need at least one subsample to do a power analysis.
            assert(length(this.subsamples) > 1);

            % Delegate to static method.
            x = PowerAnalysis.sensitivity(this.p_values_fp, p_thresh);
        end
        function x = sensitivity_fn(this, p_thresh)
            % Sensitivity to false negative scores, ranging from 0 to 1.
            % p_thresh is the alpha threshold.

            % First subsample should be full, un-bootstrapped sample.
            assert(this.subsamples(1) == 0);
            % Need at least one subsample to do a power analysis.
            assert(length(this.subsamples) > 1);

            % Delegate to static method.
            x = PowerAnalysis.sensitivity(this.p_values_fn, p_thresh);
        end
        function x = specificity_fp(this, p_thresh)
            % Specicifity to false positive scores, ranging from 0 to 1.
            % p_thresh is the alpha threshold.

            % First subsample should be full, un-bootstrapped sample.
            assert(this.subsamples(1) == 0);
            % Need at least one subsample to do a power analysis.
            assert(length(this.subsamples) > 1);

            % Delegate to static method.
            x = PowerAnalysis.specificity(this.p_values_fp, p_thresh);
        end
        function x = specificity_fn(this, p_thresh)
            % Specificity to false negative scores, ranging from 0 to 1.
            % p_thresh is the alpha threshold.

            % First subsample should be full, un-bootstrapped sample.
            assert(this.subsamples(1) == 0);
            % Need at least one subsample to do a power analysis.
            assert(length(this.subsamples) > 1);

            % Delegate to static method.
            x = PowerAnalysis.specificity(this.p_values_fn, p_thresh);
        end
    end
    methods (Static)
        function x = sensitivity(p_values, p_thresh)
            % Specifity, ranging from 0 to 1.
            % p_thresh is the alpha threshold.

            % First subsample should be full, un-bootstrapped sample.
            % Find true positives: significant scores in full sample.
            tp = p_values(1,:,:,1) < p_thresh;
            
            % Find significant scores in the bootstrap samples.
            bp = p_values(:,:,:,2:end) < p_thresh;

            % Find false negatives: significant scores in full sample that
            % were not significant in a bootstrap sample.
            fn = tp & ~bp;

            % Sensitivity = 1 - false negative rate.
            x = squeeze(1 - sum(fn)./size(fn,1))';
        end
        function x = specificity(p_values, p_thresh)
            % Specificity, ranging from 0 to 1.
            % p_thresh is the alpha threshold.

            % First subsample should be full, un-bootstrapped sample.
            % Find true negatives: not-significant scores in full sample.
            tn = p_values(1,:,:,1) >= p_thresh;
            
            % Find significant scores in the bootstrap samples.
            bp = p_values(:,:,:,2:end) < p_thresh;

            % Find false positives: significant scores in full sample that
            % were not significant in a bootstrap sample.
            fp = tn & bp;

            % Sensitivity = 1 - false positive rate.
            x = squeeze(1 - sum(fp)./size(fp,1))';
        end
    end
end