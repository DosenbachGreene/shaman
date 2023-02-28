% Generate simulated data.
% cd simulate
% simulate
% cd ..

% Create a DataProvider object that points at the simulated data.
% For analyzing real data, read the comments in DataProvider.m and
% implement your own DataProvider class.
data_provider = SimulatedDataProvider("simulate");

% Get one participant's data and see what it looks like.
% data.fmri is a 1024 time points x 394 parcels matrix
% data.motion is a 1024 time points vector of framewise displacement
data = data_provider.nextData()
% By convention, data.tbl has a column `id` for participant id.
% You will see another column `trait`. This is a simulated stand-in for a
% real trait like height, weight, IQ score, etc. The table may have
% additional columns for additional traits or covariates.
data.tbl
% We can compute a correlation matrix for this participant's data and
% visualize it using imagesc().
corrmat = corr(data.fmri);
corrmat = corrmat - diag(diag(corrmat)); % remove diagonal for visualization
figure;
imagesc(corrmat);
colormap("gray");
title("Correlation Matrix for a Single Participant");
% We can extract the unique elements of the correlation matrix as a vector.
corrmat_vectorized = corrmat_vectorize(corrmat);
size(corrmat_vectorized) % 1 x 77421
% And we can reconstitute the original correlation matrix.
corrmat = corrmat_unvectorize(corrmat_vectorized);
size(corrmat) % 394 x 394
% Now let's clean up and get ready to use Shaman.
close all
data_provider.reset();
clear data corrmat corrmat_vectorized

% Create a Shaman object. At a minimum we need a DataProvider and a vector
% of traits we would like to calculate motion impact scores for.
shaman = Shaman(data_provider, ["trait"]);
% Read the documentation for the Shaman class to learn about more advanced
% features such as including covariates.
doc Shaman

% At this point we can already look at the connnectivity matrix for our
% trait. This is the connectivity matrix from the "full" model. With our
% simulated data, we can easily see the connectivity matrix is contaminated
% by motion.
figure;
imagesc(corrmat_unvectorize(shaman.full_model_fit.t));
colormap("gray");
handle = colorbar;
set(get(handle,"label"),"string","t-value");
clear handle
title("Connectivity Matrix for Trait");
% We can also visualize the connectivity matrix from the split-half model.
% This is the raw motion impact score. We can see the motion impact score
% flags the "motion" part of the conenctivity matrix without getting any of
% the "brain."
figure;
imagesc(corrmat_unvectorize(shaman.split_model_fit.t));
colormap("gray");
handle = colorbar;
set(get(handle,"label"),"string","t-value");
clear handle
title("Raw Motion Impact Score for Trait");
% Cleanup.
close all

% To calculate a aggregate motion impact score and get p-values we first
% need to run some permutations of the split-half model. To get us started
% quickly:
shaman.permutations.nperm = 32;
% Ordinarily we would want to run thousands of permutations. The good news
% is, we can run additional permutations without losing the first 32!
%
% shaman.permutations.nperm = 1000; % Run 968 more permutations.
%
% Now we can get a single number for the entire connectivity matrix.
shaman.get_scores_as_table()
% The table shows us the motion impact score, which is a Stouffer's Z
% score, and an omnibus p-value. The omnibus p-value basically tests if
% *any* of the edges in the connectivity matrix is affected by motion. If
% we have a specific hypothesis, e.g. connectivity between nodes 18 and
% 104, then we can restrict analysis to just those nodes. Now the p-value
% is no longer significant.
shaman.get_scores_as_table("nodes", [18, 104])
% Shaman computes a two-sided score by default. This tests whether motion
% inflates or obscures connectivity. If we are only interested in avoiding
% false positive inference:
shaman.get_scores_as_table("score_type", "FalsePositive")
% In this simulated example, the motion artifact causes false positive
% inference but not false negative inference.
shaman.get_scores_as_table("score_type", "FalseNegative")

% We can also look at the motion impact score by using the NpcScores class.
% This is short for non-parametric combining score, which is the same thing
% as a motion impact score.
[npc0, npc] = shaman.get_npc_scores();
% npc0 is the npc score we are interested in. The npc are the npc scores
% from the permutations. We can visualize this as a histogram. npc0 is far
% to the right of the null distribution, so it is very significant!
figure;
histogram(npc.scores);
xline(npc0.scores);
% If we look at the npc0 object we get the same information as the table,
% plus some additional metadata. The `scores` property is the motion impact
% score and the `p_values` property is the p-value.
npc0
% Cleanup.
clear npc0
close all

% For extra detail, we can look at the u-values (pseudo p-values) at each
% edge to see where (at what edges) there is a significant motion impact.
% Note that the u-values are *not* corrected for multiple comparisons.
uvalues = shaman.get_u_values();
motion_impact = shaman.split_model_fit.t;
motion_impact(uvalues.u > 0.05) = 0; % apply threshold of u = 0.05
figure;
imagesc(corrmat_unvectorize(motion_impact));
colormap("gray");
handle = colorbar;
set(get(handle,"label"),"string","t-value");
clear handle
title("Raw Motion Impact Score for Trait, u < 0.05");
% Cleanup.
close all
clear motion_impact uvalues

% Suppose we want to render motion impact score on the brain, like in
% Figure 4. Shaman can conveniently compute motion impact score for each
% node. We can see that the output has 394 scores, one for each node.
npc0 = shaman.get_scores_by_node();
size(npc0.scores)

% I hope this example has familiarized you with Shaman!  Contact Benjamin
% Kay <benjamin.kay@wustl.edu> with questions or feedback.