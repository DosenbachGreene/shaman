function edges = nodes_to_edges(Args)
    arguments
        Args.total_nodes {mustBeInteger,mustBePositive,mustBeScalar} % the total number of nodes in a parcellation
        Args.nodes (1,:) {mustBeInteger,mustBePositive} % a vector of which nodes to select
    end

    % Generate a connectivity matrix with ones in the rows and columns for
    % the desired nodes.
    mat = false(Args.total_nodes);
    mat(Args.nodes, :) = true;
    mat(:,Args.nodes) = true;

    % Vectorize the connectivity matrix.
    vec = corrmat_vectorize(mat);
    clear mat;

    % Find the ones in the vector.
    edges = find(vec);

end