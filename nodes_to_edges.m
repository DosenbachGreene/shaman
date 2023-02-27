function edges = nodes_to_edges(Args)
    arguments
        Args.total_nodes {mustBeInteger,mustBePositive,mustBeScalar}
        Args.nodes {mustBeInteger,mustBePositive,mustBeVector}
    end
    % Convert a vector of nodes into a vector of edges.
    %
    % Takes two named arguments:
    %
    %     total_nodes: The total number of nodes in the parcellation.
    %     nodes: Vector of nodes to select.

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