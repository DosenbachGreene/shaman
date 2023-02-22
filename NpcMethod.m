classdef NpcMethod
    enumeration
        Stouffer
    end
    methods (Static)
        function obj = getDefaultValue()
           obj = NpcMethod.Stouffer;
        end
        function x = npc(u, npc_method)
            arguments
                u (1,:) {mustBeNumeric,mustBeNonempty} % vector of u-values to combine
                npc_method NpcMethod = NpcMethod.getDefaultValue()
            end

            % Dispatch to appropriate method.
            switch(npc_method)
                case NpcMethod.Stouffer
                    x = NpcMethod.stouffer(u);
            end
        end
        function z = stouffer(u, epsilon)
            arguments
                u (1,:) {mustBeNumeric,mustBeNonempty} % vector of u-values to combine
                epsilon {mustBeNumeric,mustBeScalar,mustBePositive} = eps; % a small number
            end
            
            % Replace u=0 or u=1 with epsilon.
            % Necessary to prevent getting Inf or -Inf for z-values.
            u(u < epsilon) = epsilon;
            u(u > 1-epsilon) = 1 - epsilon;
            
            % Compute Stouffer's Z score.
            z = sum(norminv(1 - u)) / sqrt(length(u));
        end
    end
end