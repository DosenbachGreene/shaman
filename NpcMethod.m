classdef NpcMethod
    % Non-parametric combining method, e.g. Stouffer's Z.
    enumeration
        Stouffer
    end
    methods
        function str = to_string(this)
            switch(this)
                case "Stouffer"
                    str = "Stouffer's Z";
                otherwise
                    error("Unknown ScoreType")
            end
        end
    end
    methods (Static)
        function obj = getDefaultValue()
            % Default method for non-parametric combining is Stouffer's Z score.
            obj = NpcMethod.Stouffer;
        end
        function x = npc(u, npc_method)
            arguments
                u {mustBeNumeric,mustBeVector,mustBeNonempty}
                npc_method NpcMethod = NpcMethod.getDefaultValue()
            end
            % Perform non-parametric combining on u-values with method.
            %
            % u: A vector of u-values to combine.
            % npc_method: An NpcMethod to use to combine the u-values.

            % Dispatch to appropriate method.
            switch(npc_method)
                case NpcMethod.Stouffer
                    x = NpcMethod.stouffer(u);
            end
        end
        function z = stouffer(u, epsilon)
            arguments
                u {mustBeNumeric,mustBeNonempty,mustBeVector}
                epsilon {mustBeNumeric,mustBeScalar,mustBePositive} = eps;
            end
            % Perform non-parametric combining to obtain a Stouffer's Z score.
            %
            % u: A vector of u-values to combine.
            % epsilon:
            %     Stouffer's method uses the inverse normal distribution on
            %     u-values. Valus of 0 or 1 will result in an infinite Z score.
            %     To avoid infinite scores, values of 0 will be changed to
            %     0 + epsilon and values of 1 will be changed to 1 - epsilon.
            %     Default: eps
            
            % Replace u=0 or u=1 with epsilon.
            % Necessary to prevent getting Inf or -Inf for z-values.
            u(u < epsilon) = epsilon;
            u(u > 1-epsilon) = 1 - epsilon;
            
            % Compute Stouffer's Z score.
            z = sum(norminv(1 - u)) / sqrt(length(u));
        end
    end
end