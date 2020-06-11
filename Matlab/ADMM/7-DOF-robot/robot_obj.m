classdef robot_obj
    %ROBOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n_links = 7;
        params  = zeros(1, 12*7);
    end
    
    methods
        function obj = robot_obj(params)
            % ROBOT Construct an instance of this class
            % Detailed explanation goes here
            obj.params = params;
           
        end
        function outputArg = mass_matrix(obj,q)
            %   METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = M_kuka(obj.params, q);
        end
        function outputArg = coriolis_matrix(obj,q,qd)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = C_kuka(obj.params, q,qd);
        end
        function outputArg = fwk(obj,q)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = FK_kuka(q);
        end
        function outputArg = jacobian(obj,q)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = Jac_kuka(q);
        end
        function outputArg = fkdyn(obj,x,u)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = fdyn(obj.params,x, u);
        end
        function outputArg = gravity(obj,q)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = G_kuka(obj.params, q);
        end
    end
end

