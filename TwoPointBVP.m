% Problem 1, PA03, MATH 4340
% Timothy Considine
classdef TwoPointBVP
    properties 
        domain;
        kappa;
        f;
        f_exists;
        boundary_values;
        gamma_values;
        type_is_Dirichlet;
        second_order_component_is_linear;
        f_is_linear;
        D2f; 
    end
    methods
        function obj = TwoPointBVP(domain, kappa, bool) %Constructor
            obj.kappa = kappa;
            obj.domain = domain;
            obj.f_exists = false;
            obj.boundary_values = zeros(2,1);
            obj.gamma_values = zeros(2,1);
            obj.type_is_Dirichlet = [false; false];
            obj.second_order_component_is_linear = bool;
        end
        function obj = set_f(obj, f, d2f, bool) %f initializer
            obj.f = f;
            obj.D2f = d2f;
            obj.f_is_linear = bool;
            obj.f_exists = true;
        end
        function obj = set_left_boundary(obj, ga, la, bool) %left boundary initializer
            obj.boundary_values(1) = ga;
            obj.gamma_values(1) = la;
            obj.type_is_Dirichlet(1) = bool;
        end
        function obj = set_right_boundary(obj, gb, lb, bool) %right boundary initializer
            obj.boundary_values(2) = gb;
            obj.gamma_values(2) = lb;
            obj.type_is_Dirichlet(2) = bool;
            
        end
    end
end
