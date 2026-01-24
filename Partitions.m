% Problem 1, PA02, MATH 4340
% Jack Nyman
% Tim Considine
classdef Partitions
    %
    properties (SetAccess=public)
        numintervals;
        x;
        z;
        h;
        dz;
    end
    %
    methods
        function obj = Partitions(x)
            obj.numintervals = size(x,1) - 1;
            obj.x = x;
            obj.z = [obj.x(1); 
                obj.x(1:obj.numintervals) + ((obj.x(2:obj.numintervals+1))-obj.x(1:obj.numintervals))./2;
                obj.x(end)];
            obj.h = obj.x(2:obj.numintervals+1) - obj.x(1:obj.numintervals);
            obj.dz = obj.z(2:obj.numintervals+2) - obj.z(1:obj.numintervals+1);
        end
    end
end
