% Problem 3, PA03, MATH 4340
% Timothy Considine
classdef FVEMSystemMap < FixedPointMap

    properties 
        fvem;
    end

    methods

        function obj = FVEMSystemMap(fvem) %Cosntructor function
            obj.fvem=fvem;
        end

        function F = Evaluate(obj,r) %Evaluate function
            [AA, BB, QQ, QQprime] = obj.fvem.Get(r);
            F = r - (AA + QQprime)\(AA*r + BB + QQ);
        end

    end

end