classdef (Abstract) FixedPointMap
    %%
    methods (Abstract)
        Evaluate(obj, u)
    end
    %%
    methods
        %%
        function [val, itercount, diff] = Iterate(obj, initguess, maxnumiter, TOL)
            prev = initguess;
            for n = 1:maxnumiter
                next = obj.Evaluate(prev);
                diff = norm(prev - next,inf);
                if diff<TOL
                    val = next;
                    itercount = n;
                    fprintf('Iterations converge in %d steps \n', itercount);
                    break;
                end
                prev = next;
            end
            if diff >= TOL
                fprintf('Iterations fail to converge in %d steps.\n', maxnumiter);
                val = next;
                itercount = maxnumiter;
            end
        end
        %%
    end
    %%
end
