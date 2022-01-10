classdef DecompositionModifiedCholeskySqrt < DecompositionModifiedCholesky

    methods
        function obj = DecompositionModifiedCholeskySqrt(A, dist, radius)
            if nargin == 1
                
                obj.SparcityPattern = A.SparcityPattern;
                obj.TMat = A.TMat;
                obj.Di = A.Di;
            else
                obj.recomputeSparcity(dist, radius, size(A, 1));

                obj.update(A);
            end
        end

        % MODIFY EVERYTHING HERE
        function X = mldivide(A, B)
            if isa(A, 'DecompositionModifiedCholeskySqrt')
                
                if isa(B, 'DecompositionModifiedCholesky')
                    %X = A.TMat.'*(A.Di.*(A.TMat*(B.TMat\(B.Di.\(B.TMat.')))));
                    error('NOT SUPPORTED')
                elseif isa(B, 'DecompositionModifiedCholeskySqrt')
                    error('NOT SUPPORTED');
                else
                    X = (sqrt(A.Di).*(A.TMat*B));
                end
            else
                %X = A\(B.TMat\(B.Di.\(B.TMat.')));
            end
        end

        function X = mrdivide(B, A)
            X = mldivide(A, B);
        end

        function X = mtimes(A, B)
            if isa(A, 'DecompositionModifiedCholeskySqrt')
                X = A.TMat\(sqrt(A.Di).\B);
            else
                error('NOT SUPPORTED');
            end
        end

        function s = size(obj)
            s = size(obj.SparcityPattern);
        end
    end
end
