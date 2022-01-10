classdef DecompositionModifiedCholesky < handle
    properties
        SparcityPattern
        TMat
        Di
    end

    methods
        function obj = DecompositionModifiedCholesky(A, dist, radius)
           if nargin == 0
               obj.SparcityPattern = [];
               obj.TMat = [];
               obj.Di = [];
           else
               obj.recomputeSparcity(dist, radius, size(A, 1));

               obj.update(A);
           end
        end

        function obj = update(obj, A)

            [n, N] = size(A);

            sparcity = obj.SparcityPattern;

            nel = nnz(sparcity) + n;

            Tis = zeros(1, nel);
            Tjs = zeros(1, nel);
            Tvs = zeros(1, nel);

            Tis(1:n) = 1:n;
            Tjs(1:n) = 1:n;
            Tvs(1:n) = ones(1, n);

            ind = n + 1;

            E = zeros(n, N);
            E(1, :) = A(1, :).';
            for i = 2:n
                js = find(sparcity(:, i));
                Zi = A(js, :);
                if size(Zi, 1) > size(Zi, 2)
                    % compute the Penrose inverse
                    bi = Zi*((Zi.'*Zi)\(A(i, :).'));
                else
                    bi = (Zi*Zi.')\(Zi*(A(i, :).'));
                end

                E(i, :) = A(i, :).' - Zi.'*bi;

                njs = numel(js);

                Tis(ind:(ind + njs - 1)) = i*ones(1, njs);
                Tjs(ind:(ind + njs - 1)) = js;
                Tvs(ind:(ind + njs - 1)) = -bi.';

                ind = ind + njs;
            end

            obj.TMat = sparse(Tis, Tjs, Tvs, n, n);

            d = sum(E.*E, 2)/(N - 1);

            obj.Di = 1./d;

        end

        function X = mldivide(A, B)
            if isa(A, 'DecompositionModifiedCholesky')
                
                if isa(B, 'DecompositionModifiedCholesky')
                    X = A.TMat.'*(A.Di.*(A.TMat*(B.TMat\(B.Di.\(B.TMat.')))));

                else
                    X = A.TMat.'*(A.Di.*(A.TMat*B));
                end
            else
                X = A\(B.TMat\(B.Di.\(B.TMat.')));
            end
        end

        function X = mrdivide(B, A)
            X = mldivide(A, B);
        end

        function X = mtimes(A, B)
            if isa(A, 'DecompositionModifiedCholesky')
                
                if isa(B, 'DecompositionModifiedCholesky')
                    X = A.TMat\(A.Di.\(A.TMat.'\(B.TMat\(B.Di.\(B.TMat.')))));
                else
                    X = A.TMat\(A.Di.\(A.TMat.'\B));
                end
            else
                X = (((A/B.TMat)./B.Di)/B.TMat.');
            end
        end

        function s = size(obj)
            s = size(obj.SparcityPattern);
        end

        function X = sqrtm(obj)
            X = DecompositionModifiedCholeskySqrt(obj);
        end


        function recomputeSparcity(obj, dist, radius, n)

            P = @(i, js) dist(i, js) <= radius;

            % I don't think that there is an easy way to pre-allocate this.
            % Fortunately, this is a function that should only be ran when
            % either the distances between the states change, or the number
            % of states changes.

            is = [];
            js = [];
            vs = [];

            for i = 2:n
                ne = find(P(i, 1:(i - 1)));
                is = [is, ne];
                js = [js, i*ones(1, numel(ne))];
                vs = [vs, ones(1, numel(ne))];
            end

            obj.SparcityPattern = sparse(is, js, vs, n, n);

        end



    end




end
