classdef MetricGraph
    % This class defines a metric graph. It requires the weighted adjacency
    % matrix A and a flag to say if it's directed or not. The main purpose
    % of using a class is to have all methods (calculate distances, build
    % distance matrices, take samples with different weights) in one place.
    properties
        A;                  % Weighted adjacency matrix. Need not be symmetric.
        is_directed;        % Boolean flag that indicates if the graph is directed
        sampling_method     % String that describes the probability distribution W
        W;                  % Matrix with edge probabilities.
        % W has to satisfy size(A)==size(W) and either sum(W, 'all')==1 if
        % G is directed or isequal(W,triu(W)) and sum(W, 'all')==1 if not.
    end
    
    methods
        % Constructor. Just supply the weighted adjacency matrix and
        % indicate if the graph is directed or not.
        % W defaults to a uniform distribution, i.e.
        % P(segment)=|segment|/sum(length of all edges)
        % If the user asked for a non-directed graph but supplies a
        % non-symmetric adjacency matrix A, we will make it symmetric by
        % copying the upper triangle into the lower triangle.
        function obj = MetricGraph(A, sampling_method, is_directed)
            if nargin<2
                sampling_method='uniform';
            end
            if nargin<3
                is_directed=false;
            end
            
            % If we are building a non-directed graph, we need to check if
            % A is symmetric. If not, we write the upper triangle into the
            % lower triangle to ensure symmetry.
            if ~is_directed && ~isequal(A,A')
                n = size(A,1);
                for i=1:n
                    for j=1:(i-1)
                        A(i,j) = A(j,i);
                    end
                end
            end
            
            % Save the values
            obj.is_directed = is_directed;
            obj.A = A;
            obj.sampling_method = sampling_method;
            obj.W = get_sampling_matrix(obj, A, sampling_method, is_directed);
            
            % We need to load Dijkstra's algorithm
            addpath('../../functions/dijkstra');
        end
        
        % Calculate the distance between the points P1 and P2. Each point
        % is a 3D vector [a, b, t] where a, b are integers that represent
        % vertices of the graph. Since a metric graph is a continuous
        % object, we represent a point in an edge with a number
        % 0 <= t <= 1.
        % To find the shortest path, we use Dijstra's algorithm to find the
        % distance between (say) a1 and a2, then add the distance of P1 to
        % a1 and P2 to a2 (using t1 and t2). Repeat that for all four
        % choices d(a1/b1, a2/b2) and pick the shortest one.
        function d = distance(obj, P1, P2)
            P1 = num2cell(P1);
            P2 = num2cell(P2);
            [a1, b1, t1] = P1{:};
            [a2, b2, t2] = P2{:};
            
            % If the two points are contained in the same edge, the
            % distance is simply abs(t2-t1) times the edge length.
            if (~obj.is_directed && max(a1,b1)==max(a2,b2) && min(a1,b1)==min(a2,b2))...
                    || (obj.is_directed && isequal([a1,b1],[a2,b2]))
                % We check two cases:
                % If the graph is non-directed, we only need the boundaries
                % of the edges to be the same. The order doesn't matter, so
                % we only compare the max and min of each edge.
                % If the graph is directed, the order matters so we compare
                % the edges as vectors.
                
                d = obj.A(a1,b1)*abs(t2-t1);
            % If the points are in different edges, we have to do more
            else
                e1 = obj.A(a1,b1);
                e2 = obj.A(a2,b2);

                % If the graph is directed, the distance between P1 and P2
                % is the sum of the distance from P1 to b1, then b1 to a2
                % (using Dijkstra) and then a2 to P2.
                if obj.is_directed
                    d = (1-t1)*e1 + dijkstra(obj.A,b1,a2) + t2*e2;
                % If the graph is non-directed, points inside of an edge
                % can choose to go to either boundary point. Thus, the
                % distance between P1 and P2 is
                %       d(P1,Bnd1)+d(Bnd1,Bnd2)+d(Bnd2,P2),
                % where each Bndi can be ai or bi. We choose Bnd1 and Bnd2 
                % to minimize the distance.
                else
                    % Find the distances between the boundaries of the edges that
                    % contain our points
                    d11 = dijkstra(obj.A, a1, a2);
                    d12 = dijkstra(obj.A, a1, b2);
                    d21 = dijkstra(obj.A, b1, a2);
                    d22 = dijkstra(obj.A, b1, b2);

                    % Now find the distances between our points by adding the
                    % distances between the boundaries of the edges to the distance
                    % of our points to the corresponding boundary
                    d11 = t1*e1 + d11 + t2*e2;
                    d12 = t1*e1 + d12 + (1-t2)*e2;
                    d21 = (1-t1)*e1 + d21 + t2*e2;
                    d22 = (1-t1)*e1 + d22 + (1-t2)*e2;

                    % The metric is the shortest path distance
                    d = min([d11, d12, d21, d22]);
                end
            end
        end
        
        % Get the distance matrix for a set of points. "points" should be
        % an n-by-3 matrix, where each row is a point in the format
        % detailed in the function distance (i.e. [a,b,t]).
        function dm = distance_matrix(obj, points)
            if nargin==1
                n = size(obj.A,1);
                points = [(1:n)', (2:n+1)', zeros(n,1)];
            end
            
            % Number of points
            n = size(points,1);
            
            dm = zeros(n);
            
            % Find all pairwise distances
            for i=1:n
                for j=(i+1):n
                    d = distance(obj, points(i,:), points(j,:));
                    
                    dm(i,j) = d;
                    dm(j,i) = d;
                end
            end
        end
        
        % Builds a probability matrix W according to the chosen
        % sampling_method.
        function W = get_sampling_matrix(obj, A, sampling_method, is_directed)
            % Every point gets the same probability. This is calculated by
            % weighting each edge by its length divided by the sum of the
            % length of all edges.
            if strcmp(sampling_method,'uniform')
                if is_directed
                    W = A/sum(A,'all');
                else
                    W = triu(A)/sum(triu(A),'all');
                end
            % Every edge gets the same probability, regardless of length.
            elseif strcmp(sampling_method,'uniform_edges')
                if is_directed
                    W = (A~=0)/sum(A(:)~=0);
                else
                    W = triu(A~=0)/sum(triu(A~=0),'all');
                end
            end
        end
        
        % Get a sample of n points.
        % We choose an edge e=[a,b] and a point in that edge by choosing a
        % t=Unif([0,1]). 
        % User can supply a new edge probability matrix W.
        % It has to satisfy size(W)=size(A) and either sum(W,'all')==1 if G
        % is directed or isequal(W,triu(W)) and sum(triu(W),'all')==1 if not.
        function [X, dm] = sample(obj, n, W)
            if nargin < 3
                W = obj.W;
            end
            
            % Get the number of vertices of the graph
            N = size(obj.A,1);
            
            % Select n random edges
            % Note: we need to reshape M into a row vector.
            E = randsrc(n,1, [1:N*N; reshape(W,1,[]) ] );
            
            % Notice that E contains numbers in 1:nEdges, but we want
            % edges [a,b]. MATLAB turns matrices into column vectors by
            % appending the columns of a matrix, which means that the index
            % 1 is the edge [1,1], 2 is [2,1] and, if A is N-by-N, the
            % index n+1 is the edge [1,2].
            % To fix that, we write an index e=q*n+r, where 1<=r<=n, and
            % turn the index e into the edge [r,q].
            A0 = N+mod(E,-N);
            B0 = (E-A0)/N+1;
            
            % Now we select a t=Unif([0,1]) to represent a point inside of
            % each edge
            T = rand(size(E));
            
            % We put everything together into X
            if obj.is_directed
                % If the graph is directed, we want to keep the direction
                % of the edges
                X = [A0,B0,T];
            else
                % If the graph is non-directed, we make the edges [a,b]
                % satisfy a<b for convenience
                X = [min(A0,B0), max(A0,B0), T];
            end
            
            % And calculate the distance matrix
            dm = distance_matrix(obj, X);
        end
        
        % Calculates the first Betti number of the graph
        % This uses the Euler characteristic:
        % B0-B1=#vertices - #edges
        % This only works for non-directed graphs. I don't know how to do
        % this for directed graphs yet.
        function B1 = Betti_1(obj)
            nVertices = size(obj.A, 1);
            nEdges = sum(triu(obj.A) ~= 0, 'all');
            
            B0 = Betti_0(obj);
            
            B1 = B0 - nVertices + nEdges;
        end
        
        % Dimension of the 0th homology of the graph, i.e. the number of
        % connected components
        function B0 = Betti_0(obj)
            % Convert my graph to a Matlab graph object without weights
            % (Could I have started this instead of having to program it myself?)
            A2 = (obj.A) ~= 0;
            if obj.is_directed
                G2 = digraph(A2);
            else
                G2 = graph(A2);
            end
            
            % For each vertex, conncomp gives the index of the connected
            % comopnent that contains the vertex. We take the max to find
            % how many components there are
            B0 = max(conncomp(G2));
        end
    end
end