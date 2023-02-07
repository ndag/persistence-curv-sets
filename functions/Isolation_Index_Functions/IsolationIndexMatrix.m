function [alpha, splits] = IsolationIndexMatrix(dm)
    n = size(dm,1);
    
    % Edge case
    if n==2
        if dm(1,2)>0
            alpha = 2*dm(1,2);
            splits(1,:) = num2cell([1,2]);
        else
            alpha = 0;
            splits = cell(0,2);
        end
    % First non-trivial case: 3 elements
    else
        % Take the first 3 elements and compute their isolation indices
        delta = zeros(1,3);
        delta(1) = dm(2,1)+dm(1,3)-dm(2,3);
        delta(2) = dm(3,2)+dm(2,1)-dm(3,1);
        delta(3) = dm(1,3)+dm(3,2)-dm(1,2);

        % Find where the non-zero indices are (there must be at least one if dm
        % is not zero)
        A = find(delta);

        % The other element of the split is the complement of each row of A in
        % [1,2,3]
        alpha = zeros(1,numel(A));
        splits = cell(numel(A), 2);
        for i=1:numel(A)
            alpha(i) = delta(A(i));
            splits{i,1} = A(i);
            splits{i,2} = setdiff([1,2,3], A(i));
        end
    end
    
    % If we only have 3 elements, we're done. If not, we call the recursive
    % function
    if n>3
        m = 4;
        while m <= n
            [alpha, splits] = UpdateIsolationIndexMatrix(m, splits, alpha, dm);
            if ~isempty(alpha)
                m=m+1;
            else
                return
            end
        end
        % When this loop ends, we will have updated the isolation index for
        % all elements of the metric space.
    end
    % Lastly, we divide the isolation indices by 2 (we didn't do this
    % all throughout the calculation)
    alpha = alpha/2;
end