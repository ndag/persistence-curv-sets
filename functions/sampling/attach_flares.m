function A = attach_flares(dm_0,I,Flr_Size)
    n = size(dm_0,1);
    % We only retain distances between consecutive points
    for i=1:n
        temp_row = dm_0(i,:);

        dm_0(i,:)=0;

        if i==1
            dm_0(i,n)=temp_row(n);
            dm_0(i,2)=temp_row(2);
        elseif i==n
            dm_0(i,n-1)=temp_row(n-1);
            dm_0(i,  1)=temp_row(1);
        else
            dm_0(i,i-1)=temp_row(i-1);
            dm_0(i,i+1)=temp_row(i+1);
        end
    end

    % We will attach flares to the points listed in the index list I by
    % attaching an identity matrix and a zero matrix
    O = zeros(n);
    Id = eye(n)*Flr_Size;

    A = [dm_0,   Id(:,I);
        Id(I,:), O(I,I)];
end