function Dn0 = D_n0_Ripser(dm)
N = size(dm,1);

Dn0 = cell(N,1);
Dn0{1} = [];

for n=2:N
    disp(n)
    bd_times = zeros(n-1,2,nchoosek(N,n));
    
    I = 1:n;
    idx = 1;
    while I
        PH0 = RipsFiltrationDM(dm(I,I), 0);
        bd_times(:,:,idx) = PH0{1}(1:n-1,:);
        
        I = nextcombi(N,I);
        idx = idx+1;
    end
    
    % Remove duplicates
    if n<N
        sizeBD = size(bd_times);
        bd_rows = reshape(bd_times,[],sizeBD(3))';
        bd_unique = unique(bd_rows,'stable','rows');
        bd_unique = reshape( bd_unique' ,sizeBD(1),sizeBD(2),[]);
    else
        bd_unique = bd_times;
    end
    
    % Store the results
    Dn0{n} = bd_unique;
end
end
