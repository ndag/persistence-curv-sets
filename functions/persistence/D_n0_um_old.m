function Dn0 = D_n0_um_old(dm, progress, prog_percent)
if nargin<2
    progress=false;
end
if nargin<3
    prog_percent = 0.10;
end

N = size(dm,1);
I0 = 1:N;

if progress
    step = round((2^N-1)*prog_percent);
    t_total = 0;
    tic;
end

% Initialize matrices to store diagrams
Dn0 = cell(N,1);
for n=1:N
    Dn0{n} = zeros(nchoosek(N,n),n-1);
end

% Count how many diagrams of each size we've found
count = ones(N,1);
for ii=1:(2^N-1)
    % Get a subset of 1:N
    idx = logical( dec2bin(ii,N)' - '0' );
    
    % We're only interested in subsets of size>=2
    n = sum(idx);
    if n>=2
        I = I0(idx);
        Dn0{n}(count(n),:) = dgm_0_um(dm(I,I));
        count(n) = count(n)+1;
    end
    
    if progress
        if mod(ii,step)==0
            dt = toc;
            t_total = t_total+dt;
            fprintf('%0.2f: %0.2f\n', ii/(2^N-1), dt);
            tic;
        end
    end
end

% Remove duplicates
Dn0 = cellfun(@(M) unique(M,'rows'), Dn0, 'UniformOutput', false);

if progress
    fprintf('Total time: %0.2f\n', t_total);
end
end