function Dn0 = D_n0_um(dm, progress)
if nargin<2
    progress=false;
end

if progress
    t_total = 0;
end

N = size(dm,1);

% Initialize matrices to store diagrams
Dn0 = cell(N,1);
for n=1:N
    Dn0{n} = zeros(nchoosek(N,n),n-1);
end

for n=2:N
    if progress
        fprintf('%i/%i: ', n,N)
        tic
    end
    
    J=1:(n-1);
    jdx = 1;
    idx = 1;
    while J
        % disp('----------')
        % disp(J)
        % disp('----------')
        D_old = Dn0{n-1}(jdx,:);
        
        p_new = max(J)+1;
        while p_new<=N
            % disp(p_new)
            dist_new = dm(J,p_new);
            % disp([D_old, min(dist_new)]);
            Dn0{n}(idx,:) = [D_old, min(dist_new)];

            p_new = p_new+1;
            idx = idx+1;
        end

        J = nextcombi(N,J);
        jdx = jdx+1;
        while max(J)==N
            J = nextcombi(N,J);
            jdx = jdx+1;
        end
    end
    
    if progress
        dt = toc;
        fprintf('%0.2f\n', dt);
        t_total = t_total+dt;
    end
end

% Sort rows and remove duplicates
Dn0 = cellfun(@(M) unique(sort(M,2),'rows'), Dn0, 'UniformOutput', false);

if progress
    fprintf('Total time: %0.2f\n', t_total);
end
end