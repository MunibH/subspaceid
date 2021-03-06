function out = my_smooth(x, N)
% operates on first dimension only
% kern = gausswin(N)

    if N==1
        out = x;
        return;
    end

    Ncol = size(x, 2);
    Nel = size(x, 1);
    if mod(N, 2)==0
        N = N+1;
    end
    
    kern = gausswin(N);
    if N>1
        kern(1:floor(N/2)) = 0; %causal
    end
    kern = kern./sum(kern);
    
    out = zeros(Nel, Ncol);
    for j = 1:Ncol
        out(:, j) = conv(x(:, j), kern, 'same');
    end

end