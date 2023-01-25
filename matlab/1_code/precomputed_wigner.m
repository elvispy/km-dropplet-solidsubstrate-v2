function W = precomputed_wigner(harmonics_qtt)
    try
        data = load('Wigner.mat');
        W = data.W;
        order = size(W, 1);
        if size(W, 1) >= harmonics_qtt && ~isempty(data.W{end, end})
            W = W(1:harmonics_qtt, 1:harmonics_qtt);
            f = @(j1, j2, j3) wigner_shifted([j1 j2 j3], W, harmonics_qtt);
            return
        end
    catch 
        W = cell(0, 1);
        order = 0;
    end
    dummyW = cell(harmonics_qtt, harmonics_qtt);
    dummyW(1:order, 1:order) = W;
    W = dummyW;
    
    for kk = (order+1):harmonics_qtt
        for ll = 1:kk
            W{kk, ll} = zeros(1, 1+kk+ll - abs(kk-ll));
            for mm = abs(kk-ll):(kk+ll)%(abs(kk-ll) + mod(abs(kk-ll), 2)):(kk+ll)
               W{kk, ll}(mm+1) = Wigner3j([kk ll mm], [0 0 0]);
            end
        end
    end    
    save('Wigner.mat', 'W');
    f = @(J) wigner_shifted(J, W, harmonics_qtt);
end


function r = wigner_shifted(J, W, maxx)
    if J(1) < J(2); r = wigner_shifted([J(2), J(1), J(3)], W); return; end
    if J(1) > maxx || J(2) > maxx || J(3) > J(1) + J(2) 
        error("Maximum allowed index is %d but got %d", maxx, max(J)); 
    end
    r = W{J(1), J(2)}(J(3) + 1);
end