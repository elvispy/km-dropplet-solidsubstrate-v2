function WS = my_wigner3j(J123)

% Compute the Wigner 3j symbol using the Racah formula. 
%
% W = Wigner3j( J123, M123 ) 
%
% J123 = [J1, J2, J3].
% M123 = [ 0,  0,  0].
% All Ji's and Mi's have to be integeres or half integers (correspondingly).
%
% According to seletion rules, W = 0 unless:
%   |Ji - Jj| <= Jk <= (Ji + Jj)    (i,j,k are permutations of 1,2,3)
%   0 <= Ji    (i = 1,2,3)
%    M1 + M2 + M3 = 0

j1 = J123(1, :); j2 = J123(2, :); j3 = J123(3, :);


end