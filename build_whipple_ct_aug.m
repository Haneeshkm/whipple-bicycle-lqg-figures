function [Aa,Bua,Ca,Da] = build_whipple_ct_aug(v,P)
% x = [phi; delta; dphi; ddelta; psi; e], u = tau_delta
[A,Bu,Cy,Du] = build_whipple_ct(v,P);   % <-- 4×4 base (no recursion)
nx = size(A,1);

Aa = blkdiag(A, zeros(2));
Aa(nx+1,2)    = v/P.l_eq;   % psi_dot = (v/leq)*delta
Aa(nx+2,nx+1) = v;          % e_dot   = v*psi

Bua = [Bu; 0; 0];

% choose: with psi measurement (add a row)...
Ca = [Cy, zeros(size(Cy,1),2);
      zeros(1,nx), 1, 0];
Da = [Du; 0];
% ...or without psi measurement, use:
% Ca = [Cy, zeros(size(Cy,1),2)];
% Da = Du;
end
