function [A,Bu,Cy,Du] = build_whipple_ct(v, P)
% 4-state Whipple: x = [phi; delta; dphi; ddelta], u = tau_delta
g = 9.81;

M  = P.M;   C1 = P.C1;   K0 = P.K0;   K2 = P.K2;  Bq = P.B;
Minv = inv(M);

A11 = zeros(2);      A12 = eye(2);
A21 = -Minv*(g*K0 + v^2*K2);
A22 = -Minv*(v*C1);
A   = [A11 A12; A21 A22];

Bu  = [zeros(2,1); Minv*Bq];

% nominal sensors: [phi; dphi; delta]
Cy  = [1 0 0 0;
       0 0 1 0;
       0 1 0 0];
Du  = zeros(size(Cy,1),1);
end
