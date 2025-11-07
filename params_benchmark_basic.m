function P = params_benchmark_basic()
% Minimal canonical Carvallo–Whipple matrices for q=[phi; delta].
% M*qdd + v*C1*qdot + (g*K0 + v^2*K2) q = B*tau_delta + E*w
% Use steer torque only; keep E empty at this stage.

g = 9.81; %#ok<NASGU> % used implicitly downstream, kept for clarity

% Canonical set (typical published benchmark-like values).
P.M  = [80.817,  2.319;
         2.319,  0.298];

P.C1 = [ 0.000, 33.866;
        -0.850,  1.685];

P.K0 = [-80.950,  -2.600;
         -2.600,  -0.803];

P.K2 = [ 0.000, 76.597;
         0.000,  2.654];

% Input: steer torque acts in steer equation only:
P.B  = [0; 1];

% no disturbances for the basic sanity checks:
%P.E  = []; % leave empty

% (We will add heading ψ and lateral deviation e later.)
% params_benchmark_basic.m  (add this line)
P.l_eq = 1.02;   % effective wheelbase / curvature gain (use your geometry)
% Add at the end of params_benchmark_basic.m
P.E = [1, 0;     % lateral/roll-equivalent push -> first column
       0, 1];    % steer-bias channel          -> second column same as B

end
