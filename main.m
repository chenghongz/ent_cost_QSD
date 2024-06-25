clear;
%% Required
d = 3; % dimension
bell2 = MaxEntangled(2)*MaxEntangled(2)';
rho = RandomDensityMatrix(d^2);
sigma = RandomDensityMatrix(d^2);


cvx_begin sdp quiet
cvx_solver sedumi
cvx_precision best
    variable W1(d^2,d^2) hermitian
    variable Q1(d^2,d^2) hermitian
    variable k nonnegative

    W2 = eye(d^2)-W1;
    Q2 = eye(d^2)-Q1;
    PT_W1 = PartialTranspose(W1, 2, [d d]);
    PT_W2 = PartialTranspose(W2, 2, [d d]);
    PT_Q1 = PartialTranspose(Q1, 2, [d d]);
    PT_Q2 = PartialTranspose(Q2, 2, [d d]);

    minimize k
    subject to
        0 <= W1 <= eye(d^2);
        0 <= Q1 <= eye(d^2);
        0 <= PT_Q1 <= eye(d^2);
        2*trace((rho-sigma) * W1) == TraceNorm(rho - sigma);
        SchattenNorm(PT_W1-PT_Q1, Inf) <= k;
        SchattenNorm(PT_W2-PT_Q2, Inf) <= k;

cvx_end
k

%% ent assisted suc probability

k1 = 2;
cvx_begin sdp quiet
cvx_solver sedumi
cvx_precision best

    variable W(d^2,d^2) hermitian
    variable Q(d^2,d^2) hermitian

    PT_W = PartialTranspose(W, 2, [d d]);
    PT_Q = PartialTranspose(Q, 2, [d d]);
    PT_W1 = PartialTranspose(eye(d^2)-W, 2, [d d]);
    PT_Q1 = PartialTranspose(eye(d^2)-Q, 2, [d d]);
    
    prob = trace(rho*W + sigma*(eye(d^2)-W))/2;

    maximize prob
    subject to
        0 <= W <= eye(d^2);
        0 <= Q <= eye(d^2);
        (1-k1)*PT_Q <= PT_W <= (1+k1)*PT_Q;
        (1-k1)*PT_Q1 <= PT_W1 <= (1+k1)*PT_Q1;
cvx_end

assisted = prob

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho1 = PermuteSystems(kron(rho,bell2), [1 3 2 4], [3 3 2 2]);
sigma1 = PermuteSystems(kron(sigma,bell2), [1 3 2 4], [3 3 2 2]);

cvx_begin sdp quiet
    variable W(d^2*4,d^2*4) hermitian

    PT_W = PartialTranspose(W, 2, [d*2 d*2]);
    prob1 = 0.5 + 0.5 * trace((rho1-sigma1)*W);
    maximize prob1
    subject to
        0 <= W <= eye(d^2*4);
        0 <= PT_W <= eye(d^2*4);
cvx_end

simplified = prob1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cvx_begin sdp quiet
cvx_solver SDPT3
    variable E(d^2,d^2) hermitian

    obj = 0.5*trace(E * rho) + 0.5*trace((eye(d^2)-E)*sigma);
    maximize real(obj)
    subject to
        eye(d^2) >= E >= 0;
        % eye(d^2) >= PartialTranspose(E, 2) >= 0; 
cvx_end
origin = obj

% cvx_begin sdp quiet
% cvx_solver sedumi
% cvx_precision best
% 
%     variable WAB(d^2,d^2) hermitian
%     variable QAB(d^2,d^2) hermitian
% 
%     PT_W = PartialTranspose(WAB, 2, [d d]);
%     PT_Q = PartialTranspose(QAB, 2, [d d]);
% 
% 
%     obj = 0.5*trace(WAB * rho) + 0.5 - 0.5*trace(WAB*sigma);
%     maximize obj
%     subject to
% 
%         eye(d^2) >= QAB >= 0;
%         eye(d^2) >= WAB >= 0;
% 
%         2* trace((rho-sigma) * WAB) == TraceNorm(rho - sigma);
% 
%         (1+k)*PT_Q - (k*eye(d^2)) <= PT_W <= k*eye(d^2) - (k-1)*PT_Q;
%         (1-k)*PT_Q <= PT_W <= (1+k) * PT_Q;
% 
% cvx_end
% 
% 
% obj
analytic = 0.5 + 0.25 * TraceNorm( rho - sigma)
