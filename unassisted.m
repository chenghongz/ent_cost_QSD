clear;

da = 3;
db = 3;

rho1 = RandomDensityMatrix(da*db,0,1);
sigma1 = RandomDensityMatrix(da*db,0,1);

cvx_begin sdp quiet
cvx_precision best
    variable E1( (da*db), (da*db)) hermitian
    variable E2( (da*db), (da*db)) hermitian

    prob = real(0.5*trace(E1*rho1) + 0.5*trace(E2*sigma1));
    maximise prob 
    subject to
        E1 >= 0;
        E2 >= 0;
        PartialTranspose(E1, 2, [da db]) >= 0;
        PartialTranspose(E2, 2, [da db]) >= 0;
        E1 + E2 == eye( (da*db) );

cvx_end
prob
