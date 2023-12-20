%=========================================================================
%                        BAYESIAN INFERENCE 
%                   DRAWS FROM POSTERIOR DISTRIBUTION
%=========================================================================



%=========================================================================
%     DEFINITION OF DATA, LAG STRUCTURE AND POSTERIOR SIMULATION
%=========================================================================

[Tdummy,n] = size(YYdum);
[Tobs,n]   = size(YYact);
X          = [XXact; XXdum];
Y          = [YYact; YYdum];
p          = spec(1);                 % Number of lags in the VAR
T          = Tobs+Tdummy;             % Number of observations
F          = zeros(n*p,n*p);          % Matrix for Companion Form
I          = eye(n);

for i=1:p-1
    F(i*n+1:(i+1)*n,(i-1)*n+1:i*n) = I;
end


%=========================================================================
%               OLS ESTIMATOR FOR PHI AND SSR (SIGMA)
%=========================================================================
inv_x     = inv(X'*X);
Phi_tilde = (inv_x)*X'*Y;
Sigma     = (Y-X*Phi_tilde)'*(Y-X*Phi_tilde);

% Matrices for collecting draws from Posterior Density

Sigmap    = zeros(nsim,n,n);
Phip      = zeros(nsim,n*p+1,n);
counter   = 0;

%=========================================================================
%            DRAWS FROM POSTERIOR DENSITY (DIRECT SAMPLING)
%=========================================================================


for j=1:nsim

    
    % Draws from the density Sigma | Y
    
    sigma   = iwishrnd(Sigma,T-n*p-1);
    
    % Draws from the density vec(Phi) |Sigma(j), Y
    
    phi_new = mvnrnd(reshape(Phi_tilde,n*(n*p+1),1),kron(sigma,inv_x));
    
    % Rearrange vec(Phi) into Phi
    
    Phi     = reshape(phi_new,n*p+1,n);
        
    Sigmap(j,:,:) = sigma;
    Phip(j,:,:)   = Phi;
    
    Phi = Phi(1:n*p,:);

    counter         = counter +1; 
     
     if counter==100    
    
    disp('                                                              ');
    disp(['       DRAW NUMBER:   ', num2str(j)]                          );
    disp('                                                              ');
    disp(['       REMAINING DRAWS:   ', num2str(nsim-j)]                 );
    disp('                                                              ');

     counter = 0;
     
     end
     
end

%=========================================================================
%            POSTERIOR MEDIAN : FORECAST
%=========================================================================

post_phi    = squeeze(median(Phip(nburn:end,:,:)));
post_sig    = squeeze(median(Sigmap(nburn:end,:,:)));
