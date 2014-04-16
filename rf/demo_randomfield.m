% RANDOM FIELD DEMO
close all; clear all;
disp('This demo runs three examples of generating realizations of a Gaussian random field.');
disp('Press any key to continue.');
pause

%--------------------------------------------------------------------------
% A simple one-dimensional Gaussian random field with an exponential
% covariance conditioned on known endpoints.
%--------------------------------------------------------------------------

% Build a 1D grid.
nx=250; 
mesh=linspace(-1,1,nx)'; % make sure it's a column vector.

% Condition on known endpoints.
data.x=[-1; 1]; data.fx=[0; 1];

% Set up the two point covariance function.
c=0.1; sigma=1;
cv=@(x1,x2) gp_exp_cov(x1,x2,c,sigma);

% Make 10 unfiltered samples.
F=randomfield(cv,mesh,'nsamples',10,'data',data);

% Plot the realizations.
figure;
plot(mesh,F);
title('Random Field Realizations');

disp('Press any key to continue.');
pause
close(gcf);

%--------------------------------------------------------------------------
% A filtered one-dimensional Gaussian process with the same conditioning.
%--------------------------------------------------------------------------

% Get filtered realizations.
F_filt=randomfield(cv,mesh,'nsamples',10,'data',data,'filter',0.6);

% Compare filtered and unfiltered.
figure;
subplot(2,1,1), plot(mesh,F); title('Unfiltered');
subplot(2,1,2), plot(mesh,F_filt); title('Filtered');

disp('Press any key to continue.');
pause
close(gcf);

%--------------------------------------------------------------------------
% A two-dimensional Gaussian process. 
%--------------------------------------------------------------------------

% Build a 2D grid.
nx=20; x=linspace(-1,1,nx)';
[X,Y]=meshgrid(x,x);
mesh=[X(:) Y(:)];

% Set up the covariance function with anisotropic correlation.
c=[0.1 1]; sigma=1;
cv=@(x1,x2) gp_exp_cov(x1,x2,c,sigma);

% Generate 30 samples with an explicit truncation and sparsity.
F=randomfield(cv,mesh,'nsamples',30,'trunc',100,'spthresh',1e-6);

% Plot four of the realizations
figure;
subplot(2,2,1), surf(X,Y,reshape(F(:,1),nx,nx));
subplot(2,2,2), surf(X,Y,reshape(F(:,2),nx,nx));
subplot(2,2,3), surf(X,Y,reshape(F(:,3),nx,nx));
subplot(2,2,4), surf(X,Y,reshape(F(:,4),nx,nx));

disp('Press any key to continue.');
pause
close(gcf);

%--------------------------------------------------------------------------
% User input snapshots and conditioned on data in 2d.
%--------------------------------------------------------------------------

% Use the previous fields as synthetic data.
D=F;

% Condition on "known" data points.
data.x=[-1 -1; 0 0; 1 1]; data.fx=[0; 1; -1];

% Set up a dummy covariance matrix. The code ignores this to infer a
% covariance function from the given data; it uses this to construct the
% data-data and data-unknowns covariance matrices.
C=eye(nx^2);

% Generate 4 samples with an explicit truncation and sparsity.
F=randomfield(C,mesh,'nsamples',4,'data',data,'snaps',D);

% Plot four of the realizations
figure;
subplot(2,2,1), surf(X,Y,reshape(F(:,1),nx,nx));
subplot(2,2,2), surf(X,Y,reshape(F(:,2),nx,nx));
subplot(2,2,3), surf(X,Y,reshape(F(:,3),nx,nx));
subplot(2,2,4), surf(X,Y,reshape(F(:,4),nx,nx));

disp('Press any key to end.');
pause
close(gcf);








