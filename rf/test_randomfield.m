function r = test_randomfield()

% This test script tests the options on the randomfield function.

% Build a 1D grid.
nx=50; 
mesh=linspace(-1,1,nx)'; % make sure it's a column vector.

% Condition on known endpoints.
data.x=[-1; 0; 1]; data.fx=[0; 1; 1];

% Set up the two point covariance function.
c=1; sigma=1;
cv=@(x1,x2) gp_exp_cov(x1,x2,c,sigma);
cv_sparse=@(x1,x2) gp_exp_cov(x1,x2,0.1,sigma);

% Generate external covariance matrices
C=covariance_matrix(cv,mesh,[],[]);
A=covariance_matrix(cv,data.x,[],[]);
B=covariance_matrix(cv,mesh,data.x,[]);

% Test options.

fprintf('Test standard\n');
tic; F=randomfield(cv,mesh); t=toc;
fprintf('Success!\t\t time=%f\n',t);
clear F


fprintf('Test given cvfun\n');
tic; F=randomfield(C,mesh); t=toc;
fprintf('Success!\t\t time=%f\n',t);
clear F


fprintf('Test given cvmat\n');
cvmat.C=C;
tic; F=randomfield(cvmat,mesh); t=toc;
fprintf('Success!\t\t time=%f\n',t);
clear F


fprintf('Test nsamples\n');
tic; F=randomfield(cv,mesh,'nsamples',10);  t=toc;
fprintf('Success!\t\t time=%f\n',t);
X=F;
clear F


fprintf('Test given data and cvfun\n');
tic; F=randomfield(cv,mesh,'data',data); t=toc;
fprintf('Success!\t\t time=%f\n',t);
clear F


fprintf('Test given data and cvmat\n');
tic; F=randomfield(C,mesh,'data',data); t=toc;
fprintf('Success!\t\t time=%f\n',t);
clear F


fprintf('Test filter\n');
tic; F=randomfield(cv,mesh,'filter',0.5); t=toc;
fprintf('Success!\t\t time=%f\n',t);
clear F


fprintf('Test truncation\n');
tic; F=randomfield(cv,mesh,'trunc',10); t=toc;
fprintf('Success!\t\t time=%f\n',t);
clear F


fprintf('Test sparse threshold\n');
tic; F=randomfield(cv_sparse,mesh,'spthresh',1e-5); t=toc;
fprintf('Success!\t\t time=%f\n',t);
clear F


fprintf('Test given mean\n');
tic; F=randomfield(cv_sparse,mesh,'mean',ones(size(mesh))); t=toc;
fprintf('Success!\t\t time=%f\n',t);
clear F


fprintf('Test given cvfun in options\n');
tic; F=randomfield(C,mesh,'cvfun',cv_sparse); t=toc;
fprintf('Success!\t\t time=%f\n',t);
clear F


fprintf('Test given cvmat in options\n');
cvmat.A=A; cvmat.B=B;
tic; F=randomfield(cv,mesh,'cvmat',cvmat); t=toc;
fprintf('Success!\t\t time=%f\n',t);
clear F


fprintf('Test given snapshots and cvfun\n');
tic; F=randomfield(cv,mesh,'snaps',X); t=toc;
fprintf('Success!\t\t time=%f\n',t);
clear F


fprintf('Test given snapshots and cvmat\n');
tic; F=randomfield(C,mesh,'snaps',X); t=toc;
fprintf('Success!\t\t time=%f\n',t);
clear F


fprintf('Test given snapshots and cvmat and data\n');
tic; F=randomfield(C,mesh,'snaps',X); t=toc;
fprintf('Success!\t\t time=%f\n',t);
clear F


% 2d test
x=linspace(-1,1,10); 
[U,V]=meshgrid(x,x); mesh=[U(:) V(:)];
data.x=[-1 -1; 0 0; 1 1]; data.fx=[0; 1; 1];

fprintf('Test 2d with data\n');
tic; F=randomfield(cv,mesh,'data',data); t=toc;
fprintf('Success!\t\t time=%f\n',t);
clear F

r=0;
