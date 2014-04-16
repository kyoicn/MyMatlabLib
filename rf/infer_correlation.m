function [c,sigma] = infer_correlation(mesh,C)
% INFER_CORRELATION infers the parameters of an exponential covariance
%
% [c,sigma] = infer_correlation(mesh,C)
%
% 'mesh' is a list of nodes in a mesh. 'C' is a given symmetric positive
% semi-definite covariance matrix. 
%
% This function uses a least squares procedure to infer the parameters 
% sigma and c_k of an exponential covariance function of the form
%
% cov(x1,x2) = sigma exp(-0.5 \sum_k (x1_k-x2_k)^2 / c_k^2 )
%
% Copyright 2010 Qiqi Wang (qiqi@mit.edu) and Paul G. Constantine 
% (pconsta@sandia.gov).

if nargin<2, error('Not enough input arguments.'); end

[n,d]=size(mesh);
A=zeros(n^2,d);

% build the least squares matrix
count=1;
for i=1:n
    x=mesh(i,:);
    for j=1:n
        y=mesh(j,:);
        for k=1:d
            A(count,k)=-0.5*(x(k)-y(k))^2;
        end
        count=count+1;
    end
end
A=[ones(n^2,1) A];

% remove covariance elements close to zero
C=C(:);
indz=find(C<sqrt(eps));
C(indz)=[];
A(indz,:)=[];

% construct the right hand side
b=log(C);

% solve the least squares problem
p=A\b;

% retrieve the parameters
sigma=exp(p(1));
c=1./sqrt(p(2:end));