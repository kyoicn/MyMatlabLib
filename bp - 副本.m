% BP.M
% 
%  Basis pursuit for nonequispaced Fourier matrix
%  A = exp(-2*pi*i*x*((-N/2):(N/2-1))) as measurement matrix.
% 
%  f_hat=bp(xx,f,N,mode,re,fth,res_tol)
% 
%  xx      indices of the nodes x_j=X(j)/N, i.e., X \subset {-N/2,...,N/2-1} or
%          nodes x \in [-1/2,1/2]
%  f       measurements, i.e., rhs of A f_hat = f
%  N       polynomial degree
% 
%  mode    'cvx'
%          'mosek'
%          'l1magic'      only for 'real'
%          'l1magic_ft'   only for 'real'
%          'linprog'      only for 'real'
% 
%  re      'real'         search for real solution vector
%          'complex'      search for complex solution vector
% 
%  fth     @rdft          coefficient restricted, nodal restricted discrete FT
%          @rfft          coefficient restricted, nodal restricted fast FT
%          @rndft         coefficient restricted, nonequispaced discrete FT
%          @rnfft         coefficient restricted, nonequispaced fast FT
% 
%  res_tol                tolerance for residual, used only in l1magic, cvx
%                         and linprog
%   
%  max_iter               only added for compatibility with routine omp
%                         does not have any effect and can be omitted in call
% 
%  Author  Holger Rauhut, Stefan Kunis
% 
% ---------------------------------------------------------------
%
% COPYRIGHT : (c) NUHAG, Dept.Math., University of Vienna, AUSTRIA
%             http://nuhag.eu/
%             Permission is granted to modify and re-distribute this
%             code in any manner as long as this notice is preserved.
%             All standard disclaimers apply.
%
function [f_hat,infos]=bp_copy(matrix, xx,f,N,mode,re,fth,res_tol,max_iter)
%
% Basis pursuit for nonequispaced Fourier matrix
% A = exp(-2*pi*i*x*((-N/2):(N/2-1))) as measurement matrix.
%
% f_hat=bp(xx,f,N,mode,re,fth,res_tol)
%
% xx      indices of the nodes x_j=X(j)/N, i.e., X \subset {-N/2,...,N/2-1} or
%         nodes x \in [-1/2,1/2]
% f       measurements, i.e., rhs of A f_hat = f
% N       polynomial degree
%
% mode    'cvx'
%         'mosek'
%         'l1magic'      only for 'real'
%         'l1magic_ft'   only for 'real'
%         'linprog'      only for 'real'
%
% re      'real'         search for real solution vector
%         'complex'      search for complex solution vector
%
% fth     @rdft          coefficient restricted, nodal restricted discrete FT
%         @rfft          coefficient restricted, nodal restricted fast FT
%         @rndft         coefficient restricted, nonequispaced discrete FT
%         @rnfft         coefficient restricted, nonequispaced fast FT
%
% res_tol                tolerance for residual, used only in l1magic, cvx
%                        and linprog
%  
% max_iter               only added for compatibility with routine omp
%                        does not have any effect and can be omitted in call
%
% Author  Holger Rauhut, Stefan Kunis

infos='';
if(strcmp(func2str(fth),'rdft')||strcmp(func2str(fth),'rfft'))
  x=xx/N;
else
  x=xx;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Guessing parameter for inner iterations and nfft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lsqr_tol=res_tol/sqrt(length(xx));
flag=0;

ft_opt.method='gaussian';
ft_opt.sigma=2;
ft_opt.m=ceil(-log10(lsqr_tol));

switch(mode)
  case 'cvx'
    t=cputime;
    cvx_begin
    cvx_quiet(true);
    cvx_precision(res_tol);

    if(strcmp(re,'real'))
      variable f_hat(N)
    else
      variable f_hat(N) complex
    end;

    minimize( norm(f_hat,1) )
    subject to
      exp(-2*pi*i*x*((-N/2):(N/2-1)))*f_hat == f;
    cvx_end
    t=cputime-t;
    infos=strcat(infos,sprintf(' Time=%1.1e',t));
  case 'mosek'
    y=[real(f);imag(f)];    
    A = exp(-2*pi*i*x*((-N/2):(N/2-1)));
    B=real(A);
    C=imag(A);
 
    % setup
    %clear prob param;
    clear prob;
    if(strcmp(re,'real'))
      prob.c=ones(2*N,1);
      prob.a=sparse([eye(2*N); [B,-B]; [C,-C]]);
      prob.blc = [zeros(2*N,1);y];
      prob.buc = [inf*ones(2*N,1);y];
      prob.blx = -inf*ones(2*N,1);
      prob.bux = inf*ones(2*N,1);

      %specify accuracy
      %param.MSK_DPAR_INTPNT_TOL_REL_GAP = res_tol;
      
      % solve the optimization problem      
      %[r,res] = mosekopt('minimize echo(0) info',prob, param);
      [r,res] = mosekopt('minimize echo(0) info',prob);
      
      if(~isfloat(res)                              && ...
	 sum(strcmp(fieldnames(res),        'sol')) && ...
         sum(strcmp(fieldnames(res.sol),    'itr')) && ...
         sum(strcmp(fieldnames(res.sol.itr),'xx')))

        xx=res.sol.itr.xx;
        infos=strcat(infos,sprintf(' Time=%1.1e',...
	    		   res.info.MSK_DINF_OPTIMIZER_CPUTIME));
        f_hat = xx(1:N)-xx((N+1):(2*N));
      else
        infos=strcat(infos,'Got no solution from mosek!');
        f_hat=zeros(N,1);
      end;
    else
      M=length(x);
      Z=zeros(M,N);
      prob.c   = [ones(N,1);zeros(2*N,1)];
      prob.a   = sparse([[Z, B, -C]; [Z, C, B]]);
      prob.blc = y;
      prob.buc = y;
      prob.blx = -inf*ones(3*N,1);
      prob.bux = inf*ones(3*N,1);
      prob.cones = cell(N,1);
      % specify all the cones  
      for k = 1:N
        prob.cones{k}.type = 'MSK_CT_QUAD';
        prob.cones{k}.sub  = [k; k+N; k+2*N];
      end
      
      %specify accuracy ?
      %param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = res_tol;
      %param.MSK_DPAR_INTPNT_TOL_REL_GAP = res_tol;
      
      % solve the optimization problem
      %[r,res] = mosekopt('minimize echo(0) info',prob,param);
      [r,res] = mosekopt('minimize echo(0) info',prob);
      
      if(~isfloat(res)                              && ...
	 sum(strcmp(fieldnames(res),        'sol')) && ...
         sum(strcmp(fieldnames(res.sol),    'itr')) && ...
         sum(strcmp(fieldnames(res.sol.itr),'xx')))

        xx=res.sol.itr.xx;
        infos=strcat(infos,sprintf(' Time=%1.1e',...
	      		   res.info.MSK_DINF_OPTIMIZER_CPUTIME));

        f_hat = xx((N+1):(2*N)) + i* xx((2*N+1):(3*N));
      else
        infos=strcat(infos,'Got no solution from mosek!');
        f_hat=zeros(N,1);
      end;
    end;
  case 'l1magic'
    y=[real(f);imag(f)];

    A = exp(-2*pi*i*x*((-N/2):(N/2-1)));
    B=real(A);
    C=imag(A);

    t=cputime;
    [f_hat,lsqrflag]=lsqr([B;C],y,lsqr_tol,length(x));
    [f_hat,l1infos]=l1eq_pd(f_hat, [B;C] , [], y);
    t=cputime-t;
    %infos=strcat(sprintf('\n\n'),l1infos);
    infos=strcat(infos,sprintf(' Time=%1.1e',t));
  case 'l1magic_ft'
    y=[real(f);imag(f)];

    t=cputime;
    % Compute the minimal 2-norm solution first (matrix free)
    [f_hat,lsqrflag]=lsqr(@(a,tflag)ft(a,1:N,xx,N,re,fth,ft_opt,tflag),y,...
			  lsqr_tol,length(x));

    % Compute minimal 1-norm solution
    [f_hat,l1infos]=l1eq_pd(f_hat, ...
		     @(a)ft(a,1:N,xx,N,re,fth,ft_opt,'notransp'),...
		     @(a)ft(a,1:N,xx,N,re,fth,ft_opt,'transp'),...
		     y, res_tol,length(x), lsqr_tol,length(x));
    t=cputime-t;
    %infos=strcat(sprintf('\n\n'),l1infos);
    infos=strcat(infos,sprintf(' LSQR-Flag for computing 2-norm'));
    infos=strcat(infos,sprintf(' solution: %d. Time=%1.1e',lsqrflag,t));
  case 'linprog'
    y=[real(f);imag(f)];

    A = matrix;
%     A = exp(-2*pi*i*x*((-N/2):(N/2-1)));
    B=real(A);
    C=imag(A);

    t=cputime;
    linprog_options= optimset('LargeScale', 'on', 'TolFun',res_tol);
    u_hat = linprog(ones(2*N,1), -eye(2*N),zeros(2*N,1), [B, -B; C, -C],y,...
		    [],[],[],linprog_options);
  
    f_hat = u_hat(1:N) - u_hat(N+1:end);
    t=cputime-t;
    infos=strcat(infos,sprintf(' Time=%1.1e',t));
  otherwise
    disp('unknown method');
end;
