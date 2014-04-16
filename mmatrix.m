function Phi = mmatrix(varargin)

%--------------------------------------------------------------------------
% arguments:
% @m: # of rows
% @n: # of cols
% @m_type: type of matrix {
%     FULL -- uniformaly random entries
%     BLOCK -- normal distribution
%     BLOCKPLUS
%     SPRASE -- guassian distribution
% }
% @entry: entry distributon {
%     
% }
%--------------------------------------------------------------------------

narginchk(2, 3);

m = vargin(1);
n = vargin(2);
m_type = '';

if (nargin == 2)
    m_type = 'FULL';
elseif (nargin == 3)
    m_type = vargin(3);
end

if (strcmp(m_type, 'FULL'))
elseif (strcmp(m_type, 'BLOCK'))
elseif (strcmp(m_type, 'BLOCKPLUS'))
elseif (strcmp(m_type, 'SPARSE'))
end