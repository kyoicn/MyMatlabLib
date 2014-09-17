function [omp_result, support, residual_error] =StOMP(Sc, dictionary_matrix, S, Pfa)
% [omp_result, support, residual_error] =StOMP(Sc, dictionary_matrix, S, Pfa)
% Sc��۲�
% dictionary_matrix���걸����
% iteration���������
% Sϡ���
residual_error=Sc;%�в�
support = [];
omp_result=zeros(size(dictionary_matrix,2),1);
Tfa = sqrt(2*log10(1/Pfa));
for k=1:10  %����10��                                
    support_0 = [];
    inner_product=dictionary_matrix'*residual_error;
    nrm = norm(inner_product);
%    nrm = norm(residual_error);
    Threshold = Tfa*nrm/sqrt(length(inner_product));
    amplitude = abs(inner_product);
    for kk=1:length(amplitude)
        if amplitude(kk) > Threshold
            support_0 = union(support_0,kk);
        end
    end
    support=union(support,support_0);
    R = pinv(dictionary_matrix(:,support))*Sc;
    residual_error=Sc-dictionary_matrix(:,support)*R;
    if isempty(support_0)
        break;
    end
end
omp_result(support)=R;
%omp_result(support)=dictionary_matrix(:,support)'*Sc';

