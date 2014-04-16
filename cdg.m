%% Compressive Data Gathering
function re=cdg(s,m,k)
n = length(s);
originData = s';

%DCT
[t, basis] = dct4(originData);

%%uniform random matrix
a = zeros(m,n);
for i = 1 : m
    for j = 1 : n
        a(i, j) = rand();
    end
end
% coh(a);

measureData = a * originData';

%% WOMP
% d=dis(n);
% r=womp(measureData, a*basis^(-1),n,k,d);
% re=basis*r';

%% OMP Recovery
estimatedSparseRecovery = omp(measureData, a * basis^(-1), n, k);
estimatedRecovery = basis^(-1) * estimatedSparseRecovery';
re=estimatedRecovery;

%% BP
% co = bp(measureData, a * basis^(-1), n);
% re=basis^(-1) *co;

% %Evaluation
% sensed_mse = 0;
% % sparse_mse = 0;
% estimatedOriginMse = 0;
% % exactOriginMse = 0;
% for i = 1 : n
% %     sensed_mse = sensed_mse + (estimatedSparseRecovery(i) - sensedData(i))^2;
% %     sparse_mse = sparse_mse + (exactSparseRecovery(i) - sparseData(i))^2;
%     estimatedOriginMse = estimatedOriginMse + (estimatedRecovery(i) - originData(i))^2;
% %     exactOriginMse = exactOriginMse + (exactRecovery(i) - originData(i))^2;
% end
% % sensed_mse = sensed_mse / n;
% % sparse_mse = sparse_mse / n;
% estimatedOriginMse = estimatedOriginMse / n;
% % exactOriginMse = exactOriginMse / n;
% % sensed_mse
% % sparse_mse
% estimatedOriginMse
% % exactOriginMse
% 
% % legend11 = ['Estimated Recovery MSE = ' num2str(estimatedOriginMse)];
% % legend12 = ['Exact Recovery MSE = ' num2str(exactOriginMse)];
% % 
% % legend21 = ['DCT sparsity = ' num2str(eltCount)];
% % legend22 = ['Estimated Sparsity = ' num2str(n * sparsity)];
% 
% %Output
% originDataPoint = [1 : n; originData];
% sparse = basis * originData';
% dctDataPoint = [1 : n; sparse'];
% % exactSparseRecoveryPoint = [1 : n; exactSparseRecovery];
% estimatedSparseRecoveryPoint = [1 : n; estimatedSparseRecovery];
% % exactRecoveryPoint = [1 : n; exactRecovery'];
% estimatedRecoveryPoint = [1 : n; estimatedRecovery'];
% % originDataPoint = [xvec; originData];
% % exactRecoveryPoint = [xvec; exactRecovery'];
% % estimatedRecoveryPoint = [xvec; estimatedRecovery'];
% figure;
% bar(1 : n, sparse');
% hold on;
% plot(1 : n, estimatedSparseRecovery, 'r+');
% 
% figure;
% plot(originDataPoint(1, 1 : n), originDataPoint(2, 1 : n), 'k.-');
% hold on;
% plot(estimatedRecoveryPoint(1, 1 : n), estimatedRecoveryPoint(2, 1 : n), 'r-');
% % plot(exactRecoveryPoint(1, 1 : basisLength), exactRecoveryPoint(2, 1 : basisLength), 'b--');
% title('Original & Recovered Data');
% % legend('Original Data', legend11, legend12, 'Location', 'SouthEast');
% % 
% % figure;
% % bar(dctDataPoint(1, 1 : basisLength), dctDataPoint(2, 1 : basisLength));
% % hold on;
% % plot(estimatedSparseRecoveryPoint(1, 1 : n), estimatedSparseRecoveryPoint(2, 1 : n), 'r+');
% % plot(exactSparseRecoveryPoint(1, 1 : n), exactSparseRecoveryPoint(2, 1 : n), 'ko');
% % title('DCT Representation & Sparse Recovery Data');
% % legend(legend21, legend22, 'Location', 'SouthEast');
