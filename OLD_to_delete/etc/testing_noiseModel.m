% Here doing it as well in a different way:
load('SC_HCP_cortexOnly.mat');load('scz_hcp_ten.mat');
for j=1:10;ICA_AROMA_2P_corr(:,:,j) = corr(time_series_hcp(:,:,j,1).');end;FC_emp = mean(ICA_AROMA_2P_corr,3);

lambda = -0.45/10;
for j=1:size(centroids,1)
    vertexCoords=centroids(j,:);
    keepPoints = 1:360;
    allDists = sqrt(((vertexCoords(1) - centroids(keepPoints,1)).^2 + (vertexCoords(2) - centroids(keepPoints,2)).^2 + (vertexCoords(3) - centroids(keepPoints,3)).^2));
    Prob = rand(size(allDists));
    ProbExponential = exp(lambda*allDists);
    matrix_rand_dist(j,keepPoints) = (Prob<ProbExponential);%.*ProbExponential;
    % matrix_rand_dist(j,keepPoints) = (Prob<ProbExponential).*allDists;
end
% figure;imagesc(log10(matrix_rand_dist))

AdjDens = matrix_rand_dist/(max(matrix_rand_dist(:))) + AdjDens;
% AdjDens = AdjCount;
C  = AdjDens/max(AdjDens(:));
C = 0.5*(C + C.');

C = circshift(C,180) + C;


figure;
hold on;
% nt2 = [100,300];
nt = 500;
% alpha = [0:0.25:5];
alpha = linspace(0,2,100);
% alpha = 0.2;
% nt2 = 500;
counter = 0;
evs = 5;
% for nt=nt2;
for a = alpha
	counter = 1 + counter;
	% nt = 300;
	D = sum(C,1);
	% GL = 0.5*(D-C).'.*(D-C);
	% [V,diagev] = eig(GL);
	% dev = diag(diagev);
	% diagev = diagev/dev(end-1);
	D2 = D;
	% D = V(:,end-3)';

	% signal = randn(nt,1);
	clear bold_simple
	% slope = linspace(-5,5,100);
	slope = 0;
	% for degWt = 1:length(slope)
		signal = randn(nt,1);
		% signal = movmean(signal,10);

		for n=1:length(D);
			% bold_simple(:,n) = (D2(n,end-evs:end-1)*(signal.')).' + a*randn(size(signal,1),1);
			bold_simple(:,n) = D(n)*signal + a*randn(size(signal,1),1);
		end
		cor_test = corr(bold_simple);
		% hh = triu(FC_emp,1);
		inds = find(triu(C,1)>0);
		% CC = corrcoef(cor_test(:),FC_emp(:));
		CC = corrcoef(cor_test(inds),FC_emp(inds));
		fitting(counter) = CC(2);
	% end
	
	% M{counter} = ['a=' num2str(a)];
end
plot(alpha,fitting,'.','MarkerSize',15)

% legend(M);
xlabel('Alpha');ylabel('Similarity (correlation)');
set(gca,'fontSize',18);


% % Fixing up C:
% matrix_rand_dist = matrix;

