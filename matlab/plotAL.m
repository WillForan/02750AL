%%%%%%%%%%%%%%%%%%
%
% PLOTS
%
%%%%%%%%%%%%%%%%%%
function plotIndividual(usemethods,genes)
    
    %folder='alResults';
    folder='PPIResults';

    if(~exist('usemethods')) 
	%load all .mat files in results folder
	cd folder;
	usemethods=ls('*mat');
	cd '..';

	usemethods=strrep(usemethods, '.mat', '');
	usemethods=textscan(usemethods,['%s' ' ']);
	usemethods=usemethods{1}.';
    end
    if(~exist('genes'))
	switch folder
	case 'PPIResults'
	    load 'PPI.mat';
	case 'alResults'
	    load 'genes.mat';
	otherwise
	    fprintf('I dont know what to do with folder %s\n',folder);
    end
    

    minimumIterations=Inf;

    %build AL
    for m=1:length(usemethods)
	alType=usemethods{m};
	temp = load([folder '/' (alType) '.mat']);
	AL.(alType)=temp.(alType);
	clear temp;

	%%get the minium number of iterterations
	m=Inf;
	for b=1:length(AL.(alType).ontrain.batch)
	    l=length(AL.(alType).ontrain.batch{b});
	    if l<m; m = l; end;
	end
	AL.(alType).ontrain.minIt = m;
	if m < minimumIterations; minimumIterations=m; end
    end

    for m=1:length(usemethods)

	clear A At NumP CumP Precision Recall;
	alType=usemethods{m};
	for fold=1:AL.(alType).validationFolds
	    TP  = AL.(alType).ontrain.TP{fold}(1:AL.(alType).ontrain.minIt);
	    FP  = AL.(alType).ontrain.FP{fold}(1:AL.(alType).ontrain.minIt);
	    TN  = AL.(alType).ontrain.TN{fold}(1:AL.(alType).ontrain.minIt);
	    FN  = AL.(alType).ontrain.FN{fold}(1:AL.(alType).ontrain.minIt);

	    TPt = AL.(alType).ontest.TP{fold}(1:AL.(alType).ontrain.minIt);
	    FPt = AL.(alType).ontest.FP{fold}(1:AL.(alType).ontrain.minIt);
	    TNt = AL.(alType).ontest.TN{fold}(1:AL.(alType).ontrain.minIt);
	    FNt = AL.(alType).ontest.FN{fold}(1:AL.(alType).ontrain.minIt);


	    A(fold,:)=(TP+TN)./(TP+FP+TN+FN);
	    At(fold,:)=(TPt+TNt)./(TPt+FPt+TNt+FNt);
	    Precision(fold,:)=(TPt)./(FPt+TPt);
	    Recall(fold,:)=(TPt)./(FNt+TPt);

	    NumP(fold,:) = AL.(alType).numPos{fold}(1:AL.(alType).ontrain.minIt);

	    totalPos=length(find(genes.labels(AL.(alType).training{fold})==1));
	    for i=1:length(NumP)
		CumP(fold,i)=sum(NumP(fold,1:i))./totalPos;
	    end
	end

	Plot{m} =  [  ...
	   sum(NumP)./AL.(alType).validationFolds ./AL.(alType).BatchSize
	   sum(CumP)./AL.(alType).validationFolds 
	   sum(A)./AL.(alType).validationFolds
	   sum(At)./AL.(alType).validationFolds
	   sum(Precision)./AL.(alType).validationFolds
	   sum(Recall)./AL.(alType).validationFolds
	  ];
	 
	halfwaypt=find(Plot{m}(2,:) > .5);
	fprintf('%s 1/2-pt\t%i\n',alType,halfwaypt(1));

	fig=figure;
	plot( [1:AL.(alType).ontrain.minIt], Plot{m}, '-s');
	xlabel('Iteration'); ylabel('%'); 
	legend({'+/batch','TP:Total', 'Train Accuracy','Test Accuracy' 'Precision' 'Recall'}, ...
	      'Location', 'NorthOutside');
	title([ AL.(alType).title ' - Batch Size ' num2str(AL.(alType).BatchSize) ...
	       '- ' num2str(AL.(alType).validationFolds) ' Folds']);

	hgexport(fig,['img/' folder '-' AL.(alType).imgfile])

    end

    for m=length(usemethods)
	minimumIterations 
	size(NumP)
	NumP(m,:)=Plot{m}(1,1:minimumIterations);
	CumP(m,:)=Plot{m}(2,1:minimumIterations);
	At(m,:)=Plot{m}(4,1:minimumIterations);
    end

    fig=figure;
    plot([1:minimumIterations],NumP)
    title('+/Batch');
    xlabel('Instances Seen'); ylabel('%'); 
    legend(usemethods,'Location', 'EastOutside');
    hgexport(fig,['img/' folder '-PosPerBatch']);

    fig=figure;
    plot([1:minimumIterations],CumP)
    title('Cumulative +');
    xlabel('Instances Seen'); ylabel(''); 
    legend(usemethods,'Location', 'EastOutside');
    hgexport(fig,['img/' folder '-Cumulative']);

    fig=figure;
    plot([1:minimumIterations],At)
    title('Test Accuracy');
    xlabel('Instances Seen'); ylabel(''); 
    legend(usemethods,'Location', 'EastOutside');
    hgexport(fig,['img/' folder '-testAccuracies']);


end

    
