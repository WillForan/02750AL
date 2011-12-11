%%%%%%%%%%%%%%%%%%
%
% Graph
%
%%%%%%%%%%%%%%%%%%
function graphs(usemethods,c)
    if(~exist('usemethods')) 
	%load all .mat files in results folder
	cd 'alResults';
	usemethods=ls('*mat');
	cd '..';

	usemethods=strrep(usemethods, '.mat', '');
	usemethods=textscan(usemethods,['%s' ' ']);
	usemethods=usemethods{1}.';
    end
    if(~exist('genes'))
	load 'genes.mat';
    end
    

    for ii=1:length(usemethods)
	%%%
	% get useful variables
	%%%
	alType=usemethods{ii};
	temp = load(['alResults/' (alType) '.mat']);
	AL.(alType)=temp.(alType);
	clear temp;

	aTP = sum(AL.(alType).ontrain.TP ) ./ AL.(alType).validationFolds;
	aTN = sum(AL.(alType).ontrain.TN ) ./ AL.(alType).validationFolds;
	aFP = sum(AL.(alType).ontrain.FP ) ./ AL.(alType).validationFolds;
	aFN = sum(AL.(alType).ontrain.FN ) ./ AL.(alType).validationFolds;


	bTP = sum(AL.(alType).ontest.TP ) ./ AL.(alType).validationFolds;
	bTN = sum(AL.(alType).ontest.TN ) ./ AL.(alType).validationFolds;
	bFP = sum(AL.(alType).ontest.FP ) ./ AL.(alType).validationFolds;
	bFN = sum(AL.(alType).ontest.FN ) ./ AL.(alType).validationFolds;

	%build True Positve ratio
	for f=1:AL.(alType).validationFolds;
	    TPr(f,:) = AL.(alType).ontrain.TP(f,:) ./ length(find(genes.labels(AL.(alType).training{f})==1));
	end

	%%%%
	%build plot data 
	%%%%
	TPr(ii,:) = sum(TPr) ./ AL.(alType).validationFolds;

	if(TPr(ii,:)>1); fprintf('%s bigger than 1!(%i)\n',alType,AL.(alType).validationFolds); end;

	trainAccuracy(ii,:) = (aTP + aTN)./(aTP + aTN + aFP + aFN);
	testAccuracy(ii,:)  = (bTP + bTN)./(bTP + bTN + bFP + bFN);
	%TP./(TP+FN) ... Sensitivity
	%TP./(TN+FP) ... Sepecificity
	%FP./(FP+TN) ... 

	%%%
	%
	% Plots
	%
	%%%%%
	instances=[ AL.(alType).numInstStart:AL.(alType).BatchSize:AL.(alType).BatchSize*(length(AL.(alType).ontest.TP)-1)+AL.(alType).numInstStart];
	fig=figure;
	%plot ([1:maxIterations] .* BatchSize,sum(AL.(alType).accuracy) ./ c.validationFolds);
	plot ( AL.(alType).ontrain.batch(1:size(TPr,2)), ...
	      [ TPr(ii,:) 
	        trainAccuracy(ii,:)
	        testAccuracy(ii,:) ...
		], ...
		'-s' ...
	     );
	

	legend({'TruePos:Total','Train Accuracy','Test Accuracy'},'Location', 'NorthOutside');
	title([ AL.(alType).title ' - Batch Size ' num2str(AL.(alType).BatchSize) '- ' num2str(AL.(alType).validationFolds) ' Folds']);
	xlabel('Instances Seen'); ylabel(''); 
	hgexport(fig,['img/' AL.(alType).imgfile])

    end

    if length(usemethods)<2
	return
	%don't overwrite comparison plots if nothing to compare
    end

    %assume all plots are scaled the same with same batch numbers.. eww
    fig=figure;
    %plot(instances,trainAccuracy)
    plot(1:length(trainAccuracy),trainAccuracy)
    title('Train Accuracies');
    xlabel('Instances Seen'); ylabel(''); 
    legend(usemethods,'Location', 'NorthOutside');
    hgexport(fig,'img/trainAccuracies');

    fig=figure;
    plot(1:length(testAccuracy),testAccuracy)
    %plot(instances,testAccuracy)
    title('Test Accuracies');
    xlabel('Instances Seen'); ylabel(''); 
    legend(usemethods,'Location', 'NorthOutside');
    hgexport(fig,'img/testAccuracies');

    fig=figure;
    plot(1:length(TPr),TPr)
    %plot(instances, TPr)
    title('True Positve Ratio');
    xlabel('Instances Seen'); ylabel(''); 
    legend(usemethods,'Location', 'NorthOutside');
    hgexport(fig,'img/truePos');

end
