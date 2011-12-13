%%%%%%%%%%%%%%%%%%
%
% PLOTS
%
%%%%%%%%%%%%%%%%%%
function plotIndividual(inputdata,usemethods)
    
    %pick what dataset to use
    if(exist('inputdata'))
	switch inputdata
	case {'PPI' 'ppi'}
	    loadPPI;
	    fprintf('using PPI')
	case {'PF','genes'}
	    loadGenes;
	otherwise
	    loadGenes;
	end
    else
	loadGenes;
    end

    if(~exist('usemethods')) 
	%load all .mat files in results folder
	cd(ResultsFolder);
	usemethods=ls('*mat');
	cd '..';

	usemethods=strrep(usemethods, '.mat', '');
	usemethods=textscan(usemethods,['%s' ' ']);
	usemethods=usemethods{1}';
	%trim off all numbers
	%assume names do not contain numbers, otherwise loading breaks everything
	for i=1:length(usemethods)
	  usemethods{i}(ismember(usemethods{i},'1234567890')) = [];
	end
	usemethods=unique(usemethods); %remove all duplicates 
    end
    

    minimumIterations=Inf;

    %build AL
    AL=struct;
    for m=1:length(usemethods)
	alType=usemethods{m};
	%%%%%%%% LOAD all different random sets
	list=ls([ResultsFolder '/*' (alType) '.mat']);

	list=textscan(list,['%s' ' ']);
	list=list{1}';
	numTrials=length(list);
	for f=1:numTrials
	   %%%Average all values
	   clear temp;
	   fprintf('loading %s\n',list{f});
	   temp = load(list{f});

	   %%%%% see values
	   %disp('TP');
	   %cat(1,temp.(alType).ontrain.TP{:})
	   %disp('TF');
	   %cat(1,temp.(alType).ontrain.FP{:})
	   %disp('TN');
	   %cat(1,temp.(alType).ontrain.TN{:})
	   %disp('FN');
	   %cat(1,temp.(alType).ontrain.FN{:})
	   %disp('numPos');
	   %cat(1,temp.(alType).numPos{:})
	   %%%%%

	   %if nothing to add, create
	   if(~isfield(AL,alType))
	     AL.(alType)=temp.(alType);
	   else
	     %check that they are comparable
	     if( AL.(alType).validationFolds ~= temp.(alType).validationFolds ...
	        || AL.(alType).BatchSize~= temp.(alType).BatchSize)
	       fprintf('Sets not comparable! (skipping  %s)\n',list{f});
	       continue;
	     end

	     %add TP,... sets together, then divide by total number of trials
	     for field={'TP','TN','FP','FN'}; for dset={'ontrain','ontest'}; for fold=1:temp.(alType).validationFolds
		 AL.(alType).(dset{1}).(field{1}){fold}=AL.(alType).(dset{1}).(field{1}){fold} + temp.(alType).(dset{1}).(field{1}){fold};
	     end;end;end;

	     %keep the batch size of the first trail
	     %should average? compile larger list? used to compute cumulative TP
	     %last will differ by however different cvpartion makes partions
	     %last TP:total will be off
	     for fold=1:temp.(alType).validationFolds
		 AL.(alType).numPos{fold}= AL.(alType).numPos{fold} + temp.(alType).numPos{fold};
	     end;

	   end
	end

	%AVERAGE RESULTS
        for field={'TP','TN','FP','FN'}; for dset={'ontrain','ontest'}; for fold=1:temp.(alType).validationFolds
            AL.(alType).(dset{1}).(field{1}){fold}=AL.(alType).(dset{1}).(field{1}){fold}./numTrials;
        end;end;end;
        for fold=1:temp.(alType).validationFolds
            AL.(alType).numPos{fold}=AL.(alType).numPos{fold}./numTrials;
        end;
	
       %%%%% see values
       %disp('TP');
       %cat(1,AL.(alType).ontrain.TP{:})
       %disp('TF');
       %cat(1,AL.(alType).ontrain.FP{:})
       %disp('TN');
       %cat(1,AL.(alType).ontrain.TN{:})
       %disp('FN');
       %cat(1,AL.(alType).ontrain.FN{:})
       %disp('numPos');
       %cat(1,AL.(alType).numPos{:})
       %%%%%
	

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
	if(length(halfwaypt)==0); halfwaypt=[0];end
	fprintf('%s 1/2-pt\t%i\n',alType,halfwaypt(1));

	%%%%Plots for individual learners
	%fig=figure;
	%plot( [1:AL.(alType).ontrain.minIt], Plot{m}, '-s');
	%xlabel('Iteration'); ylabel('%'); 
	%legend({'+/batch','TP:Total', 'Train Accuracy','Test Accuracy' 'Precision' 'Recall'}, ...
	%      'Location', 'EastOutside');
	%title([inputdata ' - ' AL.(alType).title ' - Batch Size ' num2str(AL.(alType).BatchSize) ...
	%       '- ' num2str(AL.(alType).validationFolds) ' Folds']);

	%hgexport(fig,['img/' ResultsFolder '-' AL.(alType).imgfile])

    end

    %allocate for interesting plots
    NumP = zeros(length(usemethods),minimumIterations);
    CumP = NumP;
    At   = NumP;
    for m=1:length(usemethods)
	NumP(m,:)=Plot{m}(1,1:minimumIterations);
	CumP(m,:)=Plot{m}(2,1:minimumIterations);
	At(m,:)=Plot{m}(4,1:minimumIterations);
    end
    NumP    %vaerge number of postives in each batch
    %Plot{1} %average % of postives in each batch

    fig=figure;
    plot([1:minimumIterations],NumP)
    title('+/Batch');
    xlabel('Instances Seen'); ylabel('%'); 
    legend(usemethods,'Location', 'EastOutside');
    hgexport(fig,['img/' ResultsFolder '-PosPerBatch']);

    fig=figure;
    plot([1:minimumIterations],CumP)
    title('Cumulative +');
    xlabel('Instances Seen'); ylabel(''); 
    legend(usemethods,'Location', 'EastOutside');
    hgexport(fig,['img/' ResultsFolder '-Cumulative']);

    fig=figure;
    plot([1:minimumIterations],At)
    title('Test Accuracy');
    xlabel('Instances Seen'); ylabel(''); 
    legend(usemethods,'Location', 'EastOutside');
    hgexport(fig,['img/' ResultsFolder '-testAccuracies']);


end

    
