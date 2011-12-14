%%%%%%%%%%%%%%%%%%
%
% PLOTS
%
%%%%%%%%%%%%%%%%%%
function plotIndividual(inputdata,usemethods)
    
    %Messsed up struct names, fixed in file name
    nameCorrect={}
    
    %pick what dataset to use
    if(exist('inputdata'))
	switch inputdata
	case {'PPI' 'ppi'}
	    loadPPI;
	case {'oldPF','oldgenes', 'PFnogc'}
	    loadGenes_nogc;
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
	nameCorrect{end+1}=alType;
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
	   alType=fieldnames(temp);
	   alType=alType{1};

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
	for b=1:length(AL.(alType).ontrain.batch) l=length(AL.(alType).ontrain.batch{b});
	    if l<m; m = l; end;
	end
	AL.(alType).ontrain.minIt = m;
	if m < minimumIterations; minimumIterations=m; end
    end

    usemethods=fieldnames(AL);
    for m=1:length(usemethods)

	clear A At NumP CumP Precision Recall;
	alType=usemethods{m};
	for fold=1:AL.(alType).validationFolds
	    TP  = AL.(alType).ontrain.TP{fold}(1:minimumIterations);
	    FP  = AL.(alType).ontrain.FP{fold}(1:minimumIterations);
	    TN  = AL.(alType).ontrain.TN{fold}(1:minimumIterations);
	    FN  = AL.(alType).ontrain.FN{fold}(1:minimumIterations);

	    TPt = AL.(alType).ontest.TP{fold}(1:minimumIterations);
	    FPt = AL.(alType).ontest.FP{fold}(1:minimumIterations);
	    TNt = AL.(alType).ontest.TN{fold}(1:minimumIterations);
	    FNt = AL.(alType).ontest.FN{fold}(1:minimumIterations);


	    A(fold,:)=(TP+TN)./(TP+FP+TN+FN);
	    At(fold,:)=(TPt+TNt)./(TPt+FPt+TNt+FNt);
	    Precision(fold,:)=(TPt)./(FPt+TPt);
	    Recall(fold,:)=(TPt)./(FNt+TPt);

	    %%if there ar eno positives, make fscore 0
	    Precision(fold, find(isnan(Precision(fold,:))==1))=0;

	    FScore(fold,:)=2.* Precision(fold,:).*Recall(fold,:) ./ ( Precision(fold,:) + Recall(fold,:) );

	    NumP(fold,:) = AL.(alType).numPos{fold}(1:minimumIterations);

	    %compensate for no cost in field
	    if(~isfield(AL.(alType),'cost'))
		fprintf('ERR -- no cost field in %s at fold %i\n',alType,fold);
		AL.(alType).cost{fold}=zeros(1,minimumIterations);
	    else 
	        if(size(AL.(alType).cost,2)<fold)
		    fprintf('ERR -- no cost for fold %i\n',fold);
		    AL.(alType).cost{fold}=zeros(1,minimumIterations);
		end
	    end

	    %size(AL.(alType).cost{fold})
	    Cost(fold,:) = AL.(alType).cost{fold}(1:minimumIterations);

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
	   sum(FScore)./AL.(alType).validationFolds 
	   sum(Cost)./AL.(alType).validationFolds  ...
	  ];
	 
	halfwaypt=find(Plot{m}(2,:) > .5);
	if(length(halfwaypt)==0); halfwaypt=[0];end
	fprintf('%s 1/2-pt\t%.3f\n',alType,halfwaypt(1));

	fprintf('%s $1/2pt\t%.3f\n',alType,sum(Plot{m}(8,1:halfwaypt(1))  ));
	totalCost=sum(Plot{m}(8,:));
	fprintf('%s total Cost\t%.3f\n',alType,totalCost);

	%%%%Plots for individual learners
	fig=figure;
	plot( [1:minimumIterations], Plot{m}([1 2 4 7],:), '-s');
	xlabel('Iteration'); ylabel('%'); 
	%legend({'+/batch','TP:Total', 'Train Accuracy','Test Accuracy' 'Precision' 'Recall'}, ...
	legend({'+/batch' 'TP:Total' 'Test Accuracy' 'FScore'}, ...
	      'Location', 'EastOutside');
	title([inputdata ' - ' AL.(alType).title ' - Batch Size ' num2str(AL.(alType).BatchSize) ...
	       '- ' num2str(AL.(alType).validationFolds) ' Folds']);
        
	hgexport(fig,['img/' ResultsFolder '-' AL.(alType).imgfile])



    end

    %allocate for interesting plots
    NumP   = zeros(length(usemethods),minimumIterations);
    CumP   = NumP;
    At     = NumP;
    Cost   = NumP;
    FScore = NumP;
    for m=1:length(usemethods)
	NumP(m,:)=Plot{m}(1,1:minimumIterations);
	CumP(m,:)=Plot{m}(2,1:minimumIterations);
	At(m,:)=Plot{m}(4,1:minimumIterations);
	FScore(m,:)=Plot{m}(7,1:minimumIterations);
	Cost(m,:)=Plot{m}(8,1:minimumIterations);
    end

    %NumP    %vaerge number of postives in each batch
    %Plot{1} %average % of postives in each batch

    %NumP(:,1:5)
    %At(:,1:5)
    %usemethods'
    %return

    set(0,'DefaultAxesLineStyleOrder','-|--|:')

    fig=figure;
    plot([1:minimumIterations],NumP)
    title('+/Batch');
    xlabel('Instances Seen'); ylabel('%'); 
    legend(nameCorrect,'Location', 'EastOutside');
    hgexport(fig,['img/' ResultsFolder '-PosPerBatch']);

    fig=figure;
    plot([1:minimumIterations],CumP)
    title('Cumulative +');
    xlabel('Instances Seen'); ylabel(''); 
    legend(nameCorrect,'Location', 'EastOutside');
    hgexport(fig,['img/' ResultsFolder '-Cumulative']);

    fig=figure;
    plot([1:minimumIterations],At)
    title('Test Accuracy');
    xlabel('Instances Seen'); ylabel(''); 
    legend(nameCorrect,'Location', 'EastOutside');
    hgexport(fig,['img/' ResultsFolder '-testAccuracies']);

    fig=figure;
    plot([1:minimumIterations],FScore)
    title('Test F-score');
    xlabel('Instances Seen'); ylabel(''); 
    legend(nameCorrect,'Location', 'EastOutside');
    hgexport(fig,['img/' ResultsFolder '-testAccuracies']);
    %no cost in PPI
    if(~strcmp(ResultsFolder,'PPIResults'))
	fig=figure;
	plot([1:minimumIterations],Cost)
	title('Cost');
	ylabel('$');
	xlabel('Instances Seen'); ylabel(''); 
	legend(usemethods,'Location', 'EastOutside');
	hgexport(fig,['img/' ResultsFolder '-cost']);
    end


end

    
