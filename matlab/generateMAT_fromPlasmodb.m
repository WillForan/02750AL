
clear all;
%what file has all the values
labledCSV='../plasDB/labeled.csv';

%%%%%%%load features
%from csv
genes.feats  = csvread(labledCSV,1,1);
%set lables
genes.labels = genes.feats(:,end);
%remove lables from feats
genes.feats  = genes.feats(:,1:end-1);

%load feature names
fid = fopen(labledCSV, 'r'); tline = fgetl(fid);
genes.featnames = regexp(tline, '\,', 'split');

%remove first and last field (name and label)
genes.featnames = genes.featnames(2:end-1);

%load instance names
i=1;
while(~feof(fid))
   tline = fgetl(fid);
   line = regexp(tline, '\,', 'split');
   genes.names(i) = line(1);
   i=i+1;
end

save('genes.mat', 'genes');


