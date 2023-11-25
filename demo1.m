%% Demo of Synthetic Datasets
clear;
clc;
close all;

%% load dataset
% there are 15 synthetic datasets in the 'datasets' folder
num=3;% num=1,2,3,...,15
data_name=num2str(num);
data1=load(['datasets\ds',data_name,'.mat']) ;
data=getfield (data1, 'data');
label=getfield (data1, 'label');
if(min(label(:))==0)
    label=label+1;
end
K=size(unique(label),1);

%% clustering
% Provide two methods and choose one method

% Method 1: No parameters
C=RISM(data,1);

%     % Method 2: input cluster number
%     C=RISM(data,1);

%% show results
img=figure;
gscatter(data(:,1),data(:,2),C);
legend('off');


