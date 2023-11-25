%% Demo of real Datasets
clear;
clc;
close all;
%% load dataset
% there are 9 synthetic datasets in the 'datasets' folder
name={'leuk72_3k','dermatology','eighthr','HillValley','phoneme'
    'X8D5K','libras','aggregation','unbalance',''};
name=name';
num=1;% num=1,2,3,...,9
data_name=char(name(num));
data1=load(['real\',data_name,'.mat']) ;
group1=load(['real\',data_name,'Group.mat']) ;
data=getfield (data1, 'data');
label=getfield (group1,'label');
label=clear_C(label);
K=size(unique(label),1);

%% clustering
% Provide two methods and choose one method

% Method 1: No parameters
C=RISM(data,1);

%     % Method 2: input cluster number
%     C=RISM(data,1);

%% Evaluation indicators
if(length(unique(C))<K)
    [AC,PR,RE,F1]=AC_PE_RE(C,label);
else
    [AC,PR,RE,F1]=AC_PE_RE(label, C);
end

