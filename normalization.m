function [data_n]=normalization(data)
    max_data=max(data);
    min_data=min(data);
    max_data=max_data(ones(size(data,1),1),:);
    min_data=min_data(ones(size(data,1),1),:);
    data_n=(data-min_data)./(max_data-min_data+0.00001);
end
