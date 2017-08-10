% input: filter, the static indicator for each sample obtained by a threshold,
% freq, how many samples are in one second
% output: consecutive static session recorded in QS_time_interval_info_matrix
% of size Nx3, and refined s_filter with small dynamic and static session
% outliers removed.
function [QS_time_interval_info_matrix, filter]=compute_time_interval_info_matrix(filter, freq)
l = 1;
QS_time_interval_info_matrix = zeros(1, 1 + 1 + 1); % record the static periods, unfortunately, dynamically changing size
samples = 0; % how many samples does this static period have
start = 0; % where does this static period starts
flag = 0; % what is the state of the last sample, 0 for dynamic, 1 for static

if filter(1) == 0    
    flag = 0;    
else    
    flag = 1;
    start = 1;    
end
% cycle to determine the QS_time_interval_info_matrix
for i = 1:length(filter)    
    if flag == 0 && filter(i) == 0
        0; % do nothing
    elseif flag == 1 && filter(i) == 1        
        samples = samples + 1;        
    elseif flag == 1 && filter(i) == 0 % take action when static chnages to dynamic        
        QS_time_interval_info_matrix(l, 1:3) = [start, i - 1, samples];
        l = l + 1;
        flag = 0;        
    elseif flag == 0 && filter(i) == 1 % take action when dynamic chnages to static        
        start = i;
        samples = 1;
        flag = 1;        
    end    
end
% Huai { refine QS_time_interval_info_matrix by two criteria,
% (1) if a dynamic period is too short, less than freq, set it as
% static, (2) if a static period is too short, set it as dynamic
% we prefer to set dynamic to static if (1) is met

% we do it several times, as some spikes may merge into another spike
% meeting the two criteria. There are much better ways to handle this
% problem than using loop, but for time sake, I don't do that.

moreIter=true; % do we need more iterations to remove duds
maxIter=10; % do not run too many times
howMany=0;
while(moreIter&& howMany<maxIter)
    total_sample=length(filter);
    for j = 1:size(QS_time_interval_info_matrix(:,1),1)
        if(j==1) % check first dynamic interval
            if(QS_time_interval_info_matrix(1,1)-1<freq)
                QS_time_interval_info_matrix(1,1)=1;
                QS_time_interval_info_matrix(1,3)=QS_time_interval_info_matrix(1,2);
            end
        elseif(j==size(QS_time_interval_info_matrix(:,1),1))
            if(total_sample-QS_time_interval_info_matrix(j,2)<freq)
                QS_time_interval_info_matrix(j,2)=total_sample;
                QS_time_interval_info_matrix(j,3)=QS_time_interval_info_matrix(j,2)-QS_time_interval_info_matrix(j,1)+1;
            end
        else
            if(QS_time_interval_info_matrix(j,1)-QS_time_interval_info_matrix(j-1,2)-1<freq)
                QS_time_interval_info_matrix(j-1,2)=QS_time_interval_info_matrix(j,2);
                QS_time_interval_info_matrix(j-1,3)=QS_time_interval_info_matrix(j-1,2)-QS_time_interval_info_matrix(j-1,1)+1;
                QS_time_interval_info_matrix(j,3)=0;
            end
        end
    end
    for j = 1:size(QS_time_interval_info_matrix(:,1),1)
        if(QS_time_interval_info_matrix(j,3)>0&&QS_time_interval_info_matrix(j,3)<freq)
            QS_time_interval_info_matrix(j,3)=0;
        end
    end
    % remove the periods of 0 samples
    removeIndex=QS_time_interval_info_matrix(:,3)==0; % which entry to remove
    QS_time_interval_info_matrix=QS_time_interval_info_matrix(QS_time_interval_info_matrix(:,3)~=0, :);
    howMany=howMany+1;
    if(isempty(removeIndex))
        moreIter=false;
    end
end
% refine filter and save in s_filter
filter=zeros(size(filter));
for j = 1:size(QS_time_interval_info_matrix(:,1),1)
    filter(1,QS_time_interval_info_matrix(j, 1):QS_time_interval_info_matrix(j, 2))=1;
end
% Huai }
end