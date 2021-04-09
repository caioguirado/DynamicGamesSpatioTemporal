function [Totalclusters, IntervalCluster, numberofclusters] = SpatioTemporalDetectionOFRecurrence(Signals,L,WindowLength,Overlap,ThresholdSpectralConcentration,MinNumberofPointsInaRegion,MinProminence)


% inputs:
% Signals: matrix collecting the signal/time series from all points in the
% hyperplane under analysis (rows are points, and columns are samples)
% L: sparse adjacency matrix collecting information about proximity among
% the points in the hyperplane (each row is a point, and the entries are
% the labels/positions of the points near it)


if nargin < 3 || isempty(WindowLength)
    WindowLength = 10;
end

if nargin < 4 || isempty(Overlap)
    Overlap = floor(WindowLength/2);
end

if nargin < 5 || isempty(ThresholdSpectralConcentration)
    ThresholdSpectralConcentration = 0.1;
end

if nargin < 6 || isempty(MinNumberofPointsInaRegion)
    MinNumberofPointsInaRegion = 20;
end

if nargin < 7 || isempty(MinProminence)
    MinProminence = 0.6;
end

step = WindowLength-Overlap;

Totalclusters = [];
countiteration = 1;
numberofclusters = [];
countclusters = 0;
IntervalCluster = [];

for indexwindow = 1:step:size(Signals,2)
    indexwindow
    clear B B0 U S V A Z Yz countPC pcindex positions temppos ...
        clusters positions3 positions4 nonemptycellsinclusters positionstoremove ...
        count allpositions checkzero detposition 1i i11 i111 i1111 i2 i22 i222 i2222 ...
        ii ind ind1 ind2 ind3 index index0 index1 index2 
    
    if indexwindow+WindowLength-1 <= size(Signals,2)
    B = Signals(:,indexwindow:indexwindow+WindowLength-1);
    else
        B = Signals(:,indexwindow:end);
    end
    
%     B0 = bsxfun(@minus, B, mean(B,2));
    B0 = B;
    [U, S, V] = svd(B0,'econ');
    A = U*S/sqrt(size(B,1)-1);
    Z = sqrt(size(B,1)-1)*V';
    
    Yz = abs(fft(Z'));
    countPC = 0;
    pcindex = [];
    for ii = 1:size(Yz,2)
%         [~,l,w,~] = findpeaks(Yz(:,ii)/max(Yz(:,ii)),'MinPeakProminence',0.6,'Annotate','extents');
        sy = Yz(:,ii)/max(Yz(:,ii));
        [~,ltemp,wtemp,~] = findpeaks(sy,'MinPeakProminence',MinProminence,'WidthReference','halfheight');
        [~,indmaxp] = max(sy(ltemp));
        l = ltemp(indmaxp);
        w = wtemp(indmaxp);
        if l-round(w/2)>1 & sum(Yz(l-round(w/2):l+round(w/2),ii).^2)/sum(Yz(:,ii).^2) > ThresholdSpectralConcentration %l-round(w)>1 & sum(Yz(l-round(w):l+round(w),ii).^2)/sum(Yz(:,ii).^2) > ThresholdSpectralConcentration
            countPC = countPC+1;
            pcindex = [pcindex ii];
        elseif l-round(w/2)>1 & sum(Yz(l-round(w/2):l+round(w/2),ii).^2)/sum(Yz(:,ii).^2) < ThresholdSpectralConcentration & sum(Yz(l-round(w/2):l+round(w/2),ii).^2)/sum(Yz(:,ii).^2) > 0%l-round(w)>1 & sum(Yz(l-round(w):l+round(w),ii).^2)/sum(Yz(:,ii).^2) < ThresholdSpectralConcentration & sum(Yz(l-round(w):l+round(w),ii).^2)/sum(Yz(:,ii).^2) > 0
            break
        end
    end
    
    J = countPC;
    positions = [];
    for counter = 1:J
        positions = [positions;find(abs(A(:,pcindex(counter)))/max(abs(A(:,pcindex(counter))))>3*std(abs(A(:,pcindex(counter)))/max(abs(A(:,pcindex(counter))))))];
    end
    
    temppos = unique(positions);
    clear positions
    positions = temppos;
    
    positions(diff(positions)>3) = [];
    
%     MatrixWithOnlyImportantPositions = B(positions,:);%B(positions,:);
%     for checkzero = size(MatrixWithOnlyImportantPositions,1):-1:1
%         if mean(sum(MatrixWithOnlyImportantPositions(checkzero,:))) == 0
%             MatrixWithOnlyImportantPositions(checkzero,:) = [];
%             positions(checkzero) = [];
%         end
%     end
%     m = size(MatrixWithOnlyImportantPositions,2);
    
    
    % find connected points - identify regions
    clear clusters
    clusters = cell(1);
    positions2 = positions;
    count = 1;
    for ind1 = 1:length(positions)
        [i1111,~] = ismember(positions(ind1),positions2);
        if sum(i1111) > 0
            clear tempcluster
            [i111,i222] = ismember(L(positions(ind1),:),positions);
            tempcluster = [positions(ind1) positions(i222(i111))'];
            for ind2 = ind1:length(positions)
                [i1,~] = ismember(L(positions(ind2),:),tempcluster);
                if sum(i1) > 0
                    [i11,i22] = ismember(L(positions(ind2),:),positions);
                    if sum(i11) > 0
                        tempcluster = [tempcluster positions(i22(i11))' positions(ind2)];
                    end
                end
            end
            clusters{count} = unique(tempcluster);
            clear positions2
            clear allpositions
            allpositions = [];
            for detposition = 1:count
                allpositions = [allpositions; clusters{detposition}'];
            end
            positions2 = setdiff(positions,allpositions);
            count = count + 1;
        end
    end
    
    % remove small regions - possible artifacts
    lengthclusters = [];
    for ind3 = 1:count-1
        lengthclusters(ind3) = numel(clusters{ind3});
    end
    
    for index0 = length(lengthclusters):-1:1
        if lengthclusters(index0)<MinNumberofPointsInaRegion
            clusters{index0} = [];
            lengthclusters(index0) = [];
        end
    end
    
    
    positions3 = cell(1);
    nonemptycellsinclusters = find(~cellfun('isempty', clusters));
    for index = 1:length(nonemptycellsinclusters)
        %     clear i11 i22
        %     [i11,i22] = ismember(clusters{index},positions);
        positions3{index} = clusters{nonemptycellsinclusters(index)};
    end
    
    positionstoremove = [];
    for index1 = 1:length(positions3)
        clear i1 i2
        [i1,~] = ismember(index1,positionstoremove);
        if i1 == 0
            for index2 = index1+1:length(positions3)
                clear i11 i22
                [i11,~] = ismember(positions3{index1},positions3{index2});
                if sum(i11) > 0
                    positionstoremove = [positionstoremove index2];
                    positions3{index1} = unique(union(positions3{index1},positions3{index2}));
                end
            end
        end
    end
    
    positions4 = cell(1);
    countpos = 1;
    for index = 1:length(positions3)
        clear i1 i2
        [i1,~] = ismember(index,positionstoremove);
        if i1 == 0
            positions4{countpos} = positions3{index};
            countpos = countpos+1;
        end
    end
    
    if length(positions4) == 1
        if isempty(positions4{1})
            numberofclusters(countiteration) = 0;
            else
        numberofclusters(countiteration) = length(positions4);
        end
    else
        numberofclusters(countiteration) = length(positions4);
    end
    
    
    newclusters = [];
    if numberofclusters(countiteration) > 0 & countclusters == 0
        for checkclusters0 = 1:length(positions4)
            countclusters = countclusters + 1;
            Totalclusters{countclusters} = positions4{checkclusters0};
            IntervalCluster(countclusters,:) = [indexwindow 0];
        end
    elseif numberofclusters(countiteration) > 0 & countclusters > 0
        newclusters = positions4;
        % check in the current window for clusters already identified in the previous window
        for checkclusters1 = 1:length(Totalclusters)
            flagcluster = 0;
            for checkclusters2 = length(newclusters):-1:1
                [i1,i2] = ismember(Totalclusters{checkclusters1},newclusters{checkclusters2});
                if sum(i1) > 0
                    Totalclusters{checkclusters1} = union(Totalclusters{checkclusters1},newclusters{checkclusters2});
                    newclusters{checkclusters2} = [];
                    flagcluster = 1;
                    break
                end
            end
            if flagcluster == 0 & IntervalCluster(checkclusters1,2) == 0 %%%
                IntervalCluster(checkclusters1,2) = indexwindow; %%%
            end
        end
    end
    
    if numberofclusters(countiteration) == 0 & countclusters > 0
        IntervalCluster(IntervalCluster(:,2) == 0,2) = indexwindow;
    end
        
        % check for new clusters in the current window
        
        for checkclusters3 = 1:length(newclusters)
            if ~isempty(newclusters{checkclusters3})
                countclusters = countclusters + 1;
                Totalclusters{countclusters} = newclusters{checkclusters3};
                IntervalCluster(countclusters,:) = [indexwindow 0];
            end
        end
        
        countiteration = countiteration + 1;
        
end
    
if ~isempty(IntervalCluster)
IntervalCluster(IntervalCluster(:,2) == 0,2) = size(Signals,2);
end

%%%
for index1 = 1:length(Totalclusters)
    poscluster = [];
    for index2 = index1+1:length(Totalclusters)
        clear i11
        [i11,~] = ismember(Totalclusters{index2},Totalclusters{index1});
        if sum(i11) > 0
            poscluster = [poscluster index2];
        end
    end
    for index3 = 1:length(poscluster)
        Totalclusters{index1} = unique(union(Totalclusters{index1},Totalclusters{poscluster(index3)}));
    end
    if ~isempty(poscluster)
        IntervalCluster(index1,:) = [IntervalCluster(index1,1) max(IntervalCluster([index1 poscluster],2))];
    end
    for index3 = 1:length(poscluster)
        Totalclusters{poscluster(index3)} = [];
    end  
end
         
tempTcluster = Totalclusters;
tempinterval = IntervalCluster;
clear Totalclusters IntervalCluster
countT = 1;
for index = 1:length(tempTcluster)
    if ~isempty(tempTcluster{index})
        Totalclusters{countT} = tempTcluster{index};
        IntervalCluster(countT,:) = tempinterval(index,:);
        countT = countT + 1;
    end
end
