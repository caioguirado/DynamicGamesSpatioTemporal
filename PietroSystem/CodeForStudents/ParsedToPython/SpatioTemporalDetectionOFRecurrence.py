# Generated with SMOP  0.41-beta
from libsmop import *
# SpatioTemporalDetectionOFRecurrence.m

    
@function
def SpatioTemporalDetectionOFRecurrence(Signals=None,L=None,WindowLength=None,Overlap=None,ThresholdSpectralConcentration=None,MinNumberofPointsInaRegion=None,MinProminence=None,*args,**kwargs):
    varargin = SpatioTemporalDetectionOFRecurrence.varargin
    nargin = SpatioTemporalDetectionOFRecurrence.nargin

    # inputs:
# Signals: matrix collecting the signal/time series from all points in the
# hyperplane under analysis (rows are points, and columns are samples)
# L: sparse adjacency matrix collecting information about proximity among
# the points in the hyperplane (each row is a point, and the entries are
# the labels/positions of the points near it)
    
    if nargin < 3 or isempty(WindowLength):
        WindowLength=10
# SpatioTemporalDetectionOFRecurrence.m:13
    
    if nargin < 4 or isempty(Overlap):
        Overlap=floor(WindowLength / 2)
# SpatioTemporalDetectionOFRecurrence.m:17
    
    if nargin < 5 or isempty(ThresholdSpectralConcentration):
        ThresholdSpectralConcentration=0.1
# SpatioTemporalDetectionOFRecurrence.m:21
    
    if nargin < 6 or isempty(MinNumberofPointsInaRegion):
        MinNumberofPointsInaRegion=20
# SpatioTemporalDetectionOFRecurrence.m:25
    
    if nargin < 7 or isempty(MinProminence):
        MinProminence=0.6
# SpatioTemporalDetectionOFRecurrence.m:29
    
    step=WindowLength - Overlap
# SpatioTemporalDetectionOFRecurrence.m:32
    Totalclusters=[]
# SpatioTemporalDetectionOFRecurrence.m:34
    countiteration=1
# SpatioTemporalDetectionOFRecurrence.m:35
    numberofclusters=[]
# SpatioTemporalDetectionOFRecurrence.m:36
    countclusters=0
# SpatioTemporalDetectionOFRecurrence.m:37
    IntervalCluster=[]
# SpatioTemporalDetectionOFRecurrence.m:38
    for indexwindow in arange(1,size(Signals,2),step).reshape(-1):
        indexwindow
        clear('B','B0','U','S','V','A','Z','Yz','countPC','pcindex','positions','temppos','clusters','positions3','positions4','nonemptycellsinclusters','positionstoremove','count','allpositions','checkzero','detposition','1j','i11','i111','i1111','i2','i22','i222','i2222','ii','ind','ind1','ind2','ind3','index','index0','index1','index2')
        if indexwindow + WindowLength - 1 <= size(Signals,2):
            B=Signals(arange(),arange(indexwindow,indexwindow + WindowLength - 1))
# SpatioTemporalDetectionOFRecurrence.m:48
        else:
            B=Signals(arange(),arange(indexwindow,end()))
# SpatioTemporalDetectionOFRecurrence.m:50
        #     B0 = bsxfun(@minus, B, mean(B,2));
        B0=copy(B)
# SpatioTemporalDetectionOFRecurrence.m:54
        U,S,V=svd(B0,'econ',nargout=3)
# SpatioTemporalDetectionOFRecurrence.m:55
        A=dot(U,S) / sqrt(size(B,1) - 1)
# SpatioTemporalDetectionOFRecurrence.m:56
        Z=dot(sqrt(size(B,1) - 1),V.T)
# SpatioTemporalDetectionOFRecurrence.m:57
        Yz=abs(fft(Z.T))
# SpatioTemporalDetectionOFRecurrence.m:59
        countPC=0
# SpatioTemporalDetectionOFRecurrence.m:60
        pcindex=[]
# SpatioTemporalDetectionOFRecurrence.m:61
        for ii in arange(1,size(Yz,2)).reshape(-1):
            #         [~,l,w,~] = findpeaks(Yz(:,ii)/max(Yz(:,ii)),'MinPeakProminence',0.6,'Annotate','extents');
            sy=Yz(arange(),ii) / max(Yz(arange(),ii))
# SpatioTemporalDetectionOFRecurrence.m:64
            __,ltemp,wtemp,__=findpeaks(sy,'MinPeakProminence',MinProminence,'WidthReference','halfheight',nargout=4)
# SpatioTemporalDetectionOFRecurrence.m:65
            __,indmaxp=max(sy(ltemp),nargout=2)
# SpatioTemporalDetectionOFRecurrence.m:66
            l=ltemp(indmaxp)
# SpatioTemporalDetectionOFRecurrence.m:67
            w=wtemp(indmaxp)
# SpatioTemporalDetectionOFRecurrence.m:68
            if l - round(w / 2) > logical_and(1,sum(Yz(arange(l - round(w / 2),l + round(w / 2)),ii) ** 2) / sum(Yz(arange(),ii) ** 2)) > ThresholdSpectralConcentration:
                countPC=countPC + 1
# SpatioTemporalDetectionOFRecurrence.m:70
                pcindex=concat([pcindex,ii])
# SpatioTemporalDetectionOFRecurrence.m:71
            else:
                if l - round(w / 2) > logical_and(1,sum(Yz(arange(l - round(w / 2),l + round(w / 2)),ii) ** 2) / sum(Yz(arange(),ii) ** 2)) < logical_and(ThresholdSpectralConcentration,sum(Yz(arange(l - round(w / 2),l + round(w / 2)),ii) ** 2) / sum(Yz(arange(),ii) ** 2)) > 0:
                    break
        J=copy(countPC)
# SpatioTemporalDetectionOFRecurrence.m:77
        positions=[]
# SpatioTemporalDetectionOFRecurrence.m:78
        for counter in arange(1,J).reshape(-1):
            positions=concat([[positions],[find(abs(A(arange(),pcindex(counter))) / max(abs(A(arange(),pcindex(counter)))) > dot(3,std(abs(A(arange(),pcindex(counter))) / max(abs(A(arange(),pcindex(counter)))))))]])
# SpatioTemporalDetectionOFRecurrence.m:80
        temppos=unique(positions)
# SpatioTemporalDetectionOFRecurrence.m:83
        clear('positions')
        positions=copy(temppos)
# SpatioTemporalDetectionOFRecurrence.m:85
        positions[diff(positions) > 3]=[]
# SpatioTemporalDetectionOFRecurrence.m:87
        #     MatrixWithOnlyImportantPositions = B(positions,:);#B(positions,:);
#     for checkzero = size(MatrixWithOnlyImportantPositions,1):-1:1
#         if mean(sum(MatrixWithOnlyImportantPositions(checkzero,:))) == 0
#             MatrixWithOnlyImportantPositions(checkzero,:) = [];
#             positions(checkzero) = [];
#         end
#     end
#     m = size(MatrixWithOnlyImportantPositions,2);
        # find connected points - identify regions
        clear('clusters')
        clusters=cell(1)
# SpatioTemporalDetectionOFRecurrence.m:101
        positions2=copy(positions)
# SpatioTemporalDetectionOFRecurrence.m:102
        count=1
# SpatioTemporalDetectionOFRecurrence.m:103
        for ind1 in arange(1,length(positions)).reshape(-1):
            i1111,__=ismember(positions(ind1),positions2,nargout=2)
# SpatioTemporalDetectionOFRecurrence.m:105
            if sum(i1111) > 0:
                clear('tempcluster')
                i111,i222=ismember(L(positions(ind1),arange()),positions,nargout=2)
# SpatioTemporalDetectionOFRecurrence.m:108
                tempcluster=concat([positions(ind1),positions(i222(i111)).T])
# SpatioTemporalDetectionOFRecurrence.m:109
                for ind2 in arange(ind1,length(positions)).reshape(-1):
                    i1,__=ismember(L(positions(ind2),arange()),tempcluster,nargout=2)
# SpatioTemporalDetectionOFRecurrence.m:111
                    if sum(i1) > 0:
                        i11,i22=ismember(L(positions(ind2),arange()),positions,nargout=2)
# SpatioTemporalDetectionOFRecurrence.m:113
                        if sum(i11) > 0:
                            tempcluster=concat([tempcluster,positions(i22(i11)).T,positions(ind2)])
# SpatioTemporalDetectionOFRecurrence.m:115
                clusters[count]=unique(tempcluster)
# SpatioTemporalDetectionOFRecurrence.m:119
                clear('positions2')
                clear('allpositions')
                allpositions=[]
# SpatioTemporalDetectionOFRecurrence.m:122
                for detposition in arange(1,count).reshape(-1):
                    allpositions=concat([[allpositions],[clusters[detposition].T]])
# SpatioTemporalDetectionOFRecurrence.m:124
                positions2=setdiff(positions,allpositions)
# SpatioTemporalDetectionOFRecurrence.m:126
                count=count + 1
# SpatioTemporalDetectionOFRecurrence.m:127
        # remove small regions - possible artifacts
        lengthclusters=[]
# SpatioTemporalDetectionOFRecurrence.m:132
        for ind3 in arange(1,count - 1).reshape(-1):
            lengthclusters[ind3]=numel(clusters[ind3])
# SpatioTemporalDetectionOFRecurrence.m:134
        for index0 in arange(length(lengthclusters),1,- 1).reshape(-1):
            if lengthclusters(index0) < MinNumberofPointsInaRegion:
                clusters[index0]=[]
# SpatioTemporalDetectionOFRecurrence.m:139
                lengthclusters[index0]=[]
# SpatioTemporalDetectionOFRecurrence.m:140
        positions3=cell(1)
# SpatioTemporalDetectionOFRecurrence.m:145
        nonemptycellsinclusters=find(logical_not(cellfun('isempty',clusters)))
# SpatioTemporalDetectionOFRecurrence.m:146
        for index in arange(1,length(nonemptycellsinclusters)).reshape(-1):
            #     clear i11 i22
        #     [i11,i22] = ismember(clusters{index},positions);
            positions3[index]=clusters[nonemptycellsinclusters(index)]
# SpatioTemporalDetectionOFRecurrence.m:150
        positionstoremove=[]
# SpatioTemporalDetectionOFRecurrence.m:153
        for index1 in arange(1,length(positions3)).reshape(-1):
            clear('i1','i2')
            i1,__=ismember(index1,positionstoremove,nargout=2)
# SpatioTemporalDetectionOFRecurrence.m:156
            if i1 == 0:
                for index2 in arange(index1 + 1,length(positions3)).reshape(-1):
                    clear('i11','i22')
                    i11,__=ismember(positions3[index1],positions3[index2],nargout=2)
# SpatioTemporalDetectionOFRecurrence.m:160
                    if sum(i11) > 0:
                        positionstoremove=concat([positionstoremove,index2])
# SpatioTemporalDetectionOFRecurrence.m:162
                        positions3[index1]=unique(union(positions3[index1],positions3[index2]))
# SpatioTemporalDetectionOFRecurrence.m:163
        positions4=cell(1)
# SpatioTemporalDetectionOFRecurrence.m:169
        countpos=1
# SpatioTemporalDetectionOFRecurrence.m:170
        for index in arange(1,length(positions3)).reshape(-1):
            clear('i1','i2')
            i1,__=ismember(index,positionstoremove,nargout=2)
# SpatioTemporalDetectionOFRecurrence.m:173
            if i1 == 0:
                positions4[countpos]=positions3[index]
# SpatioTemporalDetectionOFRecurrence.m:175
                countpos=countpos + 1
# SpatioTemporalDetectionOFRecurrence.m:176
        if length(positions4) == 1:
            if isempty(positions4[1]):
                numberofclusters[countiteration]=0
# SpatioTemporalDetectionOFRecurrence.m:182
            else:
                numberofclusters[countiteration]=length(positions4)
# SpatioTemporalDetectionOFRecurrence.m:184
        else:
            numberofclusters[countiteration]=length(positions4)
# SpatioTemporalDetectionOFRecurrence.m:187
        newclusters=[]
# SpatioTemporalDetectionOFRecurrence.m:191
        if numberofclusters(countiteration) > logical_and(0,countclusters) == 0:
            for checkclusters0 in arange(1,length(positions4)).reshape(-1):
                countclusters=countclusters + 1
# SpatioTemporalDetectionOFRecurrence.m:194
                Totalclusters[countclusters]=positions4[checkclusters0]
# SpatioTemporalDetectionOFRecurrence.m:195
                IntervalCluster[countclusters,arange()]=concat([indexwindow,0])
# SpatioTemporalDetectionOFRecurrence.m:196
        else:
            if numberofclusters(countiteration) > logical_and(0,countclusters) > 0:
                newclusters=copy(positions4)
# SpatioTemporalDetectionOFRecurrence.m:199
                for checkclusters1 in arange(1,length(Totalclusters)).reshape(-1):
                    flagcluster=0
# SpatioTemporalDetectionOFRecurrence.m:202
                    for checkclusters2 in arange(length(newclusters),1,- 1).reshape(-1):
                        i1,i2=ismember(Totalclusters[checkclusters1],newclusters[checkclusters2],nargout=2)
# SpatioTemporalDetectionOFRecurrence.m:204
                        if sum(i1) > 0:
                            Totalclusters[checkclusters1]=union(Totalclusters[checkclusters1],newclusters[checkclusters2])
# SpatioTemporalDetectionOFRecurrence.m:206
                            newclusters[checkclusters2]=[]
# SpatioTemporalDetectionOFRecurrence.m:207
                            flagcluster=1
# SpatioTemporalDetectionOFRecurrence.m:208
                            break
                    if flagcluster == logical_and(0,IntervalCluster(checkclusters1,2)) == 0:
                        IntervalCluster[checkclusters1,2]=indexwindow
# SpatioTemporalDetectionOFRecurrence.m:213
        if numberofclusters(countiteration) == logical_and(0,countclusters) > 0:
            IntervalCluster[IntervalCluster(arange(),2) == 0,2]=indexwindow
# SpatioTemporalDetectionOFRecurrence.m:219
        # check for new clusters in the current window
        for checkclusters3 in arange(1,length(newclusters)).reshape(-1):
            if logical_not(isempty(newclusters[checkclusters3])):
                countclusters=countclusters + 1
# SpatioTemporalDetectionOFRecurrence.m:226
                Totalclusters[countclusters]=newclusters[checkclusters3]
# SpatioTemporalDetectionOFRecurrence.m:227
                IntervalCluster[countclusters,arange()]=concat([indexwindow,0])
# SpatioTemporalDetectionOFRecurrence.m:228
        countiteration=countiteration + 1
# SpatioTemporalDetectionOFRecurrence.m:232
    
    
    if logical_not(isempty(IntervalCluster)):
        IntervalCluster[IntervalCluster(arange(),2) == 0,2]=size(Signals,2)
# SpatioTemporalDetectionOFRecurrence.m:237
    
    ###
    for index1 in arange(1,length(Totalclusters)).reshape(-1):
        poscluster=[]
# SpatioTemporalDetectionOFRecurrence.m:242
        for index2 in arange(index1 + 1,length(Totalclusters)).reshape(-1):
            clear('i11')
            i11,__=ismember(Totalclusters[index2],Totalclusters[index1],nargout=2)
# SpatioTemporalDetectionOFRecurrence.m:245
            if sum(i11) > 0:
                poscluster=concat([poscluster,index2])
# SpatioTemporalDetectionOFRecurrence.m:247
        for index3 in arange(1,length(poscluster)).reshape(-1):
            Totalclusters[index1]=unique(union(Totalclusters[index1],Totalclusters[poscluster(index3)]))
# SpatioTemporalDetectionOFRecurrence.m:251
        if logical_not(isempty(poscluster)):
            IntervalCluster[index1,arange()]=concat([IntervalCluster(index1,1),max(IntervalCluster(concat([index1,poscluster]),2))])
# SpatioTemporalDetectionOFRecurrence.m:254
        for index3 in arange(1,length(poscluster)).reshape(-1):
            Totalclusters[poscluster(index3)]=[]
# SpatioTemporalDetectionOFRecurrence.m:257
    
    
    tempTcluster=copy(Totalclusters)
# SpatioTemporalDetectionOFRecurrence.m:261
    tempinterval=copy(IntervalCluster)
# SpatioTemporalDetectionOFRecurrence.m:262
    clear('Totalclusters','IntervalCluster')
    countT=1
# SpatioTemporalDetectionOFRecurrence.m:264
    for index in arange(1,length(tempTcluster)).reshape(-1):
        if logical_not(isempty(tempTcluster[index])):
            Totalclusters[countT]=tempTcluster[index]
# SpatioTemporalDetectionOFRecurrence.m:267
            IntervalCluster[countT,arange()]=tempinterval(index,arange())
# SpatioTemporalDetectionOFRecurrence.m:268
            countT=countT + 1
# SpatioTemporalDetectionOFRecurrence.m:269
    