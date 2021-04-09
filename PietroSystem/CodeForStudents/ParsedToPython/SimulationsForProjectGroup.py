# Generated with SMOP  0.41-beta
from libsmop import *
# SimulationsForProjectGroup.m

    ## 2-D sinusoidal wave of time limited duration - 2 sinusoidal waves
    
    ## generate simulation
    clear('all')
    close_('all')
    # first wave
    W1=32
# SimulationsForProjectGroup.m:6
    H1=32
# SimulationsForProjectGroup.m:6
    Xx=concat([arange(0,W1 - 1)])
# SimulationsForProjectGroup.m:7
    Yy=concat([arange(0,H1 - 1)])
# SimulationsForProjectGroup.m:8
    X,Y=meshgrid(Xx,Yy,nargout=2)
# SimulationsForProjectGroup.m:9
    freqRGB=1 / 12
# SimulationsForProjectGroup.m:10
    angleRGB[1]=dot(dot(2,pi),1) / 7
# SimulationsForProjectGroup.m:11
    phaseRGB=0
# SimulationsForProjectGroup.m:12
    imgr=sin(multiply(dot(dot(2,pi),freqRGB(1)),((multiply(X,cos(angleRGB(1))) + multiply(Y,sin(angleRGB(1)))))) + phaseRGB)
# SimulationsForProjectGroup.m:13
    # second wave
    W2=20
# SimulationsForProjectGroup.m:16
    H2=20
# SimulationsForProjectGroup.m:16
    Xx2=concat([arange(0,W2 - 1)])
# SimulationsForProjectGroup.m:17
    Yy2=concat([arange(0,H2 - 1)])
# SimulationsForProjectGroup.m:18
    X2,Y2=meshgrid(Xx2,Yy2,nargout=2)
# SimulationsForProjectGroup.m:19
    freqRGB2=1 / 5
# SimulationsForProjectGroup.m:20
    angleRGB2[1]=dot(dot(2,pi),1) / 3
# SimulationsForProjectGroup.m:21
    phaseRGB2=0
# SimulationsForProjectGroup.m:22
    imgr2=sin(multiply(dot(dot(2,pi),freqRGB2(1)),((multiply(X2,cos(angleRGB2(1))) + multiply(Y2,sin(angleRGB2(1)))))) + phaseRGB2)
# SimulationsForProjectGroup.m:23
    W=246
# SimulationsForProjectGroup.m:25
    img1=randn(W)
# SimulationsForProjectGroup.m:26
    nx,ny,__=size(img1,nargout=3)
# SimulationsForProjectGroup.m:27
    nx1,ny1,__=size(imgr,nargout=3)
# SimulationsForProjectGroup.m:28
    nx2,ny2,__=size(imgr2,nargout=3)
# SimulationsForProjectGroup.m:29
    ##place at random location
    posx=randsample(arange(nx1,nx - nx1),1)
# SimulationsForProjectGroup.m:31
    posy=randsample(arange(ny1,ny - ny1),1)
# SimulationsForProjectGroup.m:32
    posx2=randsample(arange(nx2,nx - nx2),1)
# SimulationsForProjectGroup.m:33
    posy2=randsample(arange(ny2,ny - ny2),1)
# SimulationsForProjectGroup.m:34
    # img1(posx:posx+nx1-1,posy:posy+ny1-1,:) = imgr ;
    T[arange(),arange(),1]=img1
# SimulationsForProjectGroup.m:36
    for counter in arange(1,39).reshape(-1):
        img1=randn(W)
# SimulationsForProjectGroup.m:38
        T[arange(),arange(),counter]=img1
# SimulationsForProjectGroup.m:39
    
    for counter in arange(40,139).reshape(-1):
        if counter < 80:
            phaseRGB=phaseRGB + dot(pi / 180,20)
# SimulationsForProjectGroup.m:43
            imgr=sin(multiply(dot(dot(2,pi),freqRGB(1)),((multiply(X,cos(angleRGB(1))) + multiply(Y,sin(angleRGB(1)))))) + phaseRGB)
# SimulationsForProjectGroup.m:44
            img1=randn(W)
# SimulationsForProjectGroup.m:45
            img1[arange(posx,posx + nx1 - 1),arange(posy,posy + ny1 - 1),arange()]=imgr
# SimulationsForProjectGroup.m:46
            T[arange(),arange(),counter]=img1
# SimulationsForProjectGroup.m:47
        else:
            if counter >= logical_and(80,counter) < 100:
                phaseRGB=phaseRGB + dot(pi / 180,20)
# SimulationsForProjectGroup.m:49
                imgr=sin(multiply(dot(dot(2,pi),freqRGB(1)),((multiply(X,cos(angleRGB(1))) + multiply(Y,sin(angleRGB(1)))))) + phaseRGB)
# SimulationsForProjectGroup.m:50
                img1=randn(W)
# SimulationsForProjectGroup.m:51
                img1[arange(posx,posx + nx1 - 1),arange(posy,posy + ny1 - 1),arange()]=imgr
# SimulationsForProjectGroup.m:52
                phaseRGB2=phaseRGB2 + dot(pi / 180,20)
# SimulationsForProjectGroup.m:53
                imgr2=sin(multiply(dot(dot(2,pi),freqRGB2(1)),((multiply(X2,cos(angleRGB2(1))) + multiply(Y2,sin(angleRGB2(1)))))) + phaseRGB2)
# SimulationsForProjectGroup.m:54
                img1[arange(posx2,posx2 + nx2 - 1),arange(posy2,posy2 + ny2 - 1),arange()]=imgr2
# SimulationsForProjectGroup.m:55
                T[arange(),arange(),counter]=img1
# SimulationsForProjectGroup.m:56
            else:
                if counter >= 100:
                    img1=randn(W)
# SimulationsForProjectGroup.m:58
                    phaseRGB2=phaseRGB2 + dot(pi / 180,20)
# SimulationsForProjectGroup.m:59
                    imgr2=sin(multiply(dot(dot(2,pi),freqRGB2(1)),((multiply(X2,cos(angleRGB2(1))) + multiply(Y2,sin(angleRGB2(1)))))) + phaseRGB2)
# SimulationsForProjectGroup.m:60
                    img1[arange(posx2,posx2 + nx2 - 1),arange(posy2,posy2 + ny2 - 1),arange()]=imgr2
# SimulationsForProjectGroup.m:61
                    T[arange(),arange(),counter]=img1
# SimulationsForProjectGroup.m:62
    
    for counter in arange(140,150).reshape(-1):
        img1=randn(W)
# SimulationsForProjectGroup.m:67
        T[arange(),arange(),counter]=img1
# SimulationsForProjectGroup.m:68
    
    for counter in arange(1,size(T,3)).reshape(-1):
        imshow(squeeze(T(arange(),arange(),counter)))
        pause(0.01)
    
    ## generate B and L matrix and run function
    B=reshape(T,dot(size(T,1),size(T,2)),size(T,3))
# SimulationsForProjectGroup.m:77
    # find connected points - identify regions
    L[1,arange()]=concat([2,1 + W,2 + W,0,0,0,0,0])
# SimulationsForProjectGroup.m:80
    L[W,arange()]=concat([W - 1,dot(2,W),dot(2,W) - 1,0,0,0,0,0])
# SimulationsForProjectGroup.m:81
    L[W ** 2 - W + 1,arange()]=concat([W ** 2 - W + 2,W ** 2 - dot(2,W) + 1,W ** 2 - dot(2,W) + 2,0,0,0,0,0])
# SimulationsForProjectGroup.m:82
    L[W ** 2,arange()]=concat([W ** 2 - 1,W ** 2 - W,W ** 2 - W - 1,0,0,0,0,0])
# SimulationsForProjectGroup.m:83
    for ind in arange(2,W - 1).reshape(-1):
        L[ind,arange()]=concat([ind - 1,ind + 1,arange(ind + W - 1,ind + W + 1),0,0,0])
# SimulationsForProjectGroup.m:86
    
    for ind in arange(W ** 2 - W + 2,W ** 2 - 1).reshape(-1):
        L[ind,arange()]=concat([ind - 1,ind + 1,arange(ind - W - 1,ind - W + 1),0,0,0])
# SimulationsForProjectGroup.m:89
    
    for ind in arange(W + 1,W ** 2 - W).reshape(-1):
        if rem(ind,W) != logical_and(0,rem(ind,W)) != 1:
            L[ind,arange()]=concat([ind - 1,ind + 1,arange(ind - W - 1,ind - W + 1),arange(ind + W - 1,ind + W + 1)])
# SimulationsForProjectGroup.m:93
        else:
            if rem(ind,W) == 0:
                L[ind,arange()]=concat([ind - W - 1,ind - W,ind - 1,ind + W - 1,ind + W,0,0,0])
# SimulationsForProjectGroup.m:95
            else:
                if rem(ind,W) == 1:
                    L[ind,arange()]=concat([ind - W,ind - W + 1,ind + 1,ind + W,ind + W + 1,0,0,0])
# SimulationsForProjectGroup.m:97
    
    ## run function to detect clusters
    Totalclusters,IntervalCluster,numberofclusters=SpatioTemporalDetectionOFRecurrence(B,L,40,20,0.2,100,nargout=3)
# SimulationsForProjectGroup.m:103
    ## visualize clusters
    
    for index in arange(1,length(Totalclusters)).reshape(-1):
        structure=zeros(246)
# SimulationsForProjectGroup.m:108
        structure[Totalclusters[index]]=1
# SimulationsForProjectGroup.m:109
        imagesc(structure)
        title('Cluster 1')
        axis('square')
        colormap(concat([[1,1,1],[0,0,0]]))
        xlabel('point (pixel)')
        ylabel('point (pixel)')
        txt=cellarray([concat(['first frame: ',num2str(IntervalCluster(index,1))]),concat(['last frame: ',num2str(IntervalCluster(index,2))])])
# SimulationsForProjectGroup.m:115
        text(10,30,txt)
        set(gcf,'Position',concat([500,100,900,400]))
        pause
        close_
    