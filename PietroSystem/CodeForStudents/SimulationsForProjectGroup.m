%% 2-D sinusoidal wave of time limited duration - 2 sinusoidal waves

%% generate simulation
clear all;close all
% first wave
W1 = 32; H1 = 32;
Xx = [0:W1-1];
Yy = [0:H1-1];
[X,Y] = meshgrid(Xx,Yy);
freqRGB = 1/12;
angleRGB(1) = 2*pi*1/7;
phaseRGB = 0;
imgr = sin(2*pi*freqRGB(1).*((X.*cos(angleRGB(1) )+ Y.*sin(angleRGB(1))))+phaseRGB);

% second wave
W2 = 20; H2 = 20;
Xx2 = [0:W2-1];
Yy2 = [0:H2-1];
[X2,Y2] = meshgrid(Xx2,Yy2);
freqRGB2 = 1/5;
angleRGB2(1) = 2*pi*1/3;
phaseRGB2 = 0;
imgr2 = sin(2*pi*freqRGB2(1).*((X2.*cos(angleRGB2(1) )+ Y2.*sin(angleRGB2(1))))+phaseRGB2);

W = 246;
img1 = randn(W);
[nx,ny,~] = size(img1) ;
[nx1,ny1,~] = size(imgr) ;
[nx2,ny2,~] = size(imgr2) ;
%%place at random location
posx = randsample(nx1:nx-nx1,1) ;
posy = randsample(ny1:ny-ny1,1) ;
posx2 = randsample(nx2:nx-nx2,1) ;
posy2 = randsample(ny2:ny-ny2,1) ;
% img1(posx:posx+nx1-1,posy:posy+ny1-1,:) = imgr ;
T(:,:,1) = img1;
for counter = 1:39
    img1 = randn(W);
    T(:,:,counter) = img1;
end
for counter = 40:139
    if counter <80
    phaseRGB = phaseRGB+pi/180*20;
    imgr = sin(2*pi*freqRGB(1).*((X.*cos(angleRGB(1) )+ Y.*sin(angleRGB(1))))+phaseRGB);
    img1 = randn(W);
    img1(posx:posx+nx1-1,posy:posy+ny1-1,:) = imgr ;
    T(:,:,counter) = img1;
    elseif counter >=80 & counter <100
      phaseRGB = phaseRGB+pi/180*20;
    imgr = sin(2*pi*freqRGB(1).*((X.*cos(angleRGB(1) )+ Y.*sin(angleRGB(1))))+phaseRGB);
    img1 = randn(W);
    img1(posx:posx+nx1-1,posy:posy+ny1-1,:) = imgr ;
    phaseRGB2 = phaseRGB2+pi/180*20;
    imgr2 = sin(2*pi*freqRGB2(1).*((X2.*cos(angleRGB2(1) )+ Y2.*sin(angleRGB2(1))))+phaseRGB2);
    img1(posx2:posx2+nx2-1,posy2:posy2+ny2-1,:) = imgr2 ;
    T(:,:,counter) = img1;
    elseif counter >=100
    img1 = randn(W);
    phaseRGB2 = phaseRGB2+pi/180*20;
    imgr2 = sin(2*pi*freqRGB2(1).*((X2.*cos(angleRGB2(1) )+ Y2.*sin(angleRGB2(1))))+phaseRGB2);
    img1(posx2:posx2+nx2-1,posy2:posy2+ny2-1,:) = imgr2 ;
    T(:,:,counter) = img1;
    end
end

for counter = 140:150
    img1 = randn(W);
    T(:,:,counter) = img1;
end

for counter = 1:size(T,3)
    imshow(squeeze(T(:,:,counter)))
    pause(0.01)
end

%% generate B and L matrix and run function 
B = reshape(T,size(T,1)*size(T,2),size(T,3));

% find connected points - identify regions
L(1,:) = [2 1+W 2+W 0 0 0 0 0];
L(W,:) = [W-1 2*W 2*W-1 0 0 0 0 0];
L(W^2-W+1,:) = [W^2-W+2 W^2-2*W+1 W^2-2*W+2 0 0 0 0 0];
L(W^2,:) = [W^2-1 W^2-W W^2-W-1 0 0 0 0 0];
    
for ind = 2:W-1
    L(ind,:) = [ind-1 ind+1 ind+W-1:ind+W+1 0 0 0]; end

for ind = W^2-W+2:W^2-1
    L(ind,:) = [ind-1 ind+1 ind-W-1:ind-W+1 0 0 0]; end

for ind = W+1:W^2-W
    if rem(ind,W) ~= 0 & rem(ind,W) ~= 1
        L(ind,:) = [ind-1 ind+1 ind-W-1:ind-W+1 ind+W-1:ind+W+1];
    elseif rem(ind,W) == 0
        L(ind,:) = [ind-W-1 ind-W ind-1 ind+W-1 ind+W 0 0 0];
    elseif rem(ind,W) == 1
        L(ind,:) = [ind-W ind-W+1 ind+1 ind+W ind+W+1 0 0 0];
    end
end


%% run function to detect clusters
[Totalclusters, IntervalCluster, numberofclusters] = SpatioTemporalDetectionOFRecurrence(B,L,40,20,0.2,100);

%% visualize clusters

for index = 1:length(Totalclusters)
structure = zeros(246);
structure(Totalclusters{index}) = 1;
imagesc(structure)
title('Cluster 1')
axis square
colormap([1 1 1;0 0 0])
xlabel('point (pixel)'), ylabel('point (pixel)') 
txt = {['first frame: ' num2str(IntervalCluster(index,1))], ['last frame: ' num2str(IntervalCluster(index,2))]};
text(10,30,txt)
set(gcf, 'Position',  [500, 100, 900, 400]) 
pause 
close 
end


