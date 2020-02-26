cd ~/Documents/tensors/processedfiles/redo/
readfilenamel='~/Documents/tensors/solarfiles/rhoPhcpR25_distres.csv'
bigfilenamel='~/Documents/tensors/HCPrhoP10mmrightDistres.csv';
writefilenamel='rhoP_neighbors_rightv3'
orientwritefilenamel=[writefilenamel '_orient.csv']
meanwritefilenamel=[writefilenamel '_mean.csv']
maxwritefilenamel=[writefilenamel '_max.csv']
aniwritefilenamel=[writefilenamel '_axialVStransverse.csv']
meandiffwritefilenamel=[writefilenamel '_meandif.csv']
ang0writefilenamel=[writefilenamel '_ang0.csv']
ang90writefilenamel=[writefilenamel '_ang90.csv']
stdwritefilenamel=[writefilenamel '_std.csv']
ang0ishwritefilenamel=[writefilenamel '_ang0ish30.csv']
ang90ishwritefilenamel=[writefilenamel '_ang90ish30.csv']


%_meandif.csv mean difference across all angles
%_meanani max-90degreesfrommax
%_std std of the values

sulcdepthlh = csvread('rhsulcdepth.csv');
format long %so there are enough decimal places
rhogfile=csvread(readfilenamel, 1);
bigrhogfile=csvread(bigfilenamel, 1);

newrhogfile=[];
for i=1:size(rhogfile,1)

index1=find(bigrhogfile(:,1) == rhogfile(i,1) & bigrhogfile(:,2) == rhogfile(i,2));
index2=find(bigrhogfile(:,2) == rhogfile(i,1) & bigrhogfile(:,1) == rhogfile(i,2));
index=[index1 index2];
newrhogfile=[newrhogfile; [rhogfile(i,1:2) bigrhogfile(index(1),3) rhogfile(i,3)]];
	end
rhogfile=newrhogfile;
sulcdepth=sulcdepthlh;
spherefile=read_surf2('/Applications/freesurfer/subjects/fsaverage5/surf/rh.sphere');
spherefile.coords = spherefile.coords ./ 100
%trait1, trait2, dist, rhog

%getdistfrombigfile


bigmean=[]
biganisotropy=[]
bigorient=[]
bigmax=[]
big90=[]
big0=[]
big90ish=[]
big0ish=[]
bigmeandif=[]
bigsd=[]
for i = 1:10242

%get matrix of correlation and change in sulcal depth 
%for node i and nieghbors
%order of columns is arbitrary/redundant for consistency with old code
index1=find(rhogfile(:,1) == i);
index2=find(rhogfile(:,2) == i);
minirhog=rhogfile(index1,:);
minirhog=[minirhog; rhogfile(index2,[2 1 3 4])];
sulcdepthchange=sulcdepth(minirhog(:,1))-sulcdepth(minirhog(:,2));
sucldepthchangepermm=sulcdepthchange./minirhog(:,3);
minirhog=[minirhog sucldepthchangepermm];
minirhog=[minirhog minirhog(:,4) minirhog(:,4)];
minirhog2=minirhog;
%only proceed if there are 3 neighbors, otherwas NaN
if size(minirhog,1) > 2 

nnei=size(minirhog,1);
vertices=spherefile.coords(i,:);
miniorient=[];
minianisotropy=[];
minimean=[];
minimax=[];
minimeandiff=[];
mini0=[];
mini90=[];
mini0ish=[];
mini90ish=[];
minisd=[];

%iterate through all neighbors as reference point, 
%to ensure robustness to numerical problems
for whichr=1:nnei
	minirhog=minirhog2;
	bigangx = [];
bigangz = [];
bigir=[];
ni = minirhog(whichr,2);
reference=spherefile.coords(ni,:);

%iterate through all neighbors as final point in triangle
for j = 1:nnei
ni = minirhog(j,2);
neighbors=spherefile.coords(ni,:);
angles=rad2deg(Spherical_angles(vertices,neighbors,reference));
angz=angles(1);
[~, ~, ir]=Spherical_angles(vertices,neighbors,reference);
bigangz = [bigangz angz];
bigir=[bigir all(ir)];
end

%bigangz gives angle in degrees but it will always be between 0 and 180
%below makes it go between -180 and 180 based on angle from vertex
%that is approx perpindicular to reference 
[~, perp]=min(abs(90-bigangz));
perpy=90-bigangz(perp);
ni = minirhog(perp,2);
reference=spherefile.coords(ni,:);
for j = 1:nnei
ni = minirhog(j,2);
neighbors=spherefile.coords(ni,:);
angles=rad2deg(Spherical_angles(vertices,neighbors,reference));
angz=angles(1);
bigangx = [bigangx angz];
end

indexper=find(((bigangx > (90-abs(perpy))) & (bigangz > (90+perpy))) | (bigangx > (90+abs(perpy))));

bigangz(indexper)= bigangz(indexper) .* -1;


bigangz(whichr)=0;
%bigangz=bigangz*-1

%code below is to deal with troublesome cases
%where multiple points are same angle
if any(bigangz == -180)
indexer=find(bigangz == -180);
bigangz(indexer)=180;
end
thingclose=abs(bigangz - bigangz');
thingclose(find(eye(size(thingclose))))=5;

thresh=1;
if sum(sum(thingclose < thresh)) >0
repvals=find(sum(thingclose < thresh));
while length(repvals)>0
repvalsi=repvals(1);
findrep2=find(thingclose(:,repvalsi) <thresh);
findrep=[findrep2; repvalsi];
minirhog(findrep,5)=mean(minirhog(findrep,5));
minirhog(findrep,7)=mean(minirhog(findrep,7));
minirhog(findrep2,:)=[];
bigangz(findrep2)=[];
thingclose(findrep2,:)=[];
thingclose(:,findrep2)=[];
repvals=find(sum(thingclose < thresh));
end
end

%interpolate change in sulcal depth across 360 by degree 
a=(-179:180);
tensorrot = [];
for kk=1:360
normDeg = mod(a(kk)-bigangz,360);
normDegx = abs(a(kk)-bigangz);

absDiffDeg = min(360-normDeg, normDeg);
absDiffDegx = min(360-normDegx, normDegx);

[val1, indexthing1]=min(360-normDeg);
[val2 ,indexthing2]=min(normDeg);
if any(val1==0 | val2 == 0)
[~,indexmin]=min([val1 val2]);
listindex=[indexthing1 indexthing2];
indexthing=listindex(indexmin);
tensor = minirhog(indexthing,5);
else
tensor = ((1/absDiffDeg(indexthing1)*(minirhog(indexthing1,5))) + (1/absDiffDeg(indexthing2)*(minirhog(indexthing2,5))))/((1/absDiffDeg(indexthing1)) + (1/absDiffDeg(indexthing2)));
end
tensorrot = [tensorrot tensor];
end

%find axis of minimize change in sulcal depth = sulcal axis
tomin=[];
for kk=1:180
tomin=[tomin abs(tensorrot(kk))+abs(tensorrot(kk+180))];
end
[~, orient]=min(tomin);

%for troubleshooting plots below
%orientmin=tomin;
%orientrot=tensorrot;

%interpolate correlation across 360 by degree
tensorrot = [];
for kk=1:360
normDeg = mod(a(kk)-bigangz,360);
absDiffDeg = min(360-normDeg, normDeg);
[val1, indexthing1]=min(360-normDeg);
[val2 ,indexthing2]=min(normDeg);
if any(val1==0 | val2 == 0)
[~,indexmin]=min([val1 val2]);
listindex=[indexthing1 indexthing2];
indexthing=listindex(indexmin);
tensor = minirhog(indexthing,7);
else
tensor = ((1/absDiffDeg(indexthing1)*(minirhog(indexthing1,7))) + (1/absDiffDeg(indexthing2)*(minirhog(indexthing2,7))))/((1/absDiffDeg(indexthing1)) + (1/absDiffDeg(indexthing2)));
end
tensorrot = [tensorrot tensor];
end

%find axis of maximum correlation
tomin=[];
for kk=1:180
tomin=[tomin (tensorrot(kk)+tensorrot(kk+180))];
end
[maxcor, orientcor]=max(tomin);

%troubleshooting plots
%cormin=tomin;
%corrot=tensorrot;
%figure
%plot([corrot; orientrot]')
%title(['360' ', whichr=' num2str(whichr)])
%figure
%plot([cormin; orientmin]')
%title(['180' ', whichr=' num2str(whichr)])

%angle between axis of maximum correlation
%and sulcal axis
minorientcor=mod(orientcor+90,180);
if minorientcor==0 
minorientcor=180;
end
mincor=tomin(minorientcor);
normDeg = mod(orient-orientcor,180);
absDiffDeg = min(180-normDeg, normDeg);
miniorient=[miniorient absDiffDeg];

meandiff=[];
for tomini=1:180
tominnoi=tomin;
tominnoi(tomini)=[];
meandiff=[meandiff mean(abs(tomin(tomini)-tominnoi))];
end

avemeandiff=mean(meandiff);
axVtr=tomin(90)-tomin(180);


minianisotropy=[minianisotropy axVtr];
minimean=[minimean mean(tomin)];
minimax=[minimax maxcor];
minimeandiff=[minimeandiff avemeandiff];
mini0=[mini0 tomin(180)];
mini90=[mini90 tomin(90)];

mini90ish=[mini90ish mean(tomin(60:120))];
mini0ish=[mini0ish mean(tomin([1:31 151:180]))];
minisd=[minisd std(tomin)];
end

%take median across nieghbors
%for robustness to numerical issues
bigorient=[bigorient median(miniorient)];
biganisotropy=[biganisotropy median(minianisotropy)];
bigmean=[bigmean median(minimean)];
bigmax=[bigmax median(minimax)];
big90=[big90 median(mini90)];
big0=[big0 median(mini0)];
big90ish=[big90ish median(mini90ish)];
big0ish=[big0ish median(mini0ish)];
bigsd=[bigsd median(minisd)];

bigmeandif=[bigmeandif median(minimeandiff)];
else
bigorient=[bigorient 0/0];
biganisotropy=[biganisotropy 0/0];
bigmean=[bigmean 0/0];
bigmax=[bigmax 0/0];
big90=[big90 0/0];
big0=[big0 0/0];
big90ish=[big90ish 0/0];
big0ish=[big0ish 0/0];
bigmeandif=[bigmeandif 0/0];
bigsd=[bigsd 0/0];

end
length(bigorient)
%max(max(abs(miniorient - miniorient')))
end

csvwrite(orientwritefilenamel,bigorient)
csvwrite(meanwritefilenamel,bigmean)
csvwrite(maxwritefilenamel,bigmax)
csvwrite(aniwritefilenamel,biganisotropy)
csvwrite(meandiffwritefilenamel, bigmeandif)
csvwrite(ang0writefilenamel, big0)
csvwrite(ang90writefilenamel,big90)
csvwrite(ang0ishwritefilenamel, big0ish)
csvwrite(ang90ishwritefilenamel,big90ish)
csvwrite(stdwritefilenamel,bigsd)




%0-90 where 90 is perp and 0 is par