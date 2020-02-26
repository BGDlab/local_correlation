% Code by Aaron Alexander-Bloch, Simon Vandekar and Zhixin Lu
% Alexander-Bloch et al, "Imaging local genetic influences on cortical folding", PNAS, 2020
% rhsulcdepth.csv contains average sulcal depth in fsaverage5 space
% rh.sphere is fsaverage5 spherical surface
% corrfile.csv is a csv file with 4 columns in order from left to right
% 1) index of vertex 1 in fsaverage5
% 2) index of vertex 2 in fsaverage5
% 3) geodesic distance between vertex 1 and vertex 2, calculated on mid-gray surface
% 4) residual of structural covariance between vertex 1 and vertex 2 after regressing out non-linear distance effect
% Note that corrfile includes all and only the pairs of vertices within a 10-mm geodesic
% after exluding noncortical vertices along medial wall: ~27,500 vertices
% Note that according to usage agreements data cannot be shared directly
% please refer to paper for methodological steps to generate corrfile.csv
% depends on Spherical_angles.m available at 

format long 

sulcdepth = csvread('corrfile.csv');
rhogfile=csvread(readfilenamel, 1);
spherefile=read_surf2('/path/to/fsaverage5/surf/rh.sphere');
spherefile.coords = spherefile.coords ./ 100


big90ish=[]
big0ish=[]
bigsd=[]
for i = 1:10242

%get matrix of correlation and change in sulcal depth 
%for node i and nieghbors
index1=find(rhogfile(:,1) == i);
index2=find(rhogfile(:,2) == i);
minirhog=rhogfile(index1,:);
minirhog=[minirhog; rhogfile(index2,[2 1 3 4])];
sulcdepthchange=sulcdepth(minirhog(:,1))-sulcdepth(minirhog(:,2));
sucldepthchangepermm=sulcdepthchange./minirhog(:,3);
minirhog=[minirhog sucldepthchangepermm];
minirhog=[minirhog minirhog(:,4) minirhog(:,4)]; %two identical columns added for consistency with old code
minirhog2=minirhog;

%only proceed if there are 3 neighbors, otherwise set as NaN
if size(minirhog,1) > 2 

nnei=size(minirhog,1);
vertices=spherefile.coords(i,:);

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
%below converts to between -180 and 180 based on angle from vertex
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

%deal with troublesome cases
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

%angle between axis of maximum correlation
%and sulcal axis
minorientcor=mod(orientcor+90,180);
if minorientcor==0 
minorientcor=180;
end
mincor=tomin(minorientcor);
normDeg = mod(orient-orientcor,180);
absDiffDeg = min(180-normDeg, normDeg);

mini90ish=[mini90ish mean(tomin(60:120))];
mini0ish=[mini0ish mean(tomin([1:31 151:180]))];
minisd=[minisd std(tomin)];
end

%take median across nieghbors
%for robustness to numerical issues
big90ish=[big90ish median(mini90ish)];
big0ish=[big0ish median(mini0ish)];
bigsd=[bigsd median(minisd)];

else
big90ish=[big90ish 0/0];
big0ish=[big0ish 0/0];
bigsd=[bigsd 0/0];

end
length(bigorient)
end


Op=(big0ish-big90ish)./bigsd;
