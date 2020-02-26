function [angABC,arcABC, isrealABC] = Spherical_angles(ptA,ptB,ptC)
%This function calculate the angles and edge length of a tringle on the
%unit sphere, assuming that hte center of the sphere is (0,0,0). Results
%are in the unit of radus. 
%code by Zhixin Lu


% normalize the three verticies so that they lie on the unit sphere.
ptA = ptA/norm(ptA);
ptB = ptB/norm(ptB);
ptC = ptC/norm(ptC);

%calculate the length of three arcs, facing the vertices A, B, and C
arcA = atan2(norm(cross(ptB,ptC)),dot(ptB,ptC));
arcB = atan2(norm(cross(ptA,ptC)),dot(ptA,ptC));
arcC = atan2(norm(cross(ptA,ptB)),dot(ptA,ptB));

%calculate the angles at vertices A, B, and C using the spherical rule of 
% cosine.
angA = real(acos((cos(arcA)-cos(arcB)*cos(arcC))/(sin(arcB)*sin(arcC))));
angB = real(acos((cos(arcB)-cos(arcC)*cos(arcA))/(sin(arcC)*sin(arcA))));
angC = real(acos((cos(arcC)-cos(arcA)*cos(arcB))/(sin(arcA)*sin(arcB))));

irA = isreal(acos((cos(arcA)-cos(arcB)*cos(arcC))/(sin(arcB)*sin(arcC))));
irB = isreal(acos((cos(arcB)-cos(arcC)*cos(arcA))/(sin(arcC)*sin(arcA))));
irC = isreal(acos((cos(arcC)-cos(arcA)*cos(arcB))/(sin(arcA)*sin(arcB))));


angABC = [angA,angB,angC];
arcABC = [arcA,arcB,arcB];
isrealABC =[irA,irB,irC];
end

