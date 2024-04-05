function DomStudioOrientation(~)
% @ 2015 by Andrea Bistacchi, distributed under the GNU AGPL v3.0 license.
%
% Function used for the analysis of orientation data imported as CSV files with Dip / Dip Azimuth
%
% Last update by Stefano Casiraghi 8/2/2024


% initialize
clear all; close all; clc; clearvars;
rad = pi/180;
deg = 180/pi;

% file overwrite warning
disp(' ')
disp('WARNING: this script overwrites output files without asking for confirmaton.')

% load CSV file with facets data
[file, path] = uigetfile('*.csv');
filename = [path file];
[path,file,extension] = fileparts(filename);
inDataTable = readtable(filename);

% copy data to variables
DipAzimuth = inDataTable.DipDir_;
Dip = inDataTable.Dip;

% read data
Ndata = size(inDataTable,1);  % Ndata = number of data points
disp(' ')
disp(['Loaded data number: ' num2str(Ndata)])

% check
if max(Dip)>90, error('DIP ERROR'); end
if max(DipAzimuth)>360, error('AZIMUTH ERROR'); end

% strike (used by rose diagram only)
Strike = (DipAzimuth<90).*(DipAzimuth+270) + (DipAzimuth>=90).*(DipAzimuth-90);
SymmetricStrike = [Strike; (Strike<=180).*(Strike+180)+(Strike>180).*(Strike-180)];

% plunge/trend of poles to planes
Plunge = 90-Dip;
Trend = (DipAzimuth<=180).*(DipAzimuth+180)+(DipAzimuth>180).*(DipAzimuth-180);

%direction cosines of data vectors
Lpole = cos(Plunge*rad).*cos(Trend*rad);
Mpole = cos(Plunge*rad).*sin(Trend*rad);
Npole = -sin(Plunge*rad);  %we use lower hemisphere

% symmetric unit vectors -> "reflected unit vectors"
Lpole_neg = -Lpole; 
Mpole_neg = -Mpole;
Npole_neg = -Npole;

% union of unit vectors -> arrays are twice the length of input data
Lpole = [Lpole;Lpole_neg];
Mpole = [Mpole;Mpole_neg];
Npole = [Npole;Npole_neg];

% K-medoid clustering
disp(' ')
disp('Number of classes for K-medoid clustering [1]:')
nClass = input('>> ');
if nClass < 1, nClass = 1; end
nClass = round(nClass)*2;
idClass = kmedoids([Lpole Mpole Npole],nClass);
countClass = zeros(1,nClass);
for i = 1:nClass
    countClass(i) = sum(idClass==i);
end
countClassPercent = countClass./Ndata.*100;
disp(' ')
disp('K-medoid clustering done')

% Fisher's mean and K
for i = 1:nClass
    disp(["______________________" num2str(i)])
    % collect data for class
    L = Lpole(idClass==i);
    M = Mpole(idClass==i);
    N = Npole(idClass==i);
    sumL = sum(L);
    sumM = sum(M);
    sumN = sum(N);
    R = sqrt(sumL^2+sumM^2+sumN^2);
    sumLR = sumL/R;
    sumMR = sumM/R;
    sumNR = sumN/R;
    meanP = asin(-sumNR)*deg;
    meanT = atan2(sumMR,sumLR)*deg;
    meanT = meanT+(meanT<0)*360;
    fisherK(i) = (countClass(i)-1)/(countClass(i)-R);  % this is OK for n_data > 16
    confC(i) = acos( 1 - (countClass(i)-R)/R * ((1/0.01)^(1/(countClass(i)-1)) - 1) )*deg;
    spherAp(i) = asin(sqrt(2*(1-1/countClass(i))/fisherK(i)))*deg;
    meanDip(i) = 90-meanP;
    meanDir(i) = (meanT<=180).*(meanT+180)+(meanT>180).*(meanT-180);
    disp([num2str(meanDip(i)) " - " num2str(meanDir(i)) " - " num2str(fisherK(i))])
    clear L M N aL aM aN
end
disp(' ')
disp('Fisher means and Ks done')

% remove mean vectors pointing upwards
keep_class = find(meanDip<=90);
fisherK = fisherK(meanDip<=90);
confC = confC(meanDip<=90);
spherAp = spherAp(meanDip<=90);
meanDir = meanDir(meanDip<=90);
meanDip = meanDip(meanDip<=90); % keep this as the last one
nClass = nClass/2;


% plunge/trend of node vectors for contouring grid by direct calculation
n = 1;
p(n) = 90;
t(n) = 0;
for i = 1:10
    m = 6*i;
    radius = i/10;
    DeltaPhi = 360/m;
    for j = 1:m
        phi = j*DeltaPhi;
        n = n+1;
        t(n) = phi;
        theta = 2 * asin(radius/sqrt(2));
        p(n) = 90 - (theta*deg);
    end
end

% direction cosines (l,m,n) of node vectors for contouring grid
for i = 1:331
    Ln(i) = cos(p(i)*rad) * cos(t(i)*rad);
    Mn(i) = cos(p(i)*rad) * sin(t(i)*rad);
    Nn(i) = -sin(p(i)*rad);
end

% (x,y) coordinates of nodes for contouring grid from plunge/trend
for n = 1:331
    theta = (90-p(n))*rad;
    d = sqrt(2) * sin(theta/2);
    x(n) = d .* sin(t(n)*rad);
    y(n) = d .* cos(t(n)*rad);
end
disp(' ')
disp('grid nodes done')

% define counting cone for contouring grid
% key = 2.0;
% cone = (acos(Ndata/(Ndata+key^2)))*deg;
cone = 8.1096144559941786958201832872484; %%% http://en.wikipedia.org/wiki/Solid_angle

%count point density at all 331 nodes
for i = 1:331
    z(i) = 0;
    %dot product of node vectors and data vectors
    for j = 1:Ndata
        % theta = (acos(Ln(i)*Ld(j) + Mn(i)*Md(j) + Nn(i)*Ld(j)))*deg;
        theta = (acos(Ln(i)*Lpole(j) + Mn(i)*Mpole(j) + Nn(i)*Npole(j)))*deg;
        if theta <= cone
            z(i) = z(i)+1;
        end
    end
end
disp(' ')
disp('point density done') %%%

%Zmax = 0;
for j = 1:331
    %convert counts to percent
    z(j) = (z(j)/Ndata)*100;
    % if z(j) > Zmax
    %     Zmax = z(j);
    % end
end
disp(' ')
disp('density percent done')

% contouring Kalsbeek
% set parameters:
% number of x and y nodes for grid
% make these smaller if you get an error about excceding the max array size
nx = 50;
ny = nx;
% number of contours:
ci = 5;
e=x;
n=y;
hh=z;
%Grid the data linearly:
estep = abs((min(e)-max(e))/nx);
nstep = abs((min(n)-max(n))/ny);
ee = min(e):estep:max(e);
nn = min(n):nstep:max(n);
[xi,yi,zi] = griddata(e,n,hh,ee,nn');
% v = floor(min(hh)):ci:max(hh);
v = floor(linspace(min(hh),max(hh),ci));
disp(' ')
disp('gridding done') %%%

% STEREOPLOT
figure('WindowStyle','docked');
cmap =[ones(1,10)' linspace(1,0,10)' linspace(1,0,10)'];

%draw primitive circle
radius = 1;
theta = linspace(0,2*pi,180);
x = radius * cos(theta);
y = radius * sin(theta);
plot(x,y,'k-')
hold on
axis equal
axis off
imagesc(ee',nn',zi)
colormap(cmap)
plot(x,y,'k-')
plot(0,0,'k+')
plot(radius,0,'k+')
plot(0,radius,'k+')
plot(-radius,0,'k+')
plot(0,-radius,'k+')
plot(0,0,'k+')

%plot and label contours
[c,h] = contour(xi,yi,zi,v,'-k');
clabel(c,h);

% plot cluster means
for i = 1:nClass
    theta = (meanDip(i))*rad;
    d = sqrt(2) * sin(theta/2);
    xMean = d.* -sin((meanDir(i))*rad);
    yMean = d.* -cos((meanDir(i))*rad);
    plot(xMean,yMean,'kd','MarkerSize',6,'linewidth',2,'MarkerFaceColor','w');
end
t = title({['Poles in class = ' num2str(countClass)];...
    ['percent = ' num2str(countClassPercent)];...
    ['Mean Dip = ' num2str(meanDip)];...
    ['Mean Dir = ' num2str(meanDir)];...
    ['Fisher K = ' num2str(fisherK)];...
    ['99% confidence cone apical angle = ' num2str(confC)];...
    ['68.26% variability spherical aperture = ' num2str(spherAp)];...
    'Concentrations % of total per 1% area';...
    ['Maximum concentration = ' num2str(max(z))]});
set(t, 'horizontalAlignment', 'l')

% save as jpeg
print('-djpeg',[filename '_contour_Km.jpg'])
disp(' ')
disp('stereoplot done')

% ROSE DIAGRAM
figure('WindowStyle','docked');
% geologic rose plot inpired by earth_rose.m by Cameron Sparr
% in turn inspired by wind_rose.m
% The code to change the axis labels from a mathematical system to
% N,S,E,W were written by Jordan Rosenthal in a forum post:
% http://www.mathworks.com/matlabcentral/newsreader/author/4601
D = mod(90 - SymmetricStrike, 360)*pi/180;
rose(D, 36);   % rose with bins at 10° increment
hHiddenText = findall(gca,'type','text');
Angles = 0 : 30 : 330;
hObjToDelete = zeros( length(Angles)-4, 1 );
k = 0;
for ang = Angles
    hObj = findall(hHiddenText,'string',num2str(ang));
    switch ang
        case 0
            set(hObj,'string','E','HorizontalAlignment','Left');
        case 90
            set(hObj,'string','N','VerticalAlignment','Bottom');
        case 180
            set(hObj,'string','W','HorizontalAlignment','Right');
        case 270
            set(hObj,'string','S','VerticalAlignment','Top');
        otherwise
            k = k + 1;
            hObjToDelete(k) = hObj;
    end
end
delete(hObjToDelete(hObjToDelete~=0));
t = title({['Poles in class = ' num2str(countClass)];...
    ['percent = ' num2str(countClassPercent)];...
    ['Mean Dip = ' num2str(meanDip)];...
    ['Mean Dir = ' num2str(meanDir)];...
    ['Fisher K = ' num2str(fisherK)];...
    ['99% confidence cone apical angle = ' num2str(confC)];...
    ['68.26% variability spherical aperture = ' num2str(spherAp)];...
    'Concentrations % of total per 1% area';...
    ['Maximum concentration = ' num2str(max(z))]});
set(t, 'horizontalAlignment', 'left')

% save as jpeg
print('-djpeg',[filename '_rose_Km.jpg'])
disp(' ')
disp('rose diagram done')

% button = 1; xx = []; yy=[];
% while button ==1, [x,y,button] = ginput(1); xx = [xx x]; yy = [yy y]; end

% save with class
out_dip = [Dip; Dip];
out_azi = [DipAzimuth; DipAzimuth];
out_data = [out_dip out_azi idClass];
out_data = out_data(find(ismember(out_data(:,3),keep_class)),:);
csvwrite([filename '_classified.csv'],out_data)
disp(' ')
disp('output done')

% Results table
% Class_number = cell2table(keep_class);
% Mean_Dip_Dir = num2str(meanDir);
% Mean_Dip = num2str(meanDip);
% Fisher_K = num2str(fisherK);
% Results_table = table(Class_number, Mean_Dip_Dir, Mean_Dip, Fisher_K);
% writetable(Results_table, [file '_results_table']);
