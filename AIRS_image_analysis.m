close all
clear all
fname='AIRS_read/AIRS.2018.12.24.121.L1B.SUBIBRAD.v5.0.23.0.G18359100630.hdf'; % read data

%%%%specific variables are read in from the AIRS file
Latitude=hdfread(fname,'/L1B_AIRS_Science/Geolocation Fields/Latitude', 'Index',{[1 1],[1 1],[135 90]});
Longitude=hdfread(fname,'/L1B_AIRS_Science/Geolocation Fields/Longitude', 'Index',{[1 1],[1 1],[135 90]});
Channel_Subset=hdfread(fname, '/Channel_Subset','Index',{[1],[1],[130]}); %1-130
radiances=hdfread(fname, 'L1B_AIRS_Science/Data Fields/radiances', 'Index',{[1 1 1],[1 1 1],[135 90 130]});
CH75=radiances(:,:,75);
%%Max,Min
Lmax=max(max(Longitude));
Lmin=min(min(Longitude));
Lamax=max(max(Latitude));
Lamin=min(min(Latitude));
%%
%%Background Noise Removal
x=1:90;
for i=1:135
    P=polyfit(x, CH75(i,:),4);
    p_test=P(1).*x.^4+P(2).*x.^3+P(3).*x.^2+P(4).*x+P(5);
    rad75_nb(i,:)=CH75(i,:)-p_test; 
end

figure(100)
subplot(2,1,1)
plot([1:90],p_test,'-r',[1:90],CH75(135,:),'-g')
legend('Polynomial Fitting Data','AIRS Raw data')
subplot(2,1,2)
plot([1:90],rad75_nb(135,:),'-b')
legend('Background noise eliminated data')

figure(1)
pcolor(CH75)
title('With Noise Channel 75 data')
shading flat
colorbar
figure(2)
pcolor(rad75_nb)
title('Without Noise Channel 75 data')
shading flat
colorbar
%%
% D=distance(Latitude,Longitude);

%%Observe data
load coastlines
[latcells, loncells] = polysplit(coastlat, coastlon);
numel(latcells);
% %%%%plot on global projection
rad75_nb_db=double(rad75_nb);
clear Longitude_new;
clear rad75_new;
clear Latitude_new;
start_La=35;
end_La=70;
start_Lo=35;
end_Lo=70;
clear rad;
% J=1;
for i=1:end_La-start_La
rad(i,:)=rad75_nb_db(start_La-1+i,start_Lo+1:end_Lo);%35th row column to 70th row column 
Lat(i,:)=Latitude(start_La-1+i,start_Lo+1:end_Lo);
Long(i,:)=Longitude(start_La-1+i,start_Lo+1:end_Lo);
end
clear Lat_new;
clear rad75L_new 
rad75L_new=zeros(35,35);
for i=1:35
Long_new(i,:)=linspace(Long(i,35),Long(i,1),35);
Lat1_new=(linspace(Lat(35,i),Lat(1,i),35));
Lat_new(:,i)=Lat1_new;
rad75_new(i,:)=spline(Long(i,:),rad(i,:),Long_new(i,:));
end

for i=1:35
Rad_new(:,i)=spline(Lat(:,i),rad75_new(:,i),Lat_new(:,i));
end
radd=Rad_new;
figure(200)
plot([1:35],rad,'-r',[1:35],rad75_new,'-g')

figure(31)
subplot(2,1,1)
pcolor(rad)
title('Channel 75 Raw original data')
shading flat
colorbar
subplot(2,1,2)
pcolor(radd)
title('Channel 75 data after interpolation')
shading flat
colorbar

d=zeros(35,35);
D=zeros(35,35);
for j=1:34
for i=1:34
    d(j,i+1)=deg2km(distance('gc',Long(j,i),Lat_new(j,i),Long(j,i+1),Lat_new(j,i+1)))+d(j,i);
    D(j+1,i)=deg2km(distance('gc',Long(j,i),Lat_new(j,i),Long(j+1,i),Lat_new(j+1,i)))+D(j,i);
end
end
Fs1=1/(1000*abs(d(1,2)-d(1,1)));
Fs2=1/(1000*abs(D(2,1)-D(1,1)));
Lx=linspace(-Fs1,Fs1,length(radd(1,:)));
Ly=linspace(-Fs2,Fs2,length(radd(:,1)));
[Lx,Ly]=meshgrid(Lx,Ly);
RadFF=abs(fftshift(fft2(radd)));
figure(9)
mesh(Lx,Ly,RadFF);
shading flat

