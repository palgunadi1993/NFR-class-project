function twophaseflow
% Open mesh

% domain size is determined by number and size of gridcells
% ========================= Input Data ==================================
% --------------------    Grid Dimensions   -------------------------
% load fine_mesh.mat
load complex_mesh.mat

num_nodes = msh.nbNod;
matrix_pos1 = msh.POS(:,1:2);
matrix_nodes1 = msh.TRIANGLES(:,1:3);
num_matrix = length(matrix_nodes1);

[A,Index] = sort(sum(matrix_pos1,2));
matrix_pos = matrix_pos1(Index,:);

matrix_nodes = matrix_nodes1;

for i=1:length(Index)
    a = find(matrix_nodes1(:,1)==Index(i,1));
    b = find(matrix_nodes1(:,2)==Index(i,1));
    c = find(matrix_nodes1(:,3)==Index(i,1));
    matrix_nodes(a,1) = i;
    matrix_nodes(b,2) = i;
    matrix_nodes(c,3) = i;
end

sumnodes = sum(matrix_nodes,2);
[A,I] = sort(sumnodes);
matrix_nodes = matrix_nodes(I,:);

% Finding the center of the triangles
[center_matrix, matrix_id, matrix_area] = centerNodes(num_matrix, matrix_pos, matrix_nodes);

% Finding the center of connectivities and 1D-area
[center_fracture, fracture_id, frac_area, fracture_nodes_pairs] = centerLines(num_matrix, matrix_nodes, matrix_pos, matrix_id);

% Initiation of fracture
% Find fracture
% x_bor = 5;
% y_bor = [2,8];
% ab = 1;
% a = zeros(200,1);
% for i=1:length(center_fracture)
%     if center_fracture(i,1) <= x_bor + 0.025 && center_fracture(i,1) >= x_bor - 0.025
%         if center_fracture(i,2) <= y_bor(2) && center_fracture(i,2) >= y_bor(1)
%             a(ab) = fracture_id(i);
%             ab = ab + 1;
%         end
%     end
% end
% a(a==0) = [];
% y_bor = 5;
% x_bor = [2.0,8.0];
% ab = 1;
% b = zeros(200,1);
% for i=1:length(center_fracture)
%     if center_fracture(i,1) <= x_bor(2) && center_fracture(i,1) >= x_bor(1)
%         if center_fracture(i,2) <= y_bor + 0.025 && center_fracture(i,2) >= y_bor - 0.025
%             b(ab) = fracture_id(i);
%             ab = ab + 1;
%         end
%     end
% end
% b(b==0) = [];
% c = [a;b];

% x_bor = 5;
% y_bor = [2,8];
% ab = 1;
% c = zeros(200,1);
% for i=1:length(center_fracture)
%     if center_fracture(i,1) <= x_bor + 0.025 && center_fracture(i,1) >= x_bor - 0.025
%         if center_fracture(i,2) <= y_bor(2) && center_fracture(i,2) >= y_bor(1)
%             c(ab) = fracture_id(i);
%             ab = ab + 1;
%         end
%     end
% end
% c(c==0) = [];

% frac = [];
% frac = [161, 182, 201, 223, 251, 278, 299, 314];
% frac = [210, 231,250,273,306,325,344,188,214,244, 272,304,341,362,380];
% frac = [1148,1100,1050,1006,974,948,931,924,900,916,952,971,993,1025,1064,1119,1164,1215,1263,...
%     1282,1315,1321,1310,1300,1156,1083,1040,982,929,876,835,796,764,741];
% frac = [7792	4421	7700	10433	6399	9933	6440	7180	4579	6904	9058	4786	8309	4648	7130	9464	9493	4929	7462	9106	7433	4323	8500	5472	8396	8909	4876	4775	9913];
% frac = [411,444,490,526,559];
% frac = c;% cmg
frac = fracture_line_id();
init_fracture = frac - num_matrix;

% change fracture_id
new_fracture_nodes_pairs = zeros(length(init_fracture),2);
new_center_fracture = zeros(length(init_fracture),2);
new_fracture_id = zeros(length(init_fracture),1);

new_fracture_id(:,1) = linspace(1+num_matrix, length(init_fracture)+num_matrix,length(init_fracture));
new_frac_area = frac_area(init_fracture,1);
new_center_fracture(:,1) = center_fracture(init_fracture,1);
new_center_fracture(:,2) = center_fracture(init_fracture,2);
new_fracture_nodes_pairs(:,1) = fracture_nodes_pairs(init_fracture,1);
new_fracture_nodes_pairs(:,2) = fracture_nodes_pairs(init_fracture,2);

% Merge_matrix_fracture id
matrix_frac_id = [matrix_id;new_fracture_id];
matrix_frac_center = [center_matrix;new_center_fracture];
matrix_frac_area = [matrix_area; new_frac_area];
aperture = 0.01; %in m
new_frac_areaXaperture = new_frac_area*aperture;
matrix_frac_aper_area = [matrix_area; new_frac_areaXaperture];
con_frac = [new_fracture_nodes_pairs, zeros(length(new_fracture_nodes_pairs),1)];
matrix_frac_nodes = [matrix_nodes;con_frac];

% Define maximum intersection for fracture
max_intersection = 5;

tot_id = length(matrix_frac_id);
% ---------------------- units ----------------
% Iunit = 'FIELD'
Alpha=0.001127;
vol_unit=5.615; % (ft3 -> bbl)
% Time -> days
% distance  -> foot
% rates -> bbl/day
% volume -> bbl
% permeability -> mD
% viscoity -> cP
% pressure -> psi
% --------------------   Discretizations in x-, y- and z- directions --------------------
% Dx=zeros(Nx,1);
% Dy=zeros(Ny,1);
%Dz=zeros(Nz,1);
% --------------------   permeabilities  --------------------
Permx=zeros(tot_id,1);
% --------------------   Wells   --------------------
Wells(1:tot_id,1)=struct('id',0,'bhp',0,'rate',0);
% Note for wells
% if id = 1 => the well is an injector and rate should be specified
% if id = -1 => the well is a producer and the bhp should be specified

% ------------------ case 1 ---------
% input grid dimensions fracture
frac_id = matrix_frac_id(num_matrix+1:end);
x_frac_init = matrix_frac_center(frac_id,1);
y_frac_init = matrix_frac_center(frac_id,2);

% input permeabilities (mD)
Permx(1:end,1)=1; Permx(frac_id,1)=100000; 

% porosity
por = zeros(tot_id,1);
por(1:end,1) = 0.2;
por(num_matrix+1:end,1) = 1;

% Input for an injector (rate at res. conditions)
Wells(1,1).id=1;        %Coordinate x = 
Wells(1,1).rate=0.15;


% Input for a producer 
Wells(num_matrix,1).id=-1;
Wells(num_matrix,1).bhp=3500;



% % Input for an injector (rate at res. conditions)
% Wells(637,1).id=1;
% Wells(637,1).rate=0.15;
% 
% 
% % Input for a producer 
% Wells(824,1).id=-1;
% Wells(824,1).bhp=3500;

% --------------------- fluid viscosities ---------------------
% oil viscosity (cP)
mu_o = 5;  
% water viscosity (cP)
mu_w = 1.0;
% ----------------------- relative permeabilities ---------------------
% Corey.xls
swir =	0.250 ;	sorw=	0.300 ;
wexp=	2.500 ;	oexp=	2.000 ;
krwro=	0.400 ;	krocw=	0.700 ;

% ===========================  end of input =================================

% convert to FIELD unit
Permx=Permx.*Alpha;

% well PI (don't change) transmisibility between well and block -> affects
% BHP. No issue around the well. 
PI=1000;

% Memory allocation 
%Pressure
Pij=zeros(tot_id,1);

% water saturation (old timestep)
Swij=zeros(tot_id,1);
% water saturation (new timestep)
Swij_n=zeros(tot_id,1);

% ------------- initial water saturation (could be changed), Water
% saturation is not changing in this case
Swij(:,1)= swir;
Swij_n(:,1) = Swij;


cumwi=0; %cummulative water injector
cumwp=0; %cummulative water producer
% memory allocation for mobilities at the cell interfaces 
Mob=zeros(tot_id, 4+max_intersection);
% memory allocation for mobilities at the cell centers for wells 
WellMobW=zeros(tot_id,1);
WellMobO=zeros(tot_id,1);

% memory allocation for fluxes at the cell interfaces 
Fluxw=zeros(tot_id,4+max_intersection);

% memory allocation for production rates for wells
Prodw=zeros(tot_id,1);
ProdRates=zeros(2,1);

% memory allocation for linear system assembly and solution 
Jac=sparse(tot_id,tot_id);
Rhs=zeros(tot_id,1);
P=zeros(tot_id,1);

% ----------------- time step control -----------------
% initial timestep (days)
dt=0.001;
% minimum allowed time step
dtmin=0.0000001;
% maximum allowed time step
dtmax=0.001;
% minimum sw change before increase timestep
dswmin=0.1;
% maximum sw change before cutting timestep
dswmax=0.2;
% ------------------end timestep control --------------
% simulation time
final_time =6.0;
 
time = 0;
itime=0;
iplot=0;

% for CPUtime calculation
t = cputime;

% this functions returns the volumetrics 
% OOIP : original oil in place 
% OWIP : original water in place
% totpv: total pore-volume
% PorVol: porvolume map for all gridblocks

[ooip,owip,totpv,PorVol]=CalculateVolumetrics(matrix_frac_aper_area,Swij,por,vol_unit);

%  -------  calculate transmisibilities -----------------------
[Trans]=CalculateTran(matrix_frac_center,matrix_nodes,matrix_frac_nodes, new_fracture_nodes_pairs, num_matrix,new_fracture_id,matrix_pos,aperture,Permx, max_intersection);

% figure to plot SW vs. time
figure;
hold on;

% start time loop
inum = 1;
while time < final_time
    % counter for time iterations 
    itime=itime+1;
    
    % calculate mobilities at interfaces lambda * transmisibility
    [Mob]=CalculateMob(matrix_frac_center,matrix_frac_nodes,max_intersection,Pij,Trans,num_matrix,matrix_nodes,new_fracture_nodes_pairs,new_fracture_id,Mob,Swij,swir,wexp,krwro,sorw,oexp,krocw,mu_o,mu_w);

    % calculate mobilities for producers    
    [WellMobW,WellMobO]=CalculateWellMob(matrix_frac_center,num_matrix,Swij,swir,wexp,krwro,sorw,oexp,krocw,mu_o,mu_w,Wells);
    
    % solve pressure equation
    [P,Pij]=solveP(matrix_frac_center,Mob,num_matrix,matrix_nodes,new_fracture_nodes_pairs,new_fracture_id,Wells,PI,WellMobW,WellMobO,Jac,Rhs);
    
    % solve fluxes and rates at production wells
    [Fluxw,Prodw,ProdRates]=CalculateWflux(matrix_frac_center,matrix_frac_nodes,matrix_nodes,new_fracture_nodes_pairs,num_matrix,new_fracture_id,Pij,Trans,Swij,swir,wexp,krwro,sorw,oexp,krocw,mu_o,mu_w,Wells,PI,Fluxw,Prodw,ProdRates);    
   
    
    % calculate water saturation
%     [Swij_n]=CalculateSw(matrix_frac_center,matrix_nodes,new_fracture_nodes_pairs,num_matrix,Fluxw,new_fracture_id,Prodw,Swij,Swij_n,dt,Wells,PorVol);
    [Swij_n]=CalculateSw(matrix_frac_center,Fluxw,Prodw,Swij,Swij_n,dt,Wells,PorVol);

    
    % check for stability and timestep control
    % if icut = 1 => the timestep is too large => cut time step and repeat
    % if icut = 0 => time step is ok => procced
    [dtnew,icut]=DtControl(matrix_frac_center,dt,Swij_n,Swij,dtmin,dtmax,dswmin,dswmax);
    
    % if icut=1 i.e., we decided to repeat calculations
    if(icut==1)
        % this loop will cutting timestep and repeating calcualtions unitil
        % icut becomes 0
        while(icut==1)
                % track time step
                dt=dtnew;
                % repeat calcualtion 
                [Swij_n]=CalculateSw(matrix_frac_center,Fluxw,Prodw,Swij,Swij_n,dt,Wells,PorVol);
                % check for stability again
                [dtnew,icut]=DtControl(matrix_frac_center,dt,Swij_n,Swij,dtmin,dtmax,dswmin,dswmax);
        end        
    end
    
%     save_time(inum) = time;
%     saveSw(inum) = Swij(num_matrix);
%     saveP(inum) = P(num_matrix);
%     inum = inum + 1;
    
    % advance time
    time = time +dt     ;
    
    % update water saturation
    Swij=Swij_n;
    
    % check for material balance and calculate injection and production
    % volumes
    [wip,cumwi,cumwp,err]=Materialbalance(matrix_frac_center,Swij,Wells,Prodw,PorVol,dt,cumwi,cumwp,owip);   
    err
    % print outputs every 10 iteration
    if(rem(itime,10)==0)
         iplot=iplot+1;
         plotp_Sw_time(matrix_frac_nodes,matrix_pos,Swij,time);
         alldt(iplot,1:2)=[time,dt];
         rates(iplot,1:3)=[time,ProdRates(1),ProdRates(2)];
         Bal(iplot,1:5)=[time,wip,cumwi,cumwp,err];
    end
    dt=dtnew;
    if(final_time-time<dt)dt=max(final_time-time,1e-10);end
end
% plot saturation 
% save('Saturation_lastNodes.mat', 'saveSw');
% save('Pressure_lastNodes.mat', 'saveP');
% save('time.mat', 'save_time');
plotp_Sw_time(matrix_frac_nodes,matrix_pos,Swij,time)
% plot rates
figure
plot(alldt(:,1),alldt(:,2),'-o');
title('timesteps versus time');
figure
hold on;
plot(rates(:,1),rates(:,2),'-ob');
plot(rates(:,1),rates(:,3),'-sg');
title('Production rates (RB/day)');
% display injected PV
PVI=cumwi/totpv*100
CPU_time = cputime-t
%Sw_P_vs_T
figure;
%plot(Bal(:,1),Bal(:,5))
%cumwp
%cumwi
%wip
%owip
% display total error
TotalMaterialBalanceError = Bal(iplot,5)
% plot final SW
plotp_Sw(matrix_frac_nodes, matrix_pos,matrix_frac_center,aperture,Swij,time);
% plot pressure
plotp_noBC(matrix_frac_nodes, matrix_pos,matrix_frac_center,aperture,Permx,P,time);

% *************************************** Solve the pressure equation *****************************************
function [P,Pij]=solveP_noBC(Nx,Ny,Tranx,Trany,Wells,PI,WellMobW,WellMobO,Jac,Rhs)

Nt=Nx*Ny;
ic=0;
for j=1:Ny
    for i=1:Nx
        ic=ic+1;
        Jac(ic,ic)=0;
        Rhs(ic)=0;
        ie=ic+1;
        iw=ic-1;
        in=ic-Nx;
        is=ic+Nx;
        % if producer (i.e. pressure BC)
        if Wells(i,j).id==-1
            Jac(ic,ic)=Jac(ic,ic)+PI*(WellMobW(i,j)+WellMobO(i,j));
            Rhs(ic)=Rhs(ic)+PI*(WellMobW(i,j)+WellMobO(i,j))*Wells(i,j).bhp;
            % else if there is an injector, or no well
        end
        Rhs(ic)=Rhs(ic)+Wells(i,j).rate;
        %
        if(i<Nx)  % East side
            Jac(ic,ie)=-Tranx(i+1,j); Jac(ie,ic)=Jac(ic,ie);
            Jac(ic,ic)=Jac(ic,ic)-Jac(ic,ie);
        end
        
        if(i>1) % West side
            Jac(ic,iw)=-Tranx(i,j);Jac(iw,ic)=Jac(ic,iw);
            Jac(ic,ic)=Jac(ic,ic)-Jac(ic,iw);
        end
        
        if(j>1) % north side
            Jac(ic,in)=-Trany(i,j);Jac(in,ic)=Jac(ic,in);
            Jac(ic,ic)=Jac(ic,ic)-Jac(ic,in);
        end
        if(j<Ny) % south side
            Jac(ic,is)=-Trany(i,j+1);Jac(is,ic)=Jac(ic,is);
            Jac(ic,ic)=Jac(ic,ic)-Jac(ic,is);
        end
        
    end
end
%Jac
%cond(Jac)
%Rhs
%solve the system
%size(Jac);
%size(Rhs);
%spy(Jac);

%  solve the linear system
P=Jac\Rhs;


icnt=0;
for j=1:Ny
    for i=1:Nx
        icnt=icnt+1;
        Pij(i,j)=P(icnt);
    end
end



% ----------------------------------------------
function plotp_noBC(matrix_frac_nodes, matrix_pos,matrix_frac_center,aperture,Permx,P,time)

F = scatteredInterpolant(matrix_frac_center(:,1),matrix_frac_center(:,2),P);
[xq,yq] = meshgrid(meshgrid(0:0.2:max(max(matrix_pos(:,1)),max(matrix_pos(:,2)))));
F.Method = 'linear';
vq2 = F(xq,yq);

% ---------------plot P:  Method 1 using surf -----
figure
colormap jet;
surf(xq,yq,vq2)
c=colorbar;
c.Label.String = 'Pressure';
txt=['Time = ' num2str(time) ' days'];
title(txt);

% ---------------plot P:  Method 2 using contour -----

figure
colormap jet;
contourf(xq,yq,vq2)
c=colorbar;
c.Label.String = 'Pressure';
txt=['Time = ' num2str(time) ' days'];
title(txt);

% ---------------plot P:  Method 3 using fill (used in RS) -----

figure;
hold on
for i=1:length(matrix_frac_center)
    % get cell centered pressure
    aas = matrix_frac_nodes(i,:);
    aas(aas==0) = [];
    n = length(aas);
    Xn = [];
    Yn = [];
    if n > 2
        Xn(1:n)= matrix_pos(matrix_frac_nodes(i,:),1);
        Yn(1:n)= matrix_pos(matrix_frac_nodes(i,:),2);
    else
        asn = matrix_frac_nodes(i,:);
        asn(asn==0) = [];
        gh1 = matrix_pos(asn,1);
        gh2 = matrix_pos(asn,2);
        Xn(1:4) = [gh1(1)- aperture*0.5, gh1(1) + aperture*0.5, gh1(2) + aperture*0.5,gh1(2)- aperture*0.5];
        Yn(1:4) = [gh2(1),gh2(1), gh2(2), gh2(2)];
    end
    fill(Xn,Yn,P(i,1));
end
colormap jet;
c=colorbar;
c.Label.String = 'Pressure';
txt=['Time = ' num2str(time) ' days'];
title(txt);


%  ------------------- plot velocity field ---------------
%return
% N = length(matrix_frac_center);
% Vx=zeros(N,1);
% Vy=zeros(Nx,Ny);
% for j=1:Ny
%     for i=1:Nx
%         % x-direction
%         if(i<Nx)
%             pe=(pij(i+1,j)+pij(i,j))*0.5;
%         else
%             pe=0;
%         end
%         if(i>1)
%             pw=(pij(i-1,j)+pij(i,j))*0.5;
%         else
%             pw=0;
%         end
%         if(pe>0 && pw>0)
%             Vx(i,j)=-Permx(i,j)*(pe-pw)/Dx(i);
%         end
%         % y-direction
%         if(j<Ny)
%             pn=(pij(i,j+1)+pij(i,j))*0.5;
%         else
%             pn=0;
%         end
%         if(j>1)
%             ps=(pij(i,j-1)+pij(i,j))*0.5;
%         else
%             ps=0;
%         end
%         if(ps>0 && pn>0)
%             Vy(i,j)=-Permy(i,j)*(ps-pn)/Dy(j);
%         end
%         
%     end
% end
% if(Ny>1&& Nx>1)
% h=quiver(Xv,Yv,Vx',Vy','k');
% set(h,'LineWidth',1.5)
% end

% ----------------------------------------------
function plotp_Sw_time(matrix_frac_nodes,matrix_pos,Swij,time)

T = delaunay(matrix_pos(:,1),matrix_pos(:,2));
aa = zeros(length(matrix_pos(:,1)),1);
for i=1:length(matrix_pos(:,1))
    aa(i,1) = mean(Swij(logical(sum(ismember(matrix_frac_nodes,i),2))));
end

colormap gray;
trisurf(T,matrix_pos(:,1),matrix_pos(:,2),aa)
c=colorbar;
c.Label.String = 'Sw';
c.Limits = [0.25 0.7];
caxis([0.25 0.7])
% view([30,50])
view(2)
txt=['Time = ' num2str(time) ' days'];
title(txt);
  hold off
  drawnow


% ----------------------------------------------
function plotp_Sw(matrix_frac_nodes, matrix_pos,matrix_frac_center,aperture,Swij,time)
figure;
hold on
for i=1:length(matrix_frac_center)
    % get cell centered pressure
    aas = matrix_frac_nodes(i,:);
    aas(aas==0) = [];
    n = length(aas);
    Xn = [];
    Yn = [];
    if n > 2
        Xn(1:n)= matrix_pos(matrix_frac_nodes(i,:),1);
        Yn(1:n)= matrix_pos(matrix_frac_nodes(i,:),2);
    else
        asn = matrix_frac_nodes(i,:);
        asn(asn==0) = [];
        gh1 = matrix_pos(asn,1);
        gh2 = matrix_pos(asn,2);
        Xn(1:4) = [gh1(1)- aperture*0.5, gh1(1) + aperture*0.5, gh1(2) + aperture*0.5,gh1(2)- aperture*0.5];
        Yn(1:4) = [gh2(1),gh2(1), gh2(2), gh2(2)];
    end
    fill(Xn,Yn,Swij(i,1));
end
colormap jet;
c=colorbar;
c.Label.String = 'Sw';
c.Limits = [0.25 0.7];
caxis([0.25 0.7])
txt=['Time = ' num2str(time) 'days'];
title(txt);

% figure
% plot(matrix_frac_center(:,1),Swij(:,1),'-o')

%Swij(:,1)

