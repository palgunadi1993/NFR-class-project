function twophaseflow

% domain size is determined by number and size of gridcells
% ========================= Input Data ==================================
% --------------------    Grid Dimensions   -------------------------
Nx=40;   % nb of cells in the x-direction
Ny=40;   % nb of cells in the y-direction
Nz=1;     % nb of cells in the z-direction
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
Dx=zeros(Nx,1);
Dy=zeros(Ny,1);
%Dz=zeros(Nz,1);
% --------------------   permeabilities  --------------------
Permx=zeros(Nx,Ny);
Permy=zeros(Nx,Ny);
% --------------------   Wells   --------------------
Wells(1:Nx,1:Ny)=struct('id',0,'bhp',0,'rate',0);
% Note for wells
% if id = 1 => the well is an injector and rate should be specified
% if id = -1 => the well is a producer and the bhp should be specified

% ------------------ case 1 ---------
% input grid dimensions 
Dx(1:Nx) =1; Dx(20)=0.8;  
Dy(1:Ny) =1; Dy(20)=0.8;  

% input permeabilities (mD)
Permx(1:Nx,1:Ny)=1; Permx(5:35,20)=500; 
Permy(1:Nx,1:Ny)=1; Permy(20,5:35)=500;

% porosity
por=0.2;

% Input for an injector (rate at res. conditions)
Wells(1,1).id=1;
Wells(1,1).rate=0.15;


% Input for a producer 
Wells(Nx,Ny).id=-1;
Wells(Nx,Ny).bhp=3500;

% --------------------- fluid viscosities ---------------------
% oil viscosity (cP)
mu_o = 5;  
% water viscosity (cP)
mu_w = 1.0;
% ----------------------- relative permeabilities ---------------------
%
swir =	0.250 ;	sorw=	0.300 ;
wexp=	2.500 ;	oexp=	2.000 ;
krwro=	0.400 ;	krocw=	0.700 ;

% ===========================  end of input =================================

% convert to FIELD unit
Permx=Permx.*Alpha;
Permy=Permy.*Alpha;

% well PI (don't change)
PI=1000;

% Memory allocation 
%Pressure
Pij=zeros(Nx,Ny);

% water saturation (old timestep)
Swij=zeros(Nx,Ny);
% water saturation (new timestep)
Swij_n=zeros(Nx,Ny);

% ------------- initial water saturation (could be changed)
Swij(:,:)= swir;
Swij_n(:,:) = Swij;


cumwi=0;
cumwp=0;
% memory allocation for mobilities at the cell interfaces 
Mobx=zeros(Nx+1,Ny);
Moby=zeros(Nx,Ny+1);
% memory allocation for mobilities at the cell centers for wells 
WellMobW=zeros(Nx,Ny);
WellMobO=zeros(Nx,Ny);

% memory allocation for fluxes at the cell interfaces 
Fluxw_x=zeros(Nx+1,Ny);
Fluxw_y=zeros(Nx,Ny+1);

% memory allocation for production rates for wells
Prodw=zeros(Nx,Ny);
ProdRates=zeros(2,1);

% memory allocation for linear system assembly and solution 
Jac=zeros(Nx*Ny,Nx*Ny);
Rhs=zeros(Nx*Ny,1);
P=zeros(Nx*Ny);

% ----------------- time step control -----------------
% initial timestep (days)
dt=0.15;
% minimum allowed time step
dtmin=0.001;
% maximum allowed time step
dtmax=0.8;
% minimum sw change before increase timestep
dswmin=0.1;
% maximum sw change before cutting timestep
dswmax=0.2;
% ------------------end timestep control --------------
% simulation time
final_time =80.0;
 
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

[ooip,owip,totpv,PorVol]=CalculateVolumetrics(Nx,Ny,Swij,Dx,Dy,por,vol_unit);

%  -------  calculate transmisibilities -----------------------
[Tranx, Trany]=CalculateTran(Nx,Ny,Dx,Dy,Permx,Permy);

% figure to plot SW vs. time
figure;
hold on;

% start time loop
while time < final_time
    % counter for time iterations 
    itime=itime+1;
    
    % calculate mobilities at interfaces
    [Mobx, Moby]=CalculateMob(Nx,Ny,Pij,Tranx,Trany,Mobx,Moby,Swij,swir,wexp,krwro,sorw,oexp,krocw,mu_o,mu_w);
    
    % calculate mobilities for producers    
    [WellMobW,WellMobO]=CalculateWellMob(Nx,Ny,Swij,swir,wexp,krwro,sorw,oexp,krocw,mu_o,mu_w,Wells);
    
    % solve pressure equation
    [P,Pij]=solveP(Nx,Ny,Mobx,Moby,Wells,PI,WellMobW,WellMobO,Jac,Rhs);
    
    % solve fluxes and rates at production wells
    [Fluxw_x,Fluxw_y,Prodw,ProdRates]=CalculateWflux(Nx,Ny,Pij,Tranx,Trany,Swij,swir,wexp,krwro,sorw,oexp,krocw,mu_o,mu_w,Wells,PI,WellMobW,Fluxw_x,Fluxw_y,Prodw,ProdRates);    
    
    % calculate water saturation
    [Swij_n]=CalculateSw(Nx,Ny,Fluxw_x,Fluxw_y,Prodw,Swij,Swij_n,dt,Wells,PorVol);

    % check for stability and timestep control
    % if icut = 1 => the timestep is too large => cut time step and repeat
    % if icut = 0 => time step is ok => procced
    [dtnew,icut]=DtControl(Nx,Ny,dt,Swij_n,Swij,dtmin,dtmax,dswmin,dswmax);
    
    % if icut=1 i.e., we decided to repeat calculations
    if(icut==1)
        % this loop will cutting timestep and repeating calcualtions unitil
        % icut becomes 0
        while(icut==1)
                % track time step
                dt=dtnew;
                % repeat calcualtion 
                [Swij_n]=CalculateSw(Nx,Ny,Fluxw_x,Fluxw_y,Prodw,Swij,Swij_n,dt,Wells,PorVol);
                % check for stability again
                [dtnew,icut]=DtControl(Nx,Ny,dt,Swij_n,Swij,dtmin,dtmax,dswmin,dswmax);
        end        
    end
    
    % advance time
    time = time +dt ;    
    
    % update water saturation
    Swij=Swij_n;
    
    % check for material balance and calculate injection and production
    % volumes
    [wip,cumwi,cumwp,err]=Materialbalance (Nx,Ny,Swij,Wells,Prodw,PorVol,dt,cumwi,cumwp,owip);   
    
    % print outputs every 10 iteration
    if(rem(itime,10)==0)
         iplot=iplot+1;
         plotp_Sw_time(Nx,Ny,Dx,Dy,Swij,time);
         alldt(iplot,1:2)=[time,dt];
         rates(iplot,1:3)=[time,ProdRates(1),ProdRates(2)];
         Bal(iplot,1:5)=[time,wip,cumwi,cumwp,err];
    end
    dt=dtnew;
    if(final_time-time<dt)dt=max(final_time-time,1e-10);end
end
% plot saturation 
plotp_Sw_time(Nx,Ny,Dx,Dy,Swij,time)
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
plotp_Sw(Nx,Ny,Dx,Dy,Swij,time);
% plot pressure
plotp_noBC(Nx,Ny,Dx,Dy,Permx,Permy,P,time);

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
function plotp_noBC(Nx,Ny,Dx,Dy,Permx,Permy,P,time)


% ---------------plot P:  Method 1 using surf -----
icnt=0;
for j=1:Ny
    for i=1:Nx
        icnt=icnt+1;
        pij(i,j)=P(icnt);
    end
end



% ---------------plot P:  Method 2 using contour -----

if(Ny>1&& Nx>1)
figure
colormap jet;
surf(pij)
c=colorbar;
c.Label.String = 'Pressure';
txt=['Time = ' num2str(time) ' days'];
title(txt);

figure
colormap jet;
contourf(pij,20)
c=colorbar;
c.Label.String = 'Pressure';
txt=['Time = ' num2str(time) ' days'];
title(txt);
else
pij(:,1)
end

% ---------------plot P:  Method 3 using fill (used in RS) -----

h=figure;
hold on;
icnt=0;
ly=0;
icnt=0;
Ly=sum(Dy(1:Ny));
for j=1:Ny
    ly=ly+Dy(j);
    lx=0;
    for i=1:Nx
        % get cell centered pressure
        icnt=icnt+1;
        lx=lx+Dx(i);
        Xn(1:4)=[lx-Dx(i),lx-Dx(i),lx,lx];
        Yn(1:4)=Ly-[ly-Dy(j),ly,ly,ly-Dy(j)];
        fill(Xn,Yn,P(icnt));
        %patch(Xn,Yn,P(icnt));
    end
end
colormap jet;
c=colorbar;
c.Label.String = 'Pressure';
txt=['Time = ' num2str(time) ' days'];
title(txt);


% ------- plot velocity field
% in x-direction
icnt=0;
Xv=zeros(Nx,1);
lx=0;
for i=1:Nx
    lx=lx+Dx(i);
    Xv(i)=lx-Dx(i)/2;
end

Yv=zeros(Ny,1);
ly=0;
Ly=sum(Dy(1:Ny));
for j=1:Ny
    ly=ly+Dy(j);
    Yv(j)=ly-Dy(j)/2;
end
Yv=Ly-Yv;

%  ------------------- plot velocity field ---------------
%return
Vx=zeros(Nx,Ny);
Vy=zeros(Nx,Ny);
for j=1:Ny
    for i=1:Nx
        % x-direction
        if(i<Nx)
            pe=(pij(i+1,j)+pij(i,j))*0.5;
        else
            pe=0;
        end
        if(i>1)
            pw=(pij(i-1,j)+pij(i,j))*0.5;
        else
            pw=0;
        end
        if(pe>0 && pw>0)
            Vx(i,j)=-Permx(i,j)*(pe-pw)/Dx(i);
        end
        % y-direction
        if(j<Ny)
            pn=(pij(i,j+1)+pij(i,j))*0.5;
        else
            pn=0;
        end
        if(j>1)
            ps=(pij(i,j-1)+pij(i,j))*0.5;
        else
            ps=0;
        end
        if(ps>0 && pn>0)
            Vy(i,j)=-Permy(i,j)*(ps-pn)/Dy(j);
        end
        
    end
end
if(Ny>1&& Nx>1)
h=quiver(Xv,Yv,Vx',Vy','k');
set(h,'LineWidth',1.5)
end

% ----------------------------------------------
function plotp_Sw_time(Nx,Ny,Dx,Dy,Swij,time)



icnt=0;


colormap jet;
contourf(Swij,20)
c=colorbar;
c.Label.String = 'Sw';
c.Limits = [0.25 0.7];
txt=['Time = ' num2str(time) ' days'];
title(txt);
  hold off
  drawnow

% ----------------------------------------------
function plotp_Sw(Nx,Ny,Dx,Dy,Swij,time)


hold on;
icnt=0;
ly=0;
icnt=0;
Ly=sum(Dy(1:Ny));
LX=zeros(Nx,1);
for j=1:Ny
    ly=ly+Dy(j);
    lx=0;
    for i=1:Nx
        % get cell centered pressure
        icnt=icnt+1;
        lx=lx+Dx(i);
        Xn(1:4)=[lx-Dx(i),lx-Dx(i),lx,lx];
        Yn(1:4)=Ly-[ly-Dy(j),ly,ly,ly-Dy(j)];
        fill(Xn,Yn,Swij(i,j));
        LX(i)=lx-Dx(i)/2;
        %patch(Xn,Yn,P(icnt));
    end
end
colormap jet;
c=colorbar;
c.Label.String = 'Sw';
c.Limits = [0.25 0.7];
txt=['Time = ' num2str(time) 'days'];
title(txt);


if(Ny>1&& Nx>1)
% figure
% colormap jet;
% contourf(Swij,20)
% c=colorbar;
% c.Label.String = 'Sw';
% c.Limits = [0.25 0.7];
else
    figure
    plot(LX,Swij(:,1),'-o')
end

%Swij(:,1)

