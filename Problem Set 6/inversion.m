% Code adapted from Jessica Badgely's transient inversion code

steps = [1]

%Run Options{{{
%clustername = 'discovery';
clustername = 'totten';
%clustername = 'pfe';
%clustername = oshostname();
max_run_time = 4; %discovery is in hours, pfe is in minutes

%Some of the following are overwritten depending on the chosen cluster and step
loadonly=0;
interactive=1;
waitonlock=0;

domain = 'CW'; %either a string or a list of catchment indices
domain_name = domain;

folder = 'Models'; %['Models_region_' domain]; % 

fric_type = 1; % fric_type = 1 is Budd m=1 - DEFAULT, Schoof m=1/3 is also possible now with fric_type = 2;
fric_init = 'cal'; % 'cal' for original calibration and 'Budd' for based on the Budd inversion (only applies to fric_type>1)

%icefront_type = 0; % specified with monthly data - DEFAULT
icefront_type = 5; % 1-12 means specified with annual data for that month - DEFAULT is 5 (on average most extensive Greene et al., 2014)
%icefront_type = 13; % static ice front

forecast_forcing ='MIROC5-rcp26-med'; %'control';% 'MIROC5-rcp85-med'; % 

catchment=208; %218; %221; %208; %218;

indep_type=0; %constant solution for Transient calibration - DEFAULT
%indep_type=1; %one solution per year for Transient calibration (if have non whole number of years, then have an additional one at the end)
%indep_type=2; %one solution per season (same solution every year) for Transient calirbation
%indep_type=3; %strategic seasonal solutions (same solution every year) for Transient calibration

bfc_type=0; % set basal friction coefficient to output from indep_type=0;
%bfc_type=2; % set basal friction coefficient to output from indep_type=2;

vel_min=500; %m/yr
vel_max=100000; %m/yr
dist_max=80; %km

future_SMB_test = 'ISMIP6'; %'RACMO'; %Choices are ISMIP6 and RACMO
future_Retreat_test = 'ISMIP6'; %'Greene'; %'ISMIP6'; %Choices are ISMIP6 and Greene

% THINGS THAT DEPEND ON THE DOMAIN
use_exp_for_ice_levelset = false; % if not false, provide string of path to exp file of edited ice levelset 0 contour
use_exp_for_mesh_boundary = ['Exp/' domain '_domain_extended.exp']; %false; % if not false, provide string of path to exp file of edited mesh boundary
use_exp_for_spcvxvy = false; 
%use_exp_for_snapshot_inversion = false;
constant_levelset_catchment_ids = false; % if not false, enter list of catchment ids
soften_shear_margins = []; %
mesh_buffer = 3000; % DEFAULT is 3000
num_catchments = 0;
catchment_ids = false;
if strcmp(domain,'SE')
	softer_shear_margins = [63];
	catchment_ids = [63,64,65,66,67,68,73,79,80,83,88,89,90,91,93,94,106,107,108,110,111,112,113,114,115,116,117,118,119,120,121,125,126,128,186,187,188,189,190,191,192,193,199,200,201,202,209,210,211,212,213,231,232,233,237,243,257,258,259];
	use_exp_for_spcvxvy = 'Exp/SE_spcvxvy.exp';
	num_catchments = 59;
	%use_exp_for_friction_smoothing = false;
	%use_exp_for_ice_levelset = 'Exp/SE_ice_levelset.exp';
	%constant_levelset_catchment_ids = [106]; % also potentially 94
elseif strcmp(domain,'NO')
	use_exp_for_spcvxvy = 'Exp/NO_spcvxvy.exp';
	soften_shear_margins = [223];
   catchment_ids = [39,50,51,52,53,54,55,56,57,77,84,109,215,216,146,176,223,224,];
elseif strcmp(domain,'CW')
%  I commented this to prevent reference to the Mouginot shp file which is locked (Hugh, 1/15)	
%	soften_shear_margins = [218];
%	use_exp_for_spcvxvy = 'PNAS_2024/CW_spcvxvy.exp';
	%use_exp_for_friction_smoothing = 'Exp/CW_spcvxvy_07232024.exp';
	num_catchments = 17;
   catchment_ids = [3 4 5 6 7 8 9 10 72 81 82 85 218 219 220 221 222];
elseif strcmp(domain,'NW')
	num_catchments = 58;
   catchment_ids = [1 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 41 42 43 44 45 46 47 48 49 78 86 103 104 105 127 147 148 160 162 164 165 167 172 174 175 178 180 181 183 185 195 196 197 207 208 214 225 238 239 240 242];
	use_exp_for_spcvxvy =  'Exp/NW_spcvxvy_11212024_v2.exp'; %'PNAS_2024/NW_spcvxvy.exp'; %'Exp/NW_spcvxvy_07262024.exp';
elseif strcmp(domain,'NE')
	use_exp_for_spcvxvy = 'Exp/NE_spcvxvy.exp';
   catchment_ids = [58,59,87,129,138,139,140,141,142,143,144,145,217,244,245,246,247];
elseif strcmp(domain,'CE')
	softer_shear_margins = [102];
	catchment_ids = [2,60,61,62,74,75,76,95,96,97,98,99,100,101,102,122,123,124,149,150,151,152,153,154,155,156,157,158,203,204,205,206,234,235,236,248,249,250,251,252,253,254,255];
	%mesh_buffer = 5000;
elseif strcmp(domain,'SW')
	%use_exp_for_snapshot_inversion = 'Exp/SW_snapshotInversion_keepUnchanged.exp';
	use_exp_for_spcvxvy = 'Exp/SW_spcvxvy.exp';
	num_catchments = 24;
   catchment_ids = [11,12,13,14,15,16,17,18,19,20,21,22,40,69,70,71,130,131,132,133,134,135,136,137]; 
   %use_exp_for_mesh_boundary = 'Exp/SW_domain_3000mBuffer_updated.exp';
elseif strcmp(domain,'NW_test_large')
	use_exp_for_mesh_boundary = ['Exp/NW_domain_extended_slightly.exp']; 
	use_exp_for_spcvxvy = 'PNAS_2024/NW_spcvxvy.exp';
end

% THINGS THAT DON'T CHANGE: If you change, MUST be noted

calibration_start_time = 2007; % DEFAULT is 2007
calibration_end_time = 2022; %2021.6; %2022;  % DEFAULT is 2022 %Note: This is now 2021.6 to align with when I have RACMO SMB data currently
projection_start_time = 2015; % DEFAULT is 2015
projection_end_time = 2100; % DEFAULT is 2100

niters_AD = 50; % DEFAULT is 50

indep_var_rheB = 0; %whether or not to include rheology B as an independent variable - DEFAULT is 0

%TODO: put these parameters into the model. In which steps? Perhaps dependent on when they are called in the following steps.

%}}}
%Cluster parameters{{{
if strcmpi(clustername,'totten'),
   cluster=generic('name',oshostname(),'np',30);
elseif strcmpi(clustername,'discovery'),
   cluster=discovery('numnodes',1,'cpuspernode',50,'time',max_run_time,'memory',256); %512); %256);
   interactive=0;
   waitonlock=0;
	path_noAD='/dartfs-hpc/rc/home/6/f006pz6/ISSM/bin/';
	path_AD='/dartfs-hpc/rc/home/6/f006pz6/ISSM-AD/bin/';
elseif strcmpi(clustername,'andes'),
   cluster=andes('numnodes',1,'cpuspernode',50,'time',max_run_time,'memory',128);
   interactive=0;
   waitonlock=0;
   path_noAD='/dartfs-hpc/rc/home/6/f006pz6/ISSM/bin/';
   path_AD='/dartfs-hpc/rc/home/6/f006pz6/ISSM-AD/bin/';
elseif strcmpi(clustername,'pfe'),
	%cluster=pfe('numnodes',1,'cpuspernode',28,'time',2*60,'processor','bro','queue','devel','grouplist','s2950'); %max time is 2 hours and can only run one at once
	%cluster=pfe('numnodes',2,'cpuspernode',28,'time',max_run_time,'processor','bro','queue','normal','grouplist','s2950'); %max time is 8 hours
	cluster=pfe('numnodes',2,'cpuspernode',28,'time',max_run_time,'processor','bro','queue','long','grouplist','s2950'); %max time is 120 hours
	interactive=0;
	waitonlock=0;
	path_noAD='/nobackup/jbadgele/ISSM/bin/';
   path_AD='/nobackup/jbadgele/ISSM-AD/bin/';	
else
   cluster=generic('name',oshostname(),'np',10);
end
%}}}

prefix = 'Model_';
org=organizer('repository',['./' folder],'prefix',prefix,'steps',steps); clear steps;

if perform(org,[domain_name '_Mesh']),% {{{

   disp('   -- Initializing the mesh');
   % set up the domain and initial mesh
   hmax=10e4; %10e3; %7e3; hugh changed 1/22
   hmin=10e3; %400
	if use_exp_for_mesh_boundary
		file_name = use_exp_for_mesh_boundary;
	else
      [file_name, domain_basin_IDs] = makeDomainBoundary(domain, './Exp', mesh_buffer, domain_name);
	end
   md = model();
   md=triangle(model,file_name,hmin);

   % using a priori analysis (observed velocity)
   disp('   -- Interpolating some data');
   vel_obs = interpFromMEaSUREsGeotiffs(md.mesh.x,md.mesh.y,[1994],[2016],'product',[670],'output',{'vx','vy'});
	velx = vel_obs.vx;
	vely = vel_obs.vy;
   %[velx vely] = interpJoughinCompositeGreenland(md.mesh.x,md.mesh.y);
   %M           = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask');

	% Testing this section needed to make ocean mask
   md.materials.rho_ice             = 917.;  % (917 is the value used in BedMachine)
   md.materials.rho_water           = 1027.; % ocean water (1027 is the value used in BedMachine)
   min_thickness                    = 10.; % minimum ice thickness used to setup the initial geometry
   min_surface                      = min_thickness*(1-md.materials.rho_ice/md.materials.rho_water) ;
   md.geometry.surface              = interpGimpdem(md.mesh.x, md.mesh.y) - interpBedmachineGreenland(md.mesh.x,md.mesh.y,'geoid');
   pos                              = find(md.geometry.surface<1.e-10); % ocean part or rocks
   md.geometry.surface(pos)         = min_surface; % set minimum ice surface on the ocean part
   md.geometry.bed                  = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'bed'); % interp method: linear
   md.geometry.base                 = zeros(length(md.geometry.bed),1); % initial setup
   floatation_base                  = md.geometry.surface*md.materials.rho_ice/(md.materials.rho_ice-md.materials.rho_water); % using the corrected surface
   pos                              = find(floatation_base<=md.geometry.bed); % grounded ice
   md.geometry.base(pos)            = md.geometry.bed(pos);
   pos                              = find(floatation_base>md.geometry.bed); % floating ice
   md.geometry.base(pos)            = floatation_base(pos);
   md.geometry.thickness            = md.geometry.surface-md.geometry.base;
   pos                              = find(md.geometry.thickness<min_thickness); % dealing with rocks or ocean part
   md.geometry.thickness(pos)       = min_thickness;
   md.geometry.surface(pos)         = md.geometry.thickness(pos)+md.geometry.base(pos);

   ice_masks = interpMonthlyIceMaskGreene(md.mesh.x, md.mesh.y, [calibration_start_time, calibration_start_time+1]);
   ice_levelset = ice_masks(1:end-1,1);
   %ice_levelset=-ones(size(md.mesh.x));
   %ice_levelset(find((M==0)|(M==1)|(M==4)))=1; % ocean and land from BedMachine
	%floating_ice_mask = zeros(size(md.mesh.x));
	%floating_ice_mask(find(M==3))=1;

	% Testing this section to capture fjords
   md.mask.ocean_levelset = sethydrostaticmask(md);
	%M = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask','nearest');
	ice_masks_noRock = interpMonthlyIceMaskGreene(md.mesh.x, md.mesh.y, [calibration_start_time, calibration_start_time+1],0);
	ice_levelset_noRock = ice_masks_noRock(1:end-1,1);
	pos = find((md.mask.ocean_levelset < 0) & (ice_levelset_noRock > 0)); % (M < 1));
   ice_levelset(pos) = 999999;

   pos=find(isnan(velx) | isnan(vely) | ice_levelset>0);
   velx(pos)=0; vely(pos)=0; vel=sqrt(velx.^2+vely.^2);

   md3dT = loadmodel('/totten_1/ModelData/ISMIP6/Models/ISMIP6Greenland_Thermal.mat');
   temperature = md3dT.results.ThermalSolution.Temperature;
   rheology_B_temp = cuffey(temperature);
   md.materials.rheology_B = InterpFromModel3dToMesh2d(md3dT,rheology_B_temp,md.mesh.x,md.mesh.y,NaN,258);
   md.materials.rheology_n = 3*ones(md.mesh.numberofelements,1);
   md=mechanicalproperties(md,velx,vely);
   maxEffStrain = .5;
   shear_margin_pos=md.mesh.elements(md.results.strainrate.effectivevalue>maxEffStrain);

   disp('   -- Finalizing the mesh');
   hVertices = NaN(md.mesh.numberofvertices,1);
   hminVertices = NaN(md.mesh.numberofvertices,1);
   distance2front = reinitializelevelset(md, ice_levelset); %(M==0)-.5);
   hVertices(find(distance2front<3e3 & distance2front>-5e3))=400;  % 3 km downstream, 5 km upstream
   hVertices(shear_margin_pos)=400;
	%hVertices(floating_ice_mask==1)=400;
	flags=ContourToNodes(md.mesh.x,md.mesh.y,'Exp/NoRefinementGreenland.exp',2);
   hminVertices(find(flags==1))=1000;
   md=bamg(md,'gradation',1.6,'KeepVertices',1,'hmax',hmax,'hmin',hmin,'hVertices',hVertices,'hminVertices',hminVertices,'field',vel,'err',4);

   [md.mesh.lat,md.mesh.long]=xy2ll(md.mesh.x,md.mesh.y,-1);
   md.mesh.epsg=3031;
   md.mesh.scale_factor=(1+sin(md.mesh.lat*pi/180))/(1+sin(-71*pi/180));
   md.miscellaneous.name='mesh';

   savemodel(org,md);
end %}}}
if perform(org,[domain_name '_Param']),% {{{

   md=loadmodel(org,[domain_name '_Mesh']);

   disp('   Setting up materials');
   disp('   -- Densities ');
   md.materials.rho_ice       = 917.;  % (917 is the value used in BedMachine)
   md.materials.rho_water     = 1027.; % ocean water (1027 is the value used in BedMachine)
   md.materials.rho_freshwater= 1000.; % fresh water
   md.constants.g             = 9.81;  % gravitational acceleration
  
   disp('   Setting up geometry');
   %Minimum values used in the parameterzation
   min_thickness                 = 10.; % minimum ice thickness used to setup the initial geometry
   min_surface                   = min_thickness*(1-md.materials.rho_ice/md.materials.rho_water) ;
   md.masstransport.min_thickness= min_thickness;

   %md.geometry.surface              = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'surface');
   md.geometry.surface              = interpGimpdem(md.mesh.x, md.mesh.y) - interpBedmachineGreenland(md.mesh.x,md.mesh.y,'geoid');
   pos                              = find(md.geometry.surface<1.e-10); % ocean part or rocks
   md.geometry.surface(pos)         = min_surface; % set minimum ice surface on the ocean part
   md.geometry.bed                  = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'bed'); % interp method: linear
   md.geometry.base                 = zeros(length(md.geometry.bed),1); % initial setup
   % Setting up ice thickness, ice base and grounded ice level set based on the hydrostatic equilibrium
   floatation_base                  = md.geometry.surface*md.materials.rho_ice/(md.materials.rho_ice-md.materials.rho_water); % using the corrected surface
   pos                              = find(floatation_base<=md.geometry.bed); % grounded ice
   md.geometry.base(pos)            = md.geometry.bed(pos);
   %md.mask.ocean_levelset(pos)      = 1;
   pos                              = find(floatation_base>md.geometry.bed); % floating ice
   md.geometry.base(pos)            = floatation_base(pos);
   %md.mask.ocean_levelset(pos)      = -1;
   md.geometry.thickness            = md.geometry.surface-md.geometry.base;
   pos                              = find(md.geometry.thickness<min_thickness); % dealing with rocks or ocean part
   md.geometry.thickness(pos)       = min_thickness;
   md.geometry.surface(pos)         = md.geometry.thickness(pos)+md.geometry.base(pos);

  
	disp('   Creating ice mask');
   %disp('Using Greene monthly ice front--->')
   %md.mask.ocean_levelset = zeros(md.mesh.numberofvertices,1);
   % get observed calving front
   ice_masks = interpMonthlyIceMaskGreene(md.mesh.x, md.mesh.y, [calibration_start_time, calibration_start_time+1]);
   ice_mask = reinitializelevelset(md,ice_masks(1:end-1,1));
	disp([' Using Greene ice mask from year ' num2str(ice_masks(end,1))]); 
   % initial levelset
	if use_exp_for_ice_levelset
		md.mask.ice_levelset = updateLevelset(md,ice_mask,use_exp_for_ice_levelset);
	else
		md.mask.ice_levelset = ice_mask;
	end

   %disp('   Adjusting ice mask');
   % Tricky part here: we want to offset the mask by one element so that we don't end up with a cliff at the transition
   % Find the elements in which there is at least one vertex with positive mask
   pos = find(max(md.mask.ice_levelset(md.mesh.elements),[],2)>0);
   md.mask.ice_levelset(md.mesh.elements(pos,:))= 1; % setting no ice
	%Reinitialize
   md.mask.ice_levelset = reinitializelevelset(md, md.mask.ice_levelset);

	disp('   Calculating ocean mask');
   md.mask.ocean_levelset = sethydrostaticmask(md); %calculates height above floatation as is done in C++ code
	% Ensure base is consistent with the bed
	pos = find((md.mask.ocean_levelset>0) & (md.geometry.bed ~= md.geometry.base));
	if length(pos) > 0
		disp([num2str(length(pos)) ' vertices have bases that are inconsistent with the bed.']);
		md.geometry.base(pos) = md.geometry.bed(pos);
		md.geometry.thickness(pos) = md.geometry.surface(pos)-md.geometry.base(pos);
		pos = find(md.geometry.thickness<min_thickness);
      md.geometry.thickness(pos) = min_thickness;
		md.geometry.surface(pos) = md.geometry.thickness(pos)+md.geometry.base(pos);
	end

	% Testing this section to bring back fjords to SE Greenland
	%M = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask','nearest');
   %pos = find((md.mask.ocean_levelset < 0) & (M < 1));
	ice_masks_noRock = interpMonthlyIceMaskGreene(md.mesh.x, md.mesh.y, [calibration_start_time, calibration_start_time+1],0);
	ice_levelset_noRock = ice_masks_noRock(1:end-1,1);
	pos = find((md.mask.ocean_levelset < 0) & (ice_levelset_noRock > 0)); % if there is ocean and should be no ice, then make it ocean.
   md.mask.ice_levelset(pos) = 999999;
	md.mask.ice_levelset = reinitializelevelset(md, md.mask.ice_levelset);

   disp('   Reading and interpolating velocities ');
   %[md.inversion.vx_obs md.inversion.vy_obs] = interpJoughinCompositeGreenland(md.mesh.x,md.mesh.y);
	vel_obs = interpFromMEaSUREsGeotiffs(md.mesh.x,md.mesh.y,[1994],[2016],'product',[670],'output',{'vx','vy'});
	md.inversion.vx_obs = vel_obs.vx;
	md.inversion.vy_obs = vel_obs.vy;
   % Initialize velocities 
   md.initialization.vx=md.inversion.vx_obs;
   md.initialization.vy=md.inversion.vy_obs;
   md.initialization.vel=sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2);
   % Find NaN locations	
   pos=find(isnan(md.inversion.vx_obs) | isnan(md.inversion.vy_obs));
	pos0=find(~isnan(md.inversion.vx_obs) & ~isnan(md.inversion.vy_obs));
   % Interpolate NaN locations
	md.initialization.vx(pos) = griddata(md.mesh.x(pos0),md.mesh.y(pos0),md.inversion.vx_obs(pos0),md.mesh.x(pos),md.mesh.y(pos));
	md.initialization.vy(pos) = griddata(md.mesh.x(pos0),md.mesh.y(pos0),md.inversion.vy_obs(pos0),md.mesh.x(pos),md.mesh.y(pos));
	md.initialization.vx((isnan(md.initialization.vx) | isnan(md.initialization.vy))) = 0;
	md.initialization.vy((isnan(md.initialization.vx) | isnan(md.initialization.vy))) = 0;
	md.initialization.vel(pos) = sqrt(md.initialization.vx(pos).^2+md.initialization.vy(pos).^2);
	% No vertical velocity
   md.initialization.vz  = zeros(md.mesh.numberofvertices,1);
   % Ensure obs for inversion doesn't have NaNs either
	md.inversion.vx_obs = md.initialization.vx;
	md.inversion.vy_obs = md.initialization.vy;
	md.inversion.vel_obs = md.initialization.vel;

   disp('   Setting up thermal regime');
   md3dT = loadmodel('/totten_1/ModelData/ISMIP6/Models/ISMIP6Greenland_Thermal.mat');
	temperature = DepthAverage(md3dT,md3dT.results.ThermalSolution.Temperature);
	md3dT_ice_mask = DepthAverage(md3dT, md3dT.mask.ice_levelset);
	% set Nan for no ice in ISMIP
   temperature(md3dT_ice_mask>=0) = NaN;
   % Interpolate
   T = InterpFromMeshToMesh2d(md3dT.mesh.elements2d, md3dT.mesh.x2d, md3dT.mesh.y2d, temperature, md.mesh.x, md.mesh.y);
	% interpolate average T of fast flowing regions onto the areas with no ice
	vellimit = prctile(md.initialization.vel,99);
	% compute the average temperature over the fast flowing regions
	masked = (md.initialization.vel <= vellimit) | (md.mask.ice_levelset>=0);
	[intData, avgT, areas] = integrateOverDomain(md, T, masked);
	fprintf('    -- The average temperature of %s at fast flowing region (vel>%g) is %g ^oC\n', domain_name, vellimit, avgT-273.15);
   % set no ice region in ISMIP to be the average temperature of the fast flowing retions
   temperature(md3dT_ice_mask>=0) = avgT;
	% interpolate onto model grid
   md.initialization.temperature=InterpFromMeshToMesh2d(md3dT.mesh.elements2d, md3dT.mesh.x2d, md3dT.mesh.y2d, temperature, md.mesh.x, md.mesh.y);
	% OLD CODE
	%pos = find(md.mask.ice_levelset > -1000);
   %md.initialization.temperature(pos) = mean(md.initialization.temperature); %260;
   %md.initialization.temperature      = min(0,interpSeaRISE(md.mesh.x,md.mesh.y,'temp',-1))+273.15;
   md.initialization.waterfraction  = zeros(md.mesh.numberofvertices,1);
   md.initialization.watercolumn    = zeros(md.mesh.numberofvertices,1);
   md.thermal.spctemperature        = md.initialization.temperature;
   md.thermal.isenthalpy            = 1;
   md.thermal.isdynamicbasalspc     = 1;

   disp('   -- Creating flow law parameters');
   md.materials.rheology_n    = 3*ones(md.mesh.numberofelements,1);
   md.materials.rheology_law  = 'Cuffey';
	md.materials.rheology_B    = cuffey(md.initialization.temperature);
   %temperature = md3dT.results.ThermalSolution.Temperature;
   %rheology_B_temp = cuffey(temperature);
   %md.materials.rheology_B = InterpFromModel3dToMesh2d(md3dT,rheology_B_temp,md.mesh.x,md.mesh.y,NaN,258);
	% Only soften the shear margins for certain glaciers
	for ii = soften_shear_margins 
		disp('   -- Soften shear margins for some glaciers'); 
      maxEffStrain = .5;
      if ~exist(['./Exp/basin' num2str(ii) '_domain_' num2str(mesh_buffer) 'mBuffer.exp'],'file')
         file_name = makeDomainBoundary(ii, './Exp', mesh_buffer, ['basin' num2str(ii)])
      end
      flag_nodes = ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,['./Exp/basin' num2str(ii) '_domain_' num2str(mesh_buffer) 'mBuffer.exp'],'node',1);
      pos_flag_nodes = find(flag_nodes==1);
		disp([num2str(length(pos_flag_nodes)) ' flagged nodes found for basin ' num2str(ii)]);
      flag_elements = sum(ismember(md.mesh.elements,pos_flag_nodes),2);
      flag_elements(flag_elements>0) = 1;

      mde = extract(md,flag_elements,'spccheck',0);

	   mde=mechanicalproperties(mde,mde.inversion.vx_obs,mde.inversion.vy_obs);
      pos1 = mde.mesh.elements(mde.results.strainrate.effectivevalue>maxEffStrain);
      pos2 = find(mde.inversion.vel_obs<prctile(mde.initialization.vel,90));
		disp(['shear margin velocity cutoff ' num2str(prctile(mde.initialization.vel,90))]);
      pos = intersect(pos1,pos2);
      damage=ones(mde.mesh.numberofvertices,1);
      damage(pos)=3; %2; %5;
      % exclude no ice elements
      damage(mde.mask.ice_levelset>0) = 1;
      % set damage
      md.materials.rheology_B(mde.mesh.extractedvertices)=mde.materials.rheology_B./damage;
   end
	md.materials.rheology_B=averaging(md,md.materials.rheology_B,2);

   disp('   Loading accumulation rates from RACMO (SMB_MEAN1960-1989_150m.nc)');
   smb = interpRACMO1km(md.mesh.x,md.mesh.y); % ATTENTION: ice density assumed as 917 kg/m^3)
	smb(smb==-9999)=0;
	pos0 = find(smb==0);
	pos = find(smb~=0);
	smb(pos0) = griddata(md.mesh.x(pos),md.mesh.y(pos),smb(pos),md.mesh.x(pos0),md.mesh.y(pos0),'nearest');
   md.smb.mass_balance = smb;
	% What I jused prior to 2/2/2024:
   %ocean_pos = find(md.mask.ice_levelset > -5000 | md.mask.ocean_levelset < 0);
   %md.smb.mass_balance(ocean_pos) = min(md.smb.mass_balance);

   disp('   Initializing basal friction using driving stress');
   disp('   -- Smooth the ice surface with 20 L2 projections and then compute the surface slopes');
   asurf      = averaging(md,md.geometry.surface,20); % maybe executing 20 L2 projection is ok
   athick   = averaging(md,md.geometry.thickness,20);
   [sx,sy,s]= slope(md,asurf); % slope 's' comes on elements
   sslope  = averaging(md,s,1); % average the slope once on the vertices, because 's' comes on elements, we need this data on vertices
   disp('   -- Process surface velocity data');
   vel     = md.inversion.vel_obs;
   flags      = (vel==0).*(md.mask.ice_levelset<0); % interpolate on the ice parts
   pos1    = find(flags);
   pos2    = find(~flags);
   vel(pos1)= griddata(md.mesh.x(pos2),md.mesh.y(pos2),vel(pos2),md.mesh.x(pos1),md.mesh.y(pos1)); % interpolating the velocities where vel==0
   vel     = max(vel, 0.1); % setting minimum velocity value
   disp('   -- Calculate effective pressure and the initial pressure');
   Neff                      = (md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base)*md.constants.g;
   Neff(find(Neff<=0))       = 1.; % setting minimum positve pressure
   md.initialization.pressure   = md.materials.rho_ice*athick*md.constants.g; % setting the initial pressure
   disp('   -- Deduce friction coefficient from driving stress');
   driving_stress            = md.materials.rho_ice*md.constants.g*athick.*(sslope);
   md.friction.coefficient = sqrt(driving_stress./(Neff.*vel/md.constants.yts));
   md.friction.coefficient   = min(md.friction.coefficient,400);
	md.friction.p          = ones(md.mesh.numberofelements,1);
   md.friction.q          = ones(md.mesh.numberofelements,1);
   disp('   -- Extrapolate on ice free regions (using griddata)');
   flags   = (md.mask.ice_levelset>0); % no ice
   pos1 = find(flags);
   pos2 = find(~flags);
   md.friction.coefficient(pos1) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),md.friction.coefficient(pos2),md.mesh.x(pos1),md.mesh.y(pos1));
   pos = find(isnan(md.friction.coefficient) | md.friction.coefficient <=0);
   md.friction.coefficient(pos)  = 50.;

   disp('   Loading melting rates - no basal melt currently');
   md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1); % ATTENTION: no melting on grounded ice for now
   md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices,1);

   disp('   Loading geothermal flux from Shapiro et al.');
   md.basalforcings.geothermalflux  = interpSeaRISE(md.mesh.x,md.mesh.y,'bheatflx');

   disp('   Setting boundary conditions');
   md.stressbalance.spcvx        = NaN(md.mesh.numberofvertices,1);
   md.stressbalance.spcvy        = NaN(md.mesh.numberofvertices,1);
   md.stressbalance.spcvz        = NaN(md.mesh.numberofvertices,1);
   md.stressbalance.referential  = NaN(md.mesh.numberofvertices,6);
   md.stressbalance.loadingforce = zeros(md.mesh.numberofvertices,3);
   md.masstransport.spcthickness = NaN(md.mesh.numberofvertices,1);
   pos = find((md.mask.ice_levelset).*(md.mesh.vertexonboundary)); % the contour part of the ice
   md.stressbalance.spcvx(pos)   = md.inversion.vx_obs(pos);
   md.stressbalance.spcvy(pos)   = md.inversion.vy_obs(pos);
   md.masstransport.spcthickness(pos)  = md.geometry.thickness(pos);
   if use_exp_for_spcvxvy
   	flags=ContourToNodes(md.mesh.x,md.mesh.y,use_exp_for_spcvxvy,2);
		pos = find(flags==1);
      md.stressbalance.spcvx(pos) = md.inversion.vx_obs(pos);
      md.stressbalance.spcvy(pos) = md.inversion.vy_obs(pos);
	end

   md.masstransport.stabilization = 5;

   %Use SSA for now
   md=setflowequation(md,'SSA','all');

   %Save model
   md.miscellaneous.name = 'Param';
   savemodel(org,md);
end%}}}
if perform(org,[domain_name '_InversionC_fric' num2str(fric_type) fric_init]),% {{{
   
	if fric_type > 1
		md=loadmodel(org,[domain_name '_InversionC_fric1cal']);
      disp('loading model from fric_type == 1, InversionC step');
	elseif strcmp(domain_name,'NO') | strcmp(domain_name,'NE') 
      md=loadmodel(org,[domain_name '_InversionB']);
		disp('loading model from InversionB step');
	else
		md=loadmodel(org,[domain_name '_Param']);
		disp('loading model from Param step');
	end
	
	md.cluster=cluster;
	if strcmpi(clustername,'discovery') | strcmpi(clustername,'andes')
	   md.cluster.codepath = path_noAD;
	end

   if use_exp_for_spcvxvy
   	flags=ContourToNodes(md.mesh.x,md.mesh.y,use_exp_for_spcvxvy,2);
		pos = find(flags==1);
      md.stressbalance.spcvx(pos) = md.inversion.vx_obs(pos);
      md.stressbalance.spcvy(pos) = md.inversion.vy_obs(pos);
	end

   %ratio for coeff in inversion
   ratio = [10, 2, 1];

   % set M1QN3 package
   md.inversion=m1qn3inversion(md.inversion);
   md.inversion.iscontrol=1;
   md.inversion.maxsteps=150;
   md.inversion.maxiter=150;
   md.inversion.dxmin=0.01;
   md.inversion.gttol=1.0e-6; % This should be more lik 0.01
   md.inversion.incomplete_adjoint=0; % 0: non linear viscosity, 1: linear viscosity 04/29/2019 changed to non linear
	md.inversion.dfmin_frac=0.3;

   % Find where there is no data and where cost function coefficients should be zero
   vel_obs = interpFromMEaSUREsGeotiffs(md.mesh.x,md.mesh.y,[1994],[2016],'product',[670],'output',{'vx','vy'});
	vx_obs = vel_obs.vx;
	vy_obs = vel_obs.vy;
   %[vx_obs vy_obs] = interpJoughinCompositeGreenland(md.mesh.x,md.mesh.y);
	mask_obs=ones(md.mesh.numberofvertices,1);
   mask_obs(isnan(vx_obs) & isnan(vy_obs))=0;
	%Ensure all elements that contained a NaN in the obs are fully excluded
   pos=find(min(mask_obs(md.mesh.elements)')==0);
   mask_obs(md.mesh.elements(pos,:))=0;
   mask_obs = find(mask_obs==0);

   % Friction law and cost functions definition
   if fric_type==1 % Budd - linear
      disp('   -- Setting up a Budd''s sliding law (with p=1, q=1)');
      md.friction.p=ones(md.mesh.numberofelements,1);
      md.friction.q=ones(md.mesh.numberofelements,1);
      md.inversion.cost_functions=[101 103 501];
      md.inversion.cost_functions_coefficients = ones(md.mesh.numberofvertices,3); % IMPORTANT: set as 0 again the vertices with no data (NaN)
      cost_function_coefficients = [3000, 0.52, 1.6e-8]; %[3000, 0.78, 4.77e-08];
      md.miscellaneous.dummy.cost_function_coefficients = cost_function_coefficients;
      md.inversion.cost_functions_coefficients(:,1)=cost_function_coefficients(1);
      md.inversion.cost_functions_coefficients(:,2)=cost_function_coefficients(2);
      md.inversion.cost_functions_coefficients(:,3)=cost_function_coefficients(3);
      md.inversion.cost_functions_coefficients(mask_obs,1:2)=0; % positions with NaN in the velocity data set
      pos_noice=find(md.mask.ocean_levelset<0); %find where there is no grounded ice
      md.inversion.cost_functions_coefficients(pos_noice,1:2)=0; % positions with no grounded ice

      %Initial guess
      md.friction.coefficient=EstimateFric_Budd(md,1,1); % initial guess from Driving Stress (using r=1 and s=1)
      % Extrapolate friction to ocean based on bed elevation relationship
      C = md.friction.coefficient;
      pos=find(C==0 | isnan(C) | md.mask.ocean_levelset<=0); 
		B = md.geometry.bed(pos);
		C(pos) = ((B+1500).^2)/40000 + (B/1000);
      md.friction.coefficient = C;
      
      % Controls
	   md.inversion.control_parameters={'FrictionCoefficient'};
      md.inversion.min_parameters=0.1*ones(md.mesh.numberofvertices,1);
      md.inversion.max_parameters=400*ones(md.mesh.numberofvertices,1);

		%if use_exp_for_snapshot_inversion
      %   flags=ContourToNodes(md.mesh.x,md.mesh.y,use_exp_for_snapshot_inversion,2);
      %   pos = find(flags==1);
      %   md.inversion.min_parameters(pos) = md.friction.coefficient(pos);
      %   md.inversion.max_parameters(pos) = md.friction.coefficient(pos);
		%end

	elseif fric_type==2 % Schoof
	   disp('   -- Setting up a Schoof''s sliding law (with m=1/3, Cmax=0.5)');
		md=SwitchFric_LinearBudd_to_Schoof(md); % initial guess from previous inversion

		if strcmp(fric_init,'Budd')
		   md.inversion.iscontrol=0;
		
      elseif strcmp(fric_init,'cal')
         md.inversion.cost_functions=[101 103 501];
         md.inversion.cost_functions_coefficients = ones(md.mesh.numberofvertices,3); % IMPORTANT: set as 0 again the vertices with no data (NaN)
         cost_function_coefficients = [3000, 0.52*10, (1.6e-8)/250]; %/25]; %/1650];
         md.miscellaneous.dummy.cost_function_coefficients = cost_function_coefficients;
         md.inversion.cost_functions_coefficients(:,1)=cost_function_coefficients(1);
         md.inversion.cost_functions_coefficients(:,2)=cost_function_coefficients(2);
         md.inversion.cost_functions_coefficients(:,3)=cost_function_coefficients(3);
         md.inversion.cost_functions_coefficients(mask_obs,1:2)=0; % positions with NaN in the velocity data set
         pos_noice=find(md.mask.ocean_levelset<0); %find where there is no grounded ice
         md.inversion.cost_functions_coefficients(pos_noice,1:2)=0; % positions with no grounded ice

         %Controls
         md.inversion.control_parameters={'FrictionC'};
         md.inversion.min_parameters=0.1*ones(md.mesh.numberofvertices,1); 
			pos = find((md.geometry.thickness==md.masstransport.min_thickness)|(md.mask.ocean_levelset<0));
			md.inversion.min_parameters(pos) = 10;
         md.inversion.max_parameters=12000*ones(md.mesh.numberofvertices,1); %66000*ones(md.mesh.numberofvertices,1);% 10000*ones(md.mesh.numberofvertices,1); %
		else
			error('fric_init not supported');
		end

   	pos = find(((md.geometry.thickness==md.masstransport.min_thickness)|(md.mask.ocean_levelset<0)) & md.friction.C<10);
		md.friction.C(pos) = 10;

	   pos = find(((md.geometry.thickness==md.masstransport.min_thickness)|(md.mask.ocean_levelset<0)) & md.friction.C<10);
		md.friction.C(pos) = 10;

	elseif fric_type==3 % Budd p = 6, q = 6
	   disp('   -- Setting up a Budd''s sliding law (with p=6, q=6)');
		md=SwitchFric_LinearBudd_to_ExponentialBudd(md,6,6); % initial guess from previous inversion

		if strcmp(fric_init,'Budd')
		   md.inversion.iscontrol=0;
		
      elseif strcmp(fric_init,'cal')
         md.inversion.cost_functions=[101 103 501];
         md.inversion.cost_functions_coefficients = ones(md.mesh.numberofvertices,3); % IMPORTANT: set as 0 again the vertices with no data (NaN)
         cost_function_coefficients = [3000, 0.52*30, (1.6e-8)*120000]; %/25]; %/1650];
         md.miscellaneous.dummy.cost_function_coefficients = cost_function_coefficients;
         md.inversion.cost_functions_coefficients(:,1)=cost_function_coefficients(1);
         md.inversion.cost_functions_coefficients(:,2)=cost_function_coefficients(2);
         md.inversion.cost_functions_coefficients(:,3)=cost_function_coefficients(3);
         md.inversion.cost_functions_coefficients(mask_obs,1:2)=0; % positions with NaN in the velocity data set
         pos_noice=find(md.mask.ocean_levelset<0); %find where there is no grounded ice
         md.inversion.cost_functions_coefficients(pos_noice,1:2)=0; % positions with no grounded ice

         %Controls
         md.inversion.control_parameters={'FrictionCoefficient'};
         md.inversion.min_parameters=0.1*ones(md.mesh.numberofvertices,1); 
         md.inversion.max_parameters=10*ones(md.mesh.numberofvertices,1);
		else
			error('fric_init not supported');
		end

	else
      error('fric_type not supported');
   end

   % Additional parameter for stress balance
   %md.stressbalance.restol=0.0001; % 04/29/2019
   md.stressbalance.restol=0.0001; % 04/29/2019
   md.stressbalance.reltol=0.001; % 04/29/2019
   md.stressbalance.abstol=NaN; % 04/29/2019
   md.stressbalance.maxiter=150; % 04/29/2019

   % Set cluster
   md.cluster.interactive = interactive;
   md.settings.waitonlock = waitonlock;
   md.verbose=verbose('solution',false,'control',true);
   md.miscellaneous.name= org.steps(org.currentstep).string;

   % Define solver
	if ~contains(md.cluster.codepath,'ISSM_AD') %hs changed path from ISSM-AD 1/29
      md.toolkits.DefaultAnalysis=bcgslbjacobioptions();% biconjugate gradient with block Jacobi preconditioner
	end
   md.settings.solver_residue_threshold=NaN; % 11/05/2019

   md=solve(md,'Stressbalance','runtimename',0,'loadonly',loadonly);
	
	if interactive==1
      md=loadresultsfromcluster(md);
      try
         newCoeff = recommendedCostCoeff(md, ratio);
         disp(sprintf('With the given ratio %d, %d, %d \n', ratio))
         disp(sprintf('The coefficients can be updated to %g,   %g,   %g \n', newCoeff))
      	disp(sprintf('From %g,   %g,   %g \n', cost_function_coefficients))
      	  catch
         disp('Cannot recommend cost function coefficients.');
      end
	end
  
   if (loadonly==1) | (interactive==1)
      % Update friction coefficient and velocity field
		if strcmp(fric_init,'cal')
         if fric_type==2
            md.friction.C=md.results.StressbalanceSolution.FrictionC;
		   else
            md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;
			end
      end
      md.initialization.vx=md.results.StressbalanceSolution.Vx;
      md.initialization.vy=md.results.StressbalanceSolution.Vy;
      md.initialization.vel=md.results.StressbalanceSolution.Vel;

      savemodel(org,md);
	end
end%}}}
