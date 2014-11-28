% Copyright 2010 Jakob Heyman
% This code is free software: you can redistribute it and/or modify it under the terms of the GNU
% General Public License; version 3 of the License or any later version. There is no warranty for
% this code. For a copy of the GNU Genaral Public License, see <http://www.gnu.org/licenses/>.

% Boulder exposure age script used in Heyman et al. (2010):
% Too young or too old: evaluating cosmogenic exposure dating based on an analysis of compiled
% boulder exposure ages. Earth and Planetary Science Letters.

% Monte Carlo simulation of multiple-boulder group (random deglaciation age) apparent exposure ages
% assuming prior exposure and incomplete post-depositional exposure. Each boulder group has a random
% duration of prior exposure and each individual boulder has a random prior sample depth. Each
% boulder is exhumed through and shielding by till (random depth at deglaciation) with a time-
% dependent exponential exhumation rate dx/dt=ae^(-bt). Be-10 surface prod rate due to spallation is
% taken as 98.8% of the total surface prod rate and Be-10 surface prod rate due to muon interaction
% is taken as 1.2% of the total surface prod rate. Depth-dependent spallogenic Be-10 production is
% based on Lal (1991), and muogenic Be-10 production is based on Granger and Smith (2000).

% Contact info: jakob.heyman@natgeo.su.se

clear all;

boulders = [2,3,4,5,6,7,8,9,10,11,12,14];		% nr of boulders per discrete glacial deposit
samplnr = [84,90,59,51,26,10,10,4,2,4,1,1];		% nr of boulder groups with [2,3,4 ... 14] boulders
samplnr = samplnr .* 20;						% nr of Monte Carlo iterations
maxpriordur = 50000;	% maximum duration of prior exposure (yr)
max_pri_d = 500;		% maximum prior sample depth (cm)
maxage = 450000;		% maximum deglaciation age (yr)
tstep = 50;				% time-step for boulder exhumation and Be-10 production calculation (yr)
aa = 0.011;				% time-dep exhumation curve coefficient (cm/yr)
bb = 0.000010;			% time-dep exhumation curve coefficient (yr^-1)
tdens = 2;				% till density (g/cm^3)
rdens = 2.7;			% rock density (g/cm^3)
atten = 160;			% spallogenic prod. effective attenuation length (g/cm^2)
rdc = -log(0.5)/1387000;% radioactive decay constant (yr^-1) (Chmeleff et al. 2010; Korchinek et al. 2010)
spfract = 0.988;		% spallogenic fraction of total surface production rate
mufract = 0.012;		% muogenic fraction of total surface production rate
m1 = 0.76;				% muogenic prod. approxim. coefficient
m2 = 0.11;				% muogenic prod. approxim. coefficient
m3 = 0.13;				% muogenic prod. approxim. coefficient
L1 = 738.6;				% muogenic prod. approxim. effect. atten. length (g/cm^2)
L2 = 2688;				% muogenic prod. approxim. effect. atten. length (g/cm^2)
L3 = 4360;				% muogenic prod. approxim. effect. atten. length (g/cm^2)

timev = (tstep:tstep:maxage)';				% vector for time-stepping (time steps)
decay = -exp(-rdc .* tstep);				% Be-10 decay per time step

for i = 1:numel(boulders)		% loop for boulder group data
	priordur = rand (1,samplnr (i)) .* maxpriordur;					% random duration of prior exposure
	deglacage = rand (1,samplnr (i)) .* maxage;						% random deglaciation-age for boulder group
	maxdepth = aa./bb .* (1 - exp(-bb .* deglacage));				% maximum depth for exhumation based on deglacage
	timem = repmat (timev,1,samplnr (i));							% matrix for time-stepping (time steps)
	
	for j = 1:boulders(i)			% loop for individual boulder data
		startdepth = rand (1, samplnr (i)) .* maxdepth;				% random boulder deglaciation depth based on maxdepth
		startdepthm = repmat (startdepth,(maxage/tstep),1);			% matrix for time-stepping (start depth)
		depth = startdepthm - aa./bb .* (1 - exp(-bb .* timem));	% boulder depth in time steps
		beprod = tstep .* spfract .* exp(-tdens/atten .* depth) + tstep .* mufract .* ...	% Be-10 production per time step
		(m1 .* exp(-tdens/L1 .* depth) + m2 .* exp(-tdens/L2 .* depth) + m3 .* exp(-tdens/L3 .* depth));
		beconcm = filter (1,[1,decay],beprod);						% calculate Be-10 conc through the exhumation
		exhtime = max ((depth > 0) .* timem);						% boulder exhumation duration
		beconc = max((depth > 0) .* beconcm);						% erase Be-10 change above surface and pick out conc
		pri_d = rand (1,samplnr (i)) .* max_pri_d;					% random prior sample depth
		inhconc = (spfract .* exp(-rdens/atten .* pri_d) + mufract .* (m1 .* exp(-rdens/L1 .* pri_d) + ...
		m2 .* exp(-rdens/L2 .* pri_d) + m3 .* exp(-rdens/L3 .* pri_d))) ./ rdc .* (1 - exp(-rdc .* priordur));
		inhconc = inhconc .* exp(-rdc .* exhtime);					% decay of prior exposure component during exhumation
		expage (j,:) = log (1 - (beconc + inhconc) .* rdc) ./ -rdc + deglacage - exhtime;	% apparent exposure age
		
		clear startdepth;
		clear startdepthm;
		clear depth;
		clear beprod;
		clear beconcm;
		clear exhtime;
		clear beconc;
		clear pri_d;
		clear inhconc;
    end
	
	minbould = min (expage);		% pick out boulder group min age
	meanbould = mean (expage);		% pick out boulder group mean age
	maxbould = max (expage);		% pick out boulder group max age
	standav = std (expage);			% pick out boulder group standard deviation
	
	% fill result matrix with group exposure age data
	result (sum (samplnr (1:i)) - samplnr (1,i) +1:sum (samplnr (1:i)), 1) = deglacage';	% deglac-age (yr) in column 1
	result (sum (samplnr (1:i)) - samplnr (1,i) +1:sum (samplnr (1:i)), 2) = minbould';		% min age (yr) in column 2
	result (sum (samplnr (1:i)) - samplnr (1,i) +1:sum (samplnr (1:i)), 3) = meanbould';	% mean age (yr) in column 3
	result (sum (samplnr (1:i)) - samplnr (1,i) +1:sum (samplnr (1:i)), 4) = maxbould';		% max age (yr) in column 4
	result (sum (samplnr (1:i)) - samplnr (1,i) +1:sum (samplnr (1:i)), 5) = standav';		% stdev (yr) in column 5
	
	clear expage;
	clear deglacage;
	clear maxdepth;
	clear timem;
	clear minbould;
	clear meanbould;
	clear maxbould;
	clear standav;
	clear priordur;
end

save prior_incomplete result -ascii;			% save result as ascii file

beep;