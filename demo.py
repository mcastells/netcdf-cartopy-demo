import MFT
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import math
import matplotlib.pyplot as plt
import cartopy as cart
import cartopy.crs as ccrs

ds = Dataset('wrf_output.nc', 'r')

print('ds.data_model = %s' % ds.data_model)

print('ds.groups = {0}'.format(ds.groups))

#print(f'ds.dimensions = {ds.dimensions}')

print('\nds.dimensions - ')
for dim, obj in ds.dimensions.items():
	print(f'{dim}: size = {obj.size}')

#print(f'ds.variables = {ds.variables}')

print('\nds.variables = ')
for var in ds.variables:
	print(var)

latdim = ds.dimensions['lat']
londim = ds.dimensions['lon']
nlat = len(latdim)
nlon = len(londim)

print(f'{latdim = }')

lat = ds['lat']
lon = ds['lon']
u10 = ds['U10']
v10 = ds['V10']
uc = ds['UC']
vc = ds['VC']
slp = ds['SLP']
t2 = ds['T2']
q2 = ds['Q2']
sst = ds['SST']

print(f'{u10[:] = }')

print(f'{np.ma.count(u10[:]) = }')
print(f'{np.ma.count(v10[:]) = }')

windspeeds = np.ma.empty(shape=(nlat,nlon))
curl = np.ma.empty_like(windspeeds)

u10mask = ma.getmaskarray(u10[:])
u10_masked_indices = ma.where(u10mask == True)
u10_unmasked_indices = ma.where(u10mask == False)

two_array_wind = ma.array(np.full(shape=(nlat, nlon), fill_value=2.0, dtype=float))
two_array_wind[u10_masked_indices] = ma.masked

print('Calculating windspeeds...')

windspeeds = ma.sqrt(ma.add(ma.power(u10,two_array_wind), ma.power(v10,two_array_wind)))

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
ax1.coastlines()

ax1.contourf(lon[:], lat[:], windspeeds[:], cmap='cool')
ax1.quiver(lon[::2], lat[::2], u10[::2,::2], v10[::2,::2], scale=250)

CONV_CRIT = 0.00005     #convergence critereon (fractional change)  []  
CONVECT = 0.0          #convective parameter  
warn = 1                #warning are given     
eqv_neut = 0            #output winds are winds rather than equivalent neutral winds  
z_wanted = 10.0         #height to which winds, potential temp, and humidity are adjusted                                
flux_model = 9          #BVW model=0  
Qnet = 5.0
sst_prm = 0
z0_mom_prm = 0
z0_TQ_prm = 0
stable_prm = 0
wind_ang = 0
wave_ang = 0

dyn_in_prm = 0 # Wind speed, relative to the surface current, m/s
dyn_in_val2 = 0.0
ref_ht_wind = 10.0 # Height of the wind observations, m
CONVECT = 0 # Convective parameter. Recommended value between 0.7 and 1.25. For details see TOGA NOTES #4 (recommendation comes from no capillary waves)
air_moist_prm = 0 # Relative humidity, specific humidity
sfc_moist_prm = 1 # Relative humidity, fraction
sfc_moist_val = 0.98 # 98% assumed because of salinity obstructing evaporation
salinity = 34.9 / 1000.0 # Salinity, fraction, No salinity in dataset (global average)
ss_prm = 0 # Sea state parameterization, wind-wave stability parameter
ss_val = 1.0 # No sea state data in dataset, set to 1.0 for local equilibrium
ref_ht_tq = 2 # Height of temperature and humidity observations, m
sst_prm = 0 # Designates surface temperature as skin temperature?
astab = 1 # Atmospheric stability is calculated

stress = np.ma.masked_all(shape=(nlat,nlon))

for (x,y), s in np.ndenumerate(stress):
	index = (x,y)

	if sst[index] is not ma.masked:
		dyn_in_val = float(windspeeds[index])
		pressure = float(slp[index]) * 100.0 # Atmospheric surface pressure, Pa (converted from hPa in dataset)
		air_moist_val = float(q2[index]) # (converted from g/kg to fraction)
		t_air = float(t2[index]) - 273.15 # Air temperature at the reference height of the thermometer and humidity sensor, C
		t_skin = float(sst[index]) - 273.15 # skin tempearture, C
		
		pass_by_ref_params = {'shf':0, 'lhf':0, 'tau':[0,0], 'u_star':[0,0], 't_star':0, 'q_star':0, 'z_over_L':0, 'wave_age':0, 'dom_phs_spd':0,
			'h_sig':0, 'ww_stab':0, 'zo_m':[0, 0], 'u_at_z':0, 't_at_z':0, 'q_at_z':0}

		count = MFT.ht_adj_( dyn_in_prm, dyn_in_val, dyn_in_val2, CONVECT, CONV_CRIT,
			pressure, air_moist_prm, air_moist_val, sfc_moist_prm, sfc_moist_val,
			salinity, ss_prm, ss_val, t_air, sst_prm, t_skin, ref_ht_wind, ref_ht_tq,
			z_wanted, astab, eqv_neut, Qnet, warn, flux_model, z0_mom_prm, z0_TQ_prm, stable_prm,
			pass_by_ref_params )
		if count > 0:
			stress[index] = pass_by_ref_params["tau"][0]
			#print(f'{pass_by_ref_params["tau"][0] = }')

print(f'{stress = }')

ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
ax2.coastlines()
ax2.add_feature(cart.feature.LAND, zorder=100, edgecolor='k', facecolor='#222222')

ax2.contourf(lon[:], lat[:], stress[:])

plt.show()

ds.close()
