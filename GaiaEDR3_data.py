#!/usr/bin/env python3

'''
Add a description of the code later...
Author: Raphael A. P. Oliveira (https://rapoliveira.github.io)

Gaia EDR3 docs: https://www.cosmos.esa.int/web/gaia/early-data-release-3
'''

### TO-DO LIST ###
# 1) make one flag for each criterium (Francisco)				OK!
# 2) parallax correction before calculating J2000 coordinates	OK!
###

# Standard library imports
import os, os.path, time, sys, re, warnings

# Third-party imports
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, Distance
from astropy.table import Table, Column, hstack
from astropy.utils.exceptions import AstropyWarning
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import matplotlib.pyplot as plt
import numpy as np
import inquirer as iq
import pandas as pd

def main():

	# INPUT (object name or coordinate):
	ngc = "HW79"  # "M67" or "NGC188" (Mateus)
	# ngc = "01:22:48 -75:00:06"  # units: [hms, dms]

	path = os.path.dirname(os.path.realpath(__file__))
	# source_file = os.path.basename(__file__)
	
	print()
	# Change credentials.......!!!!! --> PAREI AQUI!!!
	Gaia.login(credentials_file=f'{path}/my_credentials.txt')
	ngc, c1, rad, out, outname, lind20 = initial_inputs(ngc)  # False
	gaia3, gaia3_valid = astroquery_adql(ngc, c1, rad, outname)
	# success = upload_to_gaia(gaia3, 'tablename', outname)
	gaia3_corr = gaia_codes(gaia3, gaia3_valid, lind20, outname)

def initial_inputs(ngc='', rad=0):

	path = os.path.dirname(os.path.realpath(__file__))
	sys.path.insert(0, os.path.join(path))
	sys.path.insert(0, os.path.join(path, '..'))
	if ngc:
		if ':' not in ngc:
			t = Simbad.query_object(ngc)
			c1 = SkyCoord(t['RA'][0], t['DEC'][0], unit=(u.hourangle, u.deg))
		else:
			c1 = SkyCoord(ngc, unit=(u.hourangle, u.deg))
		RA_DEC = c1.to_string('hmsdms')
		print(ngc, c1, RA_DEC)
	# else:
	# 	import OGLE_RRLyr as ogle	# local import
	# 	ngc, c1, RA_DEC, new_gc = ogle.choose_cluster()
	# 	print(ngc, c1, RA_DEC, new_gc)

	valid_exp = r'^\d*(?:\.(\d|\d\d|\d\d\d))?$'
	if not rad:
		rad = iq.prompt([iq.Text('r', message='Insert the radius (deg)?',validate=\
						lambda _,x: (x != '' and re.match(valid_exp, x)))])['r']
	else: rad = rad/60
	msg = f'{ngc} ({c1.ra.deg:.5f}, {c1.dec.deg:.5f}), rad = {float(rad):.3f}'+\
			' deg. Confirm?'
	conf = iq.prompt({iq.Confirm('conf', message=msg, default=True)})['conf']
	cur = iq.prompt({iq.Confirm('o',message='Save output in current directory?',
				default=False)})['o']
	rad_arcmin = str(float(rad)*60).replace('.','p')
	outname = f'output/{ngc}_gaia3_c{rad_arcmin}arcmin.fits'
	if not cur:
		print(f'- Current directory: {os.getcwd()}\n')
		message = f'Write the directory for the output [Sergio_Bica/{ngc}]'
		valid_dir = lambda _,x: (x == '' or os.path.isdir(f'/{x}/'))
		while True:
			outdir = iq.prompt([iq.Text('dir', message=message, validate=\
								valid_dir)])['dir']
			if not outdir: outdir = f'/Users/rapoliveira/PhD/Sergio_Bica/{ngc}'
			conf = iq.prompt({iq.Confirm('c', message=f'Confirm: \"{outdir}\"?',\
							default=True)})['c']
			if conf: break
		outname = os.path.join(outdir, outname)
	
	lind20 = iq.prompt({iq.Confirm('l',message='Show Lindegren+2020 meshgrid'+
					' figures?', default=False)})['l']

	return (ngc, c1, float(rad), True, outname, lind20)

def astroquery_adql(ngc, c1, rad, outname, out=True):
	"""tap_astroquery(NGC, c1, rad, out)

	Based on:
	- https://astroquery.readthedocs.io/en/latest/gaia/gaia.html
	- https://cosmos.esa.int/web/gaia-users/archive/programmatic-access#PythonFlag

	Args:
		ngc ([type]): [description]
		c1 ([type]): [description]
		rad ([type]): [description]
		out ([type]): [description]
	"""

	Gaia.ROW_LIMIT = -1
	#Gaia.MAIN_GAIA_TABLE = 'gaiaedr3.gaia_source' # -> wrong....

	# 1.1) Rectangular/square projection of the sky -> Different from Vizier!
	#width = u.Quantity(5, u.si.arcmin)
	#height = u.Quantity(5, u.si.arcmin)
	#r = Gaia.query_object_async(coordinate=c1, width=width, height=height)
	#r.pprint()

	# 1.2) Cone search: circle centered at the coord -> Different from Vizier!
	# radius = u.Quantity(rad, u.deg)
	# j = Gaia.cone_search_async(coordinate=c1, radius=radius)
	# r = j.get_results()
	# r.pprint()

	# 1.3) Public tables metadata
	tables = Gaia.load_tables(only_names=True)
	#for table in (tables):
	#    print(table.get_qualified_name())
	gaiaedr3_info = Gaia.load_table('gaiaedr3.gaia_source')
	#gaiadr2_info = Gaia.load_table('gaiadr2.gaia_source')
	print(f"\ngaiaedr3_info = {gaiaedr3_info}\n")
	gaiaedr3_colnames, gaiadr2_colnames = [], []
	for column in (gaiaedr3_info.columns):
		gaiaedr3_colnames.append(column.name)
	#for column in (gaiadr2_info.columns):
	#	gaiadr2_colnames.append(column.name)

	# 1.4 e 1.5) Synchronous query --> arrumar centro!!!
	# job = Gaia.launch_job("select top 100 "
	# 				"solution_id,ref_epoch,ra_dec_corr,astrometric_n_obs_al, "
	# 				"matched_observations,duplicated_source,phot_variable_flag "
	# 				"from gaiaedr3.gaia_source order by source_id")

	def fill_cols(tbl, fill=np.nan, kind='f'):
		"""
		In-place fill of ``tbl`` columns which have dtype ``kind``
		with ``fill`` value.
		"""
		for col in tbl.itercols():
			if col.dtype.kind == kind and np.ma.isMaskedArray(col):
				# col[...] = col.filled(fill)
				tbl[col.name] = col.filled(fill)  # OK!
				# col.fill_value = fill

	# 1.6) Asynchronous query --> arrumar centro!!!
	# query examples: https://www.cosmos.esa.int/web/gaia-users/archive/writing-queries
	# use cases: https://www.cosmos.esa.int/web/gaia-users/archive/use-cases#ClusterAnalysisPythonTutorial
	# https://www.cosmos.esa.int/web/gaia-users/archive/extract-data#ADQLsyntaxTutorial
	point = "POINT('ICRS',gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec)"
	circle = f"CIRCLE('ICRS',{c1.ra.deg:.5f},{c1.dec.deg:.5f},{rad:.3f})"
	adql = f'''SELECT *
		FROM gaiaedr3.gaia_source
			WHERE 1=CONTAINS({point},{circle}))'''
	# adql = f'''SELECT *
	# 	FROM gaiaedr3.gaia_source
	# 		WHERE 1=CONTAINS({point},{circle})
	# 		AND 0=CONTAINS({point},0.01))'''
	job = Gaia.launch_job_async(adql, dump_to_file=out, output_format='fits',
	 							output_file=outname)
	gaia3 = job.get_results()	# type(r) = astropy.Table
	fill_cols(gaia3, fill=np.nan, kind='f')		### IMPORTANT!!!
	fill_cols(gaia3, fill=-99, kind='i')
	gaia3 = Table(gaia3, masked=False)

	#.Retrieving the equatorial coordinates RAJ2000 and DEJ2000 from Vizier (X)
	#..(test later with table.join, which goes wrong)
	# from astroquery.vizier import Vizier
	# v1 = Vizier(columns=['Source','RAJ2000','e_RAJ2000','DEJ2000','e_DEJ2000'],\
	# 			row_limit=1000000)
	# tab = v1.query_region(c1, radius=1.01*rad*u.deg, catalog='I/350/gaiaedr3')
	# wrong_order = list(tab[0]['Source'])
	# j2000 = tab[0][:len(gaia3)].copy()
	# j2000.remove_column('Source')
	# for (i, source) in enumerate(gaia3['source_id']):
	# 	if source in wrong_order:
	# 		idx = wrong_order.index(source)
	# 		j2000[i] = [tab[0]['RAJ2000'][idx], tab[0]['e_RAJ2000'][idx],\
	# 				   tab[0]['DEJ2000'][idx], tab[0]['e_DEJ2000'][idx]]
	# 	else: print('Not found:', i, source)
	# gaia3 = hstack([gaia3, j2000])

	#.Epoch transformation using the proper motions (two options) OK!
	# from pygaia.astrometry.coordinates import EpochPropagation # harder
	from erfa import ErfaWarning
	warnings.simplefilter("ignore", category= ErfaWarning)
	gaia3_valid = gaia3[gaia3['astrometric_params_solved'] > 3]
	t0 = Time(2016.0, format='jyear')
	parallax = np.array(gaia3_valid['parallax'])
	c = SkyCoord(gaia3_valid['ra'], gaia3_valid['dec'], frame="icrs", obstime=t0,
				pm_ra_cosdec=gaia3_valid['pmra'], pm_dec=gaia3_valid['pmdec'], 
				distance=Distance(parallax=np.abs(parallax)*u.mas))	# done!
	t1 = Time('J2000', format='jyear_str')
	new = c.apply_space_motion(new_obstime=t1)
	raj2000, dej2000, i = gaia3['ra'].copy(), gaia3['dec'].copy(), 0
	# diff = np.ones_like(gaia3_valid['ra'], dtype=float)
	for (j, star) in enumerate(gaia3['parallax']):
		if not np.isnan(star):
			raj2000[j], dej2000[j] = new.ra[i].value, new.dec[i].value
			# diff[i] = new.ra[i].value - gaia3['RAJ2000'][j] # < e-10
			i += 1
	raj2000.name, dej2000.name = 'RAJ2000', 'DEJ2000'
	raj2000.description = 'Barycentric right ascension (ICRS) at Ep=2000.0,'+\
						' calculated with Astropy\'s apply_space_motion()'
	dej2000.description = 'Barycentric declination (ICRS) at Ep=2000.0,'+\
						' calculated with Astropy\'s apply_space_motion()'
	gaia3.add_columns([raj2000, dej2000])

	# decode bytes columns to utf-8; transform >i2 to int16
	for col in gaia3.colnames:
		if gaia3[col].dtype == np.object:
			gaia3[col] = gaia3[col].apply(lambda x: x.decode('utf-8'))

	fig = plt.figure(figsize=(6.5,6))#, constrained_layout=True)
	plt.scatter(gaia3['ra'], gaia3['dec'], s=0.8, marker='.')
	plt.gca().update({'xlabel': 'RA (deg)', 'ylabel': 'DEC (deg)'})
	plt.scatter(c1.ra.deg, c1.dec.deg, s=150, marker='+', color='black', zorder=2)
	plt.gca().invert_xaxis()
	plt.tight_layout()
	plt.show()

	catalog = SkyCoord(gaia3['ra'], gaia3['dec'], unit=(u.si.deg, u.si.deg))
	sep = c1.separation(catalog)
	pmra, pmdec = gaia3['pmra'], gaia3['pmdec']
	fig = plt.figure(figsize=(6,6), constrained_layout=True)
	plt.scatter(pmra, pmdec, s=1, marker='.', label=f'r < {rad*60} arcmin', color='#CECECE')
	plt.scatter(pmra[sep < 0.05*u.deg], pmdec[sep < 0.05*u.deg], s=1, marker='.',\
				label='r < 3 arcmin')
	plt.gca().update({'xlim': [-25, 25], 'ylim': [-25, 25]})
	plt.gca().update({'xlabel': 'pmRA (mas/yr)', 'ylabel': 'pmDEC (mas/yr)'})
	plt.legend(loc='best')
	plt.show()

	return gaia3, gaia3_valid

def upload_to_gaia(gaia3, name, filename=''):
	"""upload_to_gaia(tab,name) -> upload astropy Table to Gaia user archive

	Args:
		tab (Table): [description]
		name (str): [description]

	Returns:
		[bool]: if saving is successful
	"""
	#Gaia.login_gui()
	#Gaia.login(user='userName', password='userPassword')
	path = os.path.dirname(os.path.realpath(__file__))
	Gaia.login(credentials_file=f'{path}/credentials.txt')

	#tables = Gaia.load_tables(only_names=True, include_shared_tables=True)

	# Uploading table from fits file:
	#job = Gaia.upload_table(upload_resource=filename, table_name=name,
	# 			format="fits")

	# Uploading table from an astropy Table:
	a = [1,2,3]
	b = ['a','b','c']
	table = Table([a,b], names=['col1','col2'], meta={'meta':'first table'})
	try:
		job = Gaia.delete_user_table("tablename2")
		#Gaia.upload_table(upload_resource=gaia3, table_name=name)
		Gaia.upload_table(upload_resource=filename, table_name=name)
		full_qualified_table_name = f'user_rperei01.{name}'
		query = 'select * from ' + full_qualified_table_name
		job = Gaia.launch_job(query=query)
		results = job.get_results()
		Gaia.logout()
		return True

	except ValueError:
		return False

def gaia_codes(gaia3, gaia3_valid, lind20, outname):
	
	# 1) Parallax zero-point correction:
	# https://gitlab.com/icc-ub/public/gaiadr3_zeropoint
	# 2) G-band magnitude correction (sources with 6-parameters):
	# https://github.com/agabrown/gaiaedr3-6p-gband-correction
	# 3) Corrected flux excess factor:
	# https://github.com/agabrown/gaiaedr3-flux-excess-correction

	#gaia3_valid = gaia3[gaia3['astrometric_params_solved'] > 3]
	print('- All stars with 5p or 6p, have a valid parallax? :: ', end='')
	print(len(gaia3_valid) == len(gaia3[~np.isnan(gaia3['parallax'])]), '\n')

	#.Decide which result to use (extrap or non) -> EXTRAPOLATED!
	# from 6423 (out of 16268) with valid parallax:
	# -> 414 zpvals are extrapolated (303 < 1.24, 86 > 1.72)
	# -> 25-27 stars result a NaN extrap. correction (303+86+25 = 414 stars)
	
	# zpvals: correção de ponto zero subtraída da paralaxe de cada estrela!
	zpvals_extrap, zpvals_non = gaia_zpt_correction(gaia3)
	# corr_parallax = np.copy(gaia3['parallax'][~np.isnan(gaia3['parallax'])])
	# for (idx, zpval) in enumerate(zpvals_extrap):
	# 	if not np.isnan(zpval):
	# 		corr_parallax[idx] = corr_parallax[idx] - zpval
	# 	elif not np.isnan(corr_parallax[idx]):
	# 		corr_parallax[idx] = np.nan

	corr_parallax, zp = np.copy(gaia3['parallax']), 0
	for (idx, p) in enumerate(gaia3['parallax']):
		corr_parallax[idx] = np.nan if np.isnan(p) else (p - zpvals_extrap[zp])
		if not np.isnan(p): zp += 1
	corr_parallax_clean = corr_parallax[~np.isnan(gaia3['parallax'])] # 27 nan

	if lind20: figs_lindegren20()	# not important!

	#gaia_test_gband(gaia3)			# not important!
	gmag_corr, gflux_corr, deltag, fratio = gaia_gband_correction(gaia3) # OK!
	phot_bp_rp_excess_factor_corr = gaia_flux_excess_correction(gaia3)	# OK!
	pmra_corr, pmdec_corr = gaia_pm_correction(gaia3, gmag_corr)

	#.Inserting four new columns at the end of gaia3
	c_par = Column(corr_parallax, name='corr_parallax', unit='mas')
	c_gmag = Column(gmag_corr, name='corr_phot_g_mean_mag', unit='mag')
	c_gflux = Column(gflux_corr, name='corr_phot_g_mean_flux', unit='electron s-1')
	c_excess = Column(phot_bp_rp_excess_factor_corr, name='corr_phot_bp_rp_excess_factor')
	c_pmra = Column(pmra_corr, name='corr_pmra', unit='mas / yr')
	c_pmdec = Column(pmdec_corr, name='corr_pmdec', unit='mas / yr')
	gaia3.add_columns([c_par, c_gmag, c_gflux, c_excess, c_pmra, c_pmdec])
	gaia3_valid = gaia3[gaia3['astrometric_params_solved'] > 3]
	gmag_corr2v = gaia3_valid['corr_phot_g_mean_mag']

	#.Apply quality flags: Riello+2020, Antoja+20, Fabricius+20 (Lindegren...)
	c0, c1, m = 0.0059898, 8.817841e-12, 7.618399	# Riello
	sigma_Cstar = c0 + c1*(np.array(gmag_corr)**m)
	sigma_Cstar2 = c0 + c1*(np.array(gmag_corr2v)**m)

	flag_ruwe = gaia3['ruwe'] < 1.4						# Antoja
	flag_ruwe2 = gaia3['ruwe'] < 2.0					# Antoja - LMC
	flag_tmp = abs(gaia3['corr_phot_bp_rp_excess_factor']) < 5*sigma_Cstar
	flag_exc = ((gmag_corr > 4.0) & (flag_tmp)) | (gmag_corr <= 4.0)
	flag_gmag = (gmag_corr <= 19.0)						# Fabricius
	flag_rp = gaia3['phot_rp_mean_mag'] <= 20.0			# Riello
	flag = flag_ruwe & flag_exc & flag_gmag & flag_rp
	# flag = (flag_ruwe & flag_exc & flag_gmag & flag_rp) | \
	# 		(flag_ruwe & flag_bright & flag_rp)
	#..Riello+2020: flag_exc normally excludes variable and extended sources
	flag_noExcess = flag_ruwe & flag_gmag & flag_rp	# !!!
	# flag_noExcess = (flag_ruwe & flag_gmag & flag_rp) | \
	# 				(flag_ruwe & flag_bright & flag_rp)
	# flag_noExcess = (flag_gmag & flag_rp) | \
	# 				(flag_bright & flag_rp)			# temporary -->> decidir aqui!!! 1.4 ou 2.0 ou nada?

	#.Save new table with corrected columns and a quality flag -> fits file
	c_ruwe = Column(flag_ruwe.astype(int), name='flag_ruwe_1p4')
	c_ruwe2 = Column(flag_ruwe2.astype(int), name='flag_ruwe_2p0')
	c_exc = Column(flag_exc.astype(int), name='flag_excess_factor')
	c_gmag = Column(flag_gmag.astype(int), name='flag_gmag')
	c_rp = Column(flag_rp.astype(int), name='flag_rpmag')
	c_flag = Column(flag.astype(int), name='all_flags')
	c_flag2 = Column(flag_noExcess.astype(int), name='flag_noExcess')
	gaia3.add_columns([c_ruwe, c_ruwe2, c_exc, c_gmag, c_rp, c_flag, c_flag2])
	outname = outname.replace('gaia3', 'gaia3c')

	#.Changing units to avoid Warnings when writing fits file
	warnings.simplefilter('ignore', category=AstropyWarning)
	for col in gaia3.colnames:
		if gaia3[col].unit == 'dex':
			gaia3[col].unit = u.dex()
		elif gaia3[col].unit == 'log(cm.s**-2)':
			gaia3[col].unit = u.dex(u.cm / u.s**2)
		elif gaia3[col].unit == 'electron / s':
			gaia3[col].unit = 1/u.s
	colunits = [gaia3[c].unit if gaia3[c].unit else "" for c in gaia3.colnames]
	gaia3.write(outname, format='fits', overwrite=True)
	
	#.Reapplying flags only in stars with valid astrometry
	vflag_ruwe = gaia3_valid['ruwe'] < 1.4
	vflag_tmp = abs(gaia3_valid['corr_phot_bp_rp_excess_factor']) < 5*sigma_Cstar2
	vflag_exc = ((gmag_corr2v > 4.0) & (vflag_tmp)) | (gmag_corr2v <= 4.0)
	vflag_gmag = (gmag_corr2v <= 19.0)
	vflag_rp = gaia3_valid['phot_rp_mean_mag'] <= 20.0
	vflag = vflag_ruwe & vflag_exc & vflag_gmag & vflag_rp
	# vflag = (vflag_ruwe & vflag_exc & vflag_gmag & vflag_rp) | \
	# 		(vflag_ruwe & vflag_bright & vflag_rp)

	#.Print the quantity and percentage of excluded stars 
	valid = 100*(len(gaia3_valid)/len(gaia3))
	print(f'\n\033[1m* Valid astrometry (5p, 6p): {len(gaia3_valid)} out of',\
			f'{len(gaia3)} ({valid:.2f}%) stars\033[0m')
	perc1 = 100*(len(vflag_ruwe[~vflag_ruwe])/len(gaia3_valid))
	print(f'- RUWE < 1.4: excludes {len(vflag_ruwe[~vflag_ruwe])} ({perc1:.1f}%) stars')
	perc2 = 100*(len(vflag_exc[~vflag_exc])/len(gaia3_valid))
	print(f'- Excess factor: excludes {len(vflag_exc[~vflag_exc])} ({perc2:.1f}%) stars')
	perc3 = 100*(len(vflag_gmag[~vflag_gmag])/len(gaia3_valid))
	print(f'- Gmag <= 19.0: excludes {len(vflag_gmag[~vflag_gmag])} ({perc3:.1f}%) stars')
	perc4 = 100*(len(vflag_rp[~vflag_rp])/len(gaia3_valid))
	print(f'- G_RP <= 20: excludes {len(vflag_rp[~vflag_rp])} ({perc4:.1f}%) stars')
	percf = 100*(len(vflag[vflag])/len(gaia3_valid))
	print(f'\033[1m* Combining all the flags: {len(vflag[vflag])} ({percf:.2f}%)',\
			'stars remaining!\033[0m')

	#.Sample color-magnitude diagram: G vs. G-RP
	plt.figure(figsize=(5,7.5))
	plt.scatter(gaia3['g_rp'], gaia3['corr_phot_g_mean_mag'], color='#CECECE',
				marker='.', s=0.5)
	plt.scatter(gaia3_valid[vflag]['g_rp'], gaia3_valid[vflag]['corr_phot_g_mean_mag'],
				color='red', marker='.', s=1)
	plt.gca().invert_yaxis()
	plt.gca().update({'xlabel':r'(G$-$G$_{\rm RP}$)', 'ylabel':'G'})
	plt.tight_layout()
	plt.show()

	return gaia3

def gaia_zpt_correction(gaia3):
	# Mateus: ADQL query can also be done in https://gaia.aip.de/query/ (csv)

	path = os.path.dirname(os.path.realpath(__file__))
	sys.path.insert(0, os.path.join(path, 'gaiadr3_zeropoint/'))
	from zero_point import zpt as gaia_zpt
	#.Initialises the tables containing the coefficients of the interpolations
	# for the Z5 and Z6 functions.
	gaia_zpt.load_tables()

	#.Read the data
	#data = pd.read_csv('ZPrandomG18.csv')
	data = gaia3.to_pandas()
	gmag = data['phot_g_mean_mag'].values
	nueffused = data['nu_eff_used_in_astrometry'].values
	psc = data['pseudocolour'].values	# wave number (6th parameter in 6p)
	ecl_lat = data['ecl_lat'].values
	soltype = data['astrometric_params_solved'].values #{3, 31, 95}

	#.Blindly use the get_zpt function -> ValueError
	#zpvals = gaia_zpt.get_zpt(gmag, nueffused, psc, ecl_lat, soltype)

	#.Remove the 2-p solutions (invalid astrometric solutions)
	valid = soltype > 3
	zpvals_extrap = gaia_zpt.get_zpt(gmag[valid], nueffused[valid], psc[valid],
							ecl_lat[valid], soltype[valid])
	'''As explained in Lindegren et al. 2020, the interpolations are only
		calibrated within the following intervals:
	1. G magnitude: 6 < phot_g_mean_mag < 21
	2. Colour:
		A. 1.1 < nu_eff_used_in_astrometry < 1.9 (5-p sources)
		B. 1.24 < pseudocolour < 1.72 (6-p sources)
	Outside these ranges, the zero-point obtained is an extrapolation.
	'''

	# UserWarning: The pseudocolour of some of the 6p source(s) is outside
	# the expected range (1.24-1.72 mag). The maximum corrections are reached
	# already at 1.24 and 1.72

	#.Try turning off the warnings: NaN for the sources outside the limits,
	# instead of an *extrapolated* zero point (worse)
	# (automatically deals with a mix of 5-p and 6-p solutions)
	zpvals_non = gaia_zpt.get_zpt(gmag[valid], nueffused[valid], psc[valid], 
                   		ecl_lat[valid], soltype[valid], _warnings=False)
	print(f'\nzpvals_extrap ({len(zpvals_extrap)}) = {zpvals_extrap}')
	print(f'zpvals_non ({len(zpvals_non)}) = {zpvals_non}')

	#..separate 5p and 6p sources
	fivep = (soltype == 31)
	sixp = (soltype == 95)

	#..query the zero-point for individual stars
	zpvals5p = gaia_zpt.get_zpt(gmag[fivep][0], nueffused[fivep][0], psc[fivep][0],
							ecl_lat[fivep][0], soltype[fivep][0])
	zpvals6p = gaia_zpt.get_zpt(gmag[sixp][0], nueffused[sixp][0], psc[sixp][0],
							ecl_lat[sixp][0], soltype[sixp][0])
	#print(f'zpvals5p, zpvals6p = {zpvals5p}, {zpvals6p}')

	#.Using the Pandas wrapper ->> slower alternative to get_zpt() function
	data_valid = data.query('astrometric_params_solved>3').copy()
	data_valid['zpt'] = data_valid.apply(gaia_zpt.zpt_wrapper, axis=1)
	#print(data_valid.zpt)
	data_valid.hist(column='zpt', by='astrometric_params_solved', bins=50,
					figsize=(10,5), sharex=True)#, constrained_layout=True)
	plt.show()

	return (zpvals_extrap, zpvals_non)

def figs_lindegren20():

	from zero_point import zpt as gaia_zpt
	#.Reproducing Figures 21 and 22 of Lindegren et al. (2020)
	import scipy.stats
	# define array in the space of G vs colour (nu_eff or pseudocolour)
	g_vector = np.linspace(6, 21, 100)	# (100,) -> delta = 0.15
	c_vector = np.linspace(1.1, 1.9,  100) # (100,) -> delta = 0.008

	# *meshgrid* (sparse=False, indexing='xy'): Make N-D coordinate arrays for
	# vectorized evaluations of N-D scalar/vector fields over N-D grids, given
	# one-dimensional coordinate arrays x1, x2,…, xn. Basically, it creates a
	# rectangular grid out of an array of x values and an array of y values.
	G, C = np.meshgrid(g_vector,c_vector) # (100, 100); (100, 100)

	beta_90 = np.vstack((G.flat, C.flat, -90*np.ones_like(G.flat))).T # (10000, 3)
	beta0 = np.vstack((G.flat, C.flat, 0*np.ones_like(G.flat))).T	# (10000, 3)
	beta90 = np.vstack((G.flat,C.flat,90*np.ones_like(G.flat))).T	# (10000, 3)

	''' beta_90, beta0, beta90 -> -90, 0, 90
	array([[  6.        ,   1.1       , -90.        ],
       [  6.15151515,   1.1       , -90.        ],
       [  6.3030303 ,   1.1       , -90.        ],
       ...,
       [ 20.6969697 ,   1.9       , -90.        ],
       [ 20.84848485,   1.9       , -90.        ],
       [ 21.        ,   1.9       , -90.        ]])
	'''

	#..Figure 21 -> computing z5 for beta = -90, 0, 90 with get_zpt()
	z5_beta_90 = gaia_zpt.get_zpt(G.flatten(), C.flatten(), C.flatten(),
                         	-90*np.ones_like(G.flat), 31*np.ones_like(G.flat))#,
	#						_warnings=False)
	z5_beta0 = gaia_zpt.get_zpt(G.flatten(), C.flatten(), C.flatten(),
							0*np.ones_like(G.flat), 31*np.ones_like(G.flat))
	z5_beta90 = gaia_zpt.get_zpt(G.flatten(), C.flatten(), C.flatten(),
							90*np.ones_like(G.flat), 31*np.ones_like(G.flat))
	print(z5_beta_90, z5_beta0, z5_beta90)
	
	#...Beta = -90 deg (Fig. 21)
	plt.pcolormesh(C, G, 1000*z5_beta_90.reshape(G.shape), shading='nearest',
    			cmap='jet',vmin=-80, vmax=30)
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	cbar=plt.colorbar()
	cbar.set_label(r'Z5 [$\mu$as]')
	plt.xlabel(r'Effective wavenumber [$\mu$m$^{-1}$]')
	plt.ylabel(r'Magnitude [mag]')
	plt.annotate(r'a: $\beta$=-90$^\circ$',xy=(1.85,7),fontweight='bold')
	plt.tight_layout()
	plt.show()

	#...Beta = 0 deg (Fig. 21)
	plt.pcolormesh(C, G, 1000*z5_beta0.reshape(G.shape), shading='nearest',
    			cmap='jet',vmin=-80, vmax=30)
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	cbar=plt.colorbar()
	cbar.set_label(r'Z5 [$\mu$as]')
	plt.xlabel(r'Effective wavenumber [$\mu$m$^{-1}$]')
	plt.ylabel(r'Magnitude [mag]')
	plt.annotate(r'b: $\beta$=0$^\circ$',xy=(1.85,7),fontweight='bold')
	plt.tight_layout()
	plt.show()

	#...Beta = +90 deg (Fig. 21)
	plt.pcolormesh(C, G, 1000*z5_beta90.reshape(G.shape), shading='nearest',
				cmap='jet',vmin=-80, vmax=30)
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	cbar=plt.colorbar()
	cbar.set_label(r'Z5 [$\mu$as]')
	plt.xlabel(r'Effective wavenumber [$\mu$m$^{-1}$]')
	plt.ylabel(r'Magnitude [mag]')
	plt.annotate(r'b: $\beta$=+90$^\circ$',xy=(1.85,7),fontweight='bold')
	plt.tight_layout()
	plt.show()

	#..Figure 22 -> computing z6 for beta = -90, 0, 90 with get_zpt()
	z6_beta_90 = gaia_zpt.get_zpt(G.flatten(), C.flatten(), C.flatten(),
                         	-90*np.ones_like(G.flat), 95*np.ones_like(G.flat))#,
	#						_warnings=False)
	z6_beta0 = gaia_zpt.get_zpt(G.flatten(), C.flatten(), C.flatten(),
							0*np.ones_like(G.flat), 95*np.ones_like(G.flat))
	z6_beta90 = gaia_zpt.get_zpt(G.flatten(), C.flatten(), C.flatten(),
							90*np.ones_like(G.flat), 95*np.ones_like(G.flat))
	print(z6_beta_90, z6_beta0, z6_beta90)

	#...Beta = -90 deg (Fig. 22)
	plt.pcolormesh(C, G, 1000*z6_beta_90.reshape(G.shape), shading='nearest',
        		cmap='jet', vmin=-80, vmax=30)
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	cbar=plt.colorbar()
	cbar.set_label(r'Z6 [$\mu$as]')
	plt.xlabel(r'Effective wavenumber [$\mu$m$^{-1}$]')
	plt.ylabel(r'Magnitude [mag]')
	plt.annotate(r'a: $\beta$=-90$^\circ$',xy=(1.85,7),fontweight='bold')
	plt.tight_layout()
	plt.show()

	#...Beta = 0 deg (Fig. 22)
	plt.pcolormesh(C, G, 1000*z6_beta0.reshape(G.shape), shading='nearest',
        		cmap='jet', vmin=-80, vmax=30)
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	cbar=plt.colorbar()
	cbar.set_label(r'Z6 [$\mu$as]')
	plt.xlabel(r'Effective wavenumber [$\mu$m$^{-1}$]')
	plt.ylabel(r'Magnitude [mag]')
	plt.annotate(r'b: $\beta$=0$^\circ$',xy=(1.85,7),fontweight='bold')
	plt.tight_layout()
	plt.show()

	#...Beta = 90 deg (Fig. 22)
	plt.pcolormesh(C, G, 1000*z6_beta90.reshape(G.shape), shading='nearest',
           		cmap='jet', vmin=-80, vmax=30)
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	cbar=plt.colorbar()
	cbar.set_label(r'Z6 [$\mu$as]')
	plt.xlabel(r'Effective wavenumber [$\mu$m$^{-1}$]')
	plt.ylabel(r'Magnitude [mag]')
	plt.annotate(r'b: $\beta$=+90$^\circ$',xy=(1.85,7),fontweight='bold')
	plt.tight_layout()
	plt.show()

def correct_gband(bp_rp, astrometric_params_solved, phot_g_mean_mag, phot_g_mean_flux):
		"""
		Correct the G-band fluxes/magns for the input list of Gaia EDR3 data.

		Parameters
		----------

		bp_rp: float, numpy.ndarray
			The (BP-RP) colour listed in the Gaia EDR3 archive.
		astrometric_params_solved: int, numpy.ndarray {3, 31, 95}
			The astrometric solution type listed in the Gaia EDR3 archive.
		phot_g_mean_mag: float, numpy.ndarray
			The G-band magnitude as listed in the Gaia EDR3 archive.
		phot_g_mean_flux: float, numpy.ndarray
			The G-band flux as listed in the Gaia EDR3 archive.
			
		Returns
		-------

		The corrected G-band magnitudes and fluxes. The corrections are only
		applied to sources with a 6-parameter astrometric solution fainter than
		G=13, for which a (BP-RP) colour is available.
		"""

		if np.isscalar(bp_rp) or np.isscalar(astrometric_params_solved) or \
				np.isscalar(phot_g_mean_mag) or np.isscalar(phot_g_mean_flux):
			bp_rp = np.float64(bp_rp)
			astrometric_params_solved = np.int64(astrometric_params_solved)
			phot_g_mean_mag = np.float64(phot_g_mean_mag)
			phot_g_mean_flux = np.float64(phot_g_mean_flux)

		if not (bp_rp.shape == astrometric_params_solved.shape ==
				phot_g_mean_mag.shape == phot_g_mean_flux.shape):
			raise ValueError('Function parameters must be of the same shape!')

		do_not_correct = np.isnan(bp_rp) | (phot_g_mean_mag<=13) | (astrometric_params_solved != 95)
		bright_correct = np.logical_not(do_not_correct) & (phot_g_mean_mag>13) & (phot_g_mean_mag<=16)
		faint_correct = np.logical_not(do_not_correct) & (phot_g_mean_mag>16)
		bp_rp_c = np.clip(bp_rp, 0.25, 3.0)

		# factor = A - B*bprp_c + C*bprp_c^2 - D*bprp_c^3
		corr_factor = np.ones_like(phot_g_mean_mag)
		corr_factor[faint_correct] = 1.00525 - 0.02323*bp_rp_c[faint_correct] + \
			0.01740*np.power(bp_rp_c[faint_correct],2) - 0.00253*np.power(bp_rp_c[faint_correct],3)
		corr_factor[bright_correct] = 1.00876 - 0.02540*bp_rp_c[bright_correct] + \
			0.01747*np.power(bp_rp_c[bright_correct],2) - 0.00277*np.power(bp_rp_c[bright_correct],3)

		gmag_corrected = phot_g_mean_mag - 2.5*np.log10(corr_factor)
		gflux_corrected = phot_g_mean_flux * corr_factor

		return gmag_corrected, gflux_corrected

def gaia_test_gband(gaia3):
	# Source: github.com/agabrown/gaiaedr3-6p-gband-correction/blob/main/GCorrectionCode.ipynb
	#.Code throws the expected exception when input shapes do not match?
	print('\n- gaia_test_gband:', end=' ')
	from GaiaEDR3_data import gaia_gband_correction
	try:
		correct_gband(np.float(gaia3['bp_rp'][1]), 
					gaia3['astrometric_params_solved'][1:3],
					np.float(gaia3['phot_g_mean_mag'][1]), 
					np.float(gaia3['phot_g_mean_flux'][1]))
		print('Previous line should have failed!')
	except ValueError:
		print('1', end=' ') #('Expect to land here')
		
	try:
		correct_gband(gaia3['bp_rp'][1:3], 
					np.int(gaia3['astrometric_params_solved'][1]),
					np.float(gaia3['phot_g_mean_mag'][1]), 
					np.float(gaia3['phot_g_mean_flux'][1]))
		print('Previous line should have failed!')
	except ValueError:
		print('2', end=' ') #('Expect to land here')
		
	try:
		correct_gband(np.float(gaia3['bp_rp'][1]), 
					np.int(gaia3['astrometric_params_solved'][1]),
					gaia3['phot_g_mean_mag'][1:3], 
					np.float(gaia3['phot_g_mean_flux'][1]))
		print('Previous line should have failed!')
	except ValueError:
		print('3', end=' ') #('Expect to land here')
		
	try:
		correct_gband(np.float(gaia3['bp_rp'][1]), 
					np.int(gaia3['astrometric_params_solved'][1]),
					np.float(gaia3['phot_g_mean_mag'][1]), 
					gaia3['phot_g_mean_flux'][1:3])
		print('Previous line should have failed!')
	except ValueError:
		print('4 -> OK!', end='\n\n') #('Expect to land here')

def gaia_gband_correction(gaia3):
	"""
	G-band photometry corrections for Gaia EDR3 sources with 6-parameter
	astrometric solutions.

	The Jupyter notebook in this repository presents a Python function for
	calculating the corrections to the G-band photometry for Gaia EDR3 sources
	with a 6-parameter astrometric solution. The code is listed in the appendix
	of Gaia Collaboration et al. (2020), and implements the formulae presented
	in Riello et al. (2020).

	(c) Anthony G.A. Brown, Leiden University
	"""

	# Important columns: select source_id, astrometric_params_solved, bp_rp,
	# phot_g_mean_mag, phot_g_mean_flux

	#.Verify that the code works on the example data from the archive
	#edr3data = Table.read('GCorrection_inputs.fits', format='fits')
	gmag_corr, gflux_corr = correct_gband(gaia3['bp_rp'], gaia3['astrometric_params_solved'],
								gaia3['phot_g_mean_mag'], gaia3['phot_g_mean_flux'])
	
	G_null = gaia3['phot_g_mean_mag'][np.isnan(gaia3['phot_g_mean_mag'])].size
	mag_corr_null = gmag_corr[np.isnan(gmag_corr)].size
	flux_corr_null = gflux_corr[np.isnan(gflux_corr)].size
	print(f'Number of sources with G=null: {G_null}')
	print(f'Number of sources with G_corr=null: {mag_corr_null}')
	print(f'Number of sources with Gflux_corr=null: {flux_corr_null}\n')

	not_corrected = np.isnan(gaia3['bp_rp']) | (gaia3['phot_g_mean_mag']<=13) | \
					(gaia3['astrometric_params_solved'] != 95)
	faint_corrected = np.logical_not(not_corrected) & (gaia3['phot_g_mean_mag']>16)
	bright_corrected = np.logical_not(not_corrected) & (gaia3['phot_g_mean_mag']>13) & \
						(gaia3['phot_g_mean_mag']<=16)
	
	deltag = gmag_corr - gaia3['phot_g_mean_mag']
	fratio = gflux_corr / gaia3['phot_g_mean_flux']

	print('Range in delta_G where no correction is expected (should all be zeros):',\
		f'{np.nanmin(deltag[not_corrected])}--{np.nanmax(deltag[not_corrected])}')
	print('Range in flux ratio where no correction is expected (should all be ones):',\
		f'{np.nanmin(fratio[not_corrected])}--{np.nanmax(fratio[not_corrected])}')

	fig, ax = plt.subplots(1, 1, figsize=(8,4.7))
	ax.plot(gaia3['bp_rp'][not_corrected], deltag[not_corrected], '.', label=\
			'No correction expected')
	ax.plot(gaia3['bp_rp'][bright_corrected], deltag[bright_corrected], '.', \
			label='$13 < G \leq 16$')
	ax.plot(gaia3['bp_rp'][faint_corrected], deltag[faint_corrected], '.', \
			label='$G>16$')
	ax.set_xlabel('$(G_\mathrm{BP}-G_\mathrm{RP})$')
	ax.set_ylabel('$\Delta G$')
	ax.legend(loc='lower left')
	plt.tight_layout()
	plt.show()

	#.Verify that the code works for scalar inputs
	# comp_values_mags, comp_values_flux = gmag_corr[0:10], gflux_corr[0:10]
	# print(f"{'Vector G':10s} {'Scalar G':10s}  {'Vector flux':15s} {'Scalar flux':15s}")
	# for i in range(0,10):
	# 	corrg, corrf = correct_gband(np.float(gaia3['bp_rp'][i]), 
	# 							np.int(gaia3['astrometric_params_solved'][i]),
	# 							np.float(gaia3['phot_g_mean_mag'][i]), 
 	# 							np.float(gaia3['phot_g_mean_flux'][i]))
	# 	print(f'{corrg:10.6f} {comp_values_mags[i]:10.6f}  {corrf:15.6f} '+\
	# 			f' {comp_values_flux[i]:15.6f}')
	
	return (gmag_corr, gflux_corr, deltag, fratio)

def gaia_flux_excess_correction(gaia3):
	"""
	Calculation of the corrected flux excess factor for Gaia EDR3 sources.

	This notebook contains example code for the calculation of the Gaia EDR3
	CORRECTED flux excess factor according to Eq. 6 in Riello et al. (2020,
	view on: https://doi.org/10.1051/0004-6361/202039587).

	(c) Anthony G.A. Brown, Leiden University
	"""

	#.Define a function that calculates the corrected flux excess factor
	def correct_flux_excess_factor(bp_rp, phot_bp_rp_excess_factor):
		"""
		Calculate the corrected flux excess factor for the input Gaia EDR3 data.
		
		Parameters
		----------
		
		bp_rp: float, numpy.ndarray
			The (BP-RP) colour listed in the Gaia EDR3 archive.
		phot_bp_rp_excess_factor: float, numpy.ndarray
			The flux excess factor listed in the Gaia EDR3 archive.
			
		Returns
		-------
		
		Corr_value for the flux excess factor, which is zero for "normal" stars.
		"""
		
		if np.isscalar(bp_rp) or np.isscalar(phot_bp_rp_excess_factor):
			bp_rp = np.float64(bp_rp)
			phot_bp_rp_excess_factor = np.float64(phot_bp_rp_excess_factor)
		
		if bp_rp.shape != phot_bp_rp_excess_factor.shape:
			raise ValueError('Function parameters must be of the same shape!')
			
		do_not_correct = np.isnan(bp_rp)
		bluerange = np.logical_not(do_not_correct) & (bp_rp < 0.5)
		greenrange = np.logical_not(do_not_correct) & (bp_rp >= 0.5) & (bp_rp < 4.0)
		redrange = np.logical_not(do_not_correct) & (bp_rp > 4.0)
		
		correction = np.zeros_like(bp_rp)
		correction[bluerange] = 1.154360 + 0.033772*bp_rp[bluerange] + \
								0.032277*np.power(bp_rp[bluerange], 2)
		correction[greenrange] = 1.162004 + 0.011464*bp_rp[greenrange] + \
								0.049255*np.power(bp_rp[greenrange], 2) \
								- 0.005879*np.power(bp_rp[greenrange], 3)
		correction[redrange] = 1.057572 + 0.140537*bp_rp[redrange]
		
		return phot_bp_rp_excess_factor - correction

	#.Test on the example data downloaded from the archive
	phot_bp_rp_excess_factor_corr = correct_flux_excess_factor(gaia3['bp_rp'],
									gaia3['phot_bp_rp_excess_factor'])
	fig, (axa, axb) = plt.subplots(1, 2, figsize=(9,4.3))

	axa.hexbin(gaia3['bp_rp'], gaia3['phot_bp_rp_excess_factor'], bins='log',
			extent=[-4,8,-1.5,5], cmap='cividis', mincnt=1, gridsize=(240,130))
	axa.set_xlabel('$(G_\mathrm{BP}-G_\mathrm{RP}$)')
	axa.set_ylabel('$C$')
	axa.set_title('Raw flux excess factor')

	axb.hexbin(gaia3['bp_rp'], phot_bp_rp_excess_factor_corr, bins='log',
			extent=[-4,8,-1.5,5], cmap='cividis', mincnt=1, gridsize=(240,130))
	axb.set_xlabel('$(G_\mathrm{BP}-G_\mathrm{RP}$)')
	axb.set_ylabel('$C^*$')
	axb.set_title('Corrected flux excess factor')
	plt.tight_layout()
	plt.show()

	#.Verify that the code works on scalar inputs
	# comparison_values = phot_bp_rp_excess_factor_corr[0:10]
	# print(f"\n{'Vector':10s} {'Scalar':10s}")
	# for i in range(0,10):
	# 	corr = correct_flux_excess_factor(np.float(gaia3['bp_rp'][i]), 
	# 					np.float(gaia3['phot_bp_rp_excess_factor'][i]))
	# 	print(f'{corr:10.6f} {comparison_values[i]:10.6f}')
	
	#.Verify that the code throws an exception when input shapes do not match
	# try:
	# 	correct_flux_excess_factor(np.float(gaia3['bp_rp'][1]), gaia3['phot_bp_rp_excess_factor'][1:3])
	# 	print('Previous line should have failed!')
	# except ValueError:
	# 	print('\nExpect to land here')
		
	# try:
	# 	correct_flux_excess_factor(gaia3['bp_rp'][1:3], np.float(gaia3['phot_bp_rp_excess_factor'][1]))
	# 	print('Previous line should have failed!')
	# except ValueError:
	# 	print('Expect to land here\n')

	return phot_bp_rp_excess_factor_corr

def gaia_pm_correction(gaia3, gmag_corr):
	"""gaia_pm_correction(gaia3) -> Cantat-Gaudin & Brandt (2021, A&A)
	- Correcting the proper motion bias of the bright Gaia EDR3 sources.

	Input: gaia3 catalogue from Gaia EDR3, with format astropy.Table
	Output: corrected proper motions.
	"""

	pmra = np.array(gaia3['pmra'])	# 16268 (6423 not nan)
	pmdec = np.array(gaia3['pmdec'])
	ra, dec = np.array(gaia3['ra']), np.array(gaia3['dec'])
	#G = np.array(gaia3['corr_phot_g_mean_mag'])	# já corrigido!
	G = np.array(gmag_corr)
	
	def sind(x):
		return np.sin(np.radians(x))
	def cosd(x):
		return np.cos(np.radians(x))

	table1 = """ 0.0	9.0		18.4	33.8	-11.3
				 9.0	9.5		14.0	30.7	-19.4
				 9.5	10.0	12.8	31.4	-11.8
				 10.0	10.5	13.6	35.7	-10.5
				 10.5	11.0	16.2	50.0	2.1
				 11.0	11.5	19.4	59.9	0.2
				 11.5	11.75	21.8	64.2	1.0
				 11.75	12.0	17.7	65.6	-1.9
				 12.0	12.25	21.3	74.8	2.1
				 12.25	12.5	25.7	73.6	1.0 
				 12.5	12.75	27.3	76.6	0.5
				 12.75	13.0	34.9	68.9	-2.9 """

	table1 = np.fromstring(table1, sep='\t').reshape((12,5)).T
	Gmin, Gmax = table1[0], table1[1]

	pmraCorr, pmdecCorr = np.ones_like(pmra), np.ones_like(pmdec)	# /1000
	for i in range(len(pmra)):
		if G[i] >= 13 or np.isnan(G[i]):
			pmraCorr[i], pmdecCorr[i] = pmra[i], pmdec[i]
		else:
			#.pick the appropriate omegaXYZ for the source’s magnitude:
			omegaX = table1[2][(Gmin <= G[i]) & (Gmax > G[i])][0]
			omegaY = table1[3][(Gmin <= G[i]) & (Gmax > G[i])][0]
			omegaZ = table1[4][(Gmin <= G[i]) & (Gmax > G[i])][0]
			pmraCorr[i] = -1*sind(dec[i])*cosd(ra[i])*omegaX - \
						sind(dec[i])*sind(ra[i])*omegaY + cosd(dec[i])*omegaZ
			pmdecCorr[i] = sind(ra[i])*omegaX -cosd(ra[i])*omegaY

	return pmra-pmraCorr/1000., pmdec-pmdecCorr/1000.

### Xmatch:::
# full_qualified_table_name = 'user_<your_login_name>.my_sources'
# xmatch_table_name = 'xmatch_table'
#Gaia.cross_match(full_qualified_table_name_a=full_qualified_table_name,
#					full_qualified_table_name_b='gaiadr2.gaia_source',
#					results_table_name=xmatch_table_name, radius=1.0)
# VER CONTINUAÇÃO NO DOC

if __name__ == '__main__':
	main()
