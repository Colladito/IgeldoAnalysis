import cdsapi
import os    

c = cdsapi.Client()

cdsGribcode = dict()
## cdsGribcode['SFC_CAPE.128'] = "convective_available_potential_energy"
## cdsGribcode['SFC_151.128'] = "mean_sea_level_pressure"
## cdsGribcode['SFC_CIN.128'] = "convective_inhibition"
## cdsGribcode['SFC_K.128'] = "k_index"
cdsGribcode['SFC_TT.128'] = "total_totals_index"
## https://confluence.ecmwf.int/display/FUG/Section+9.6.2+Instability+Indices
months = ["01","02","03","04","05","06","07","08","09","10","11","12"]
for i in cdsGribcode : 
	for cdsInitYear in range(2019,2024):
		## ## Month
		## for cdsInitMonth in months:
		nameFile = "./ECMWF_ERA5_" + str(cdsInitYear) + "_" + i + ".nc"
		c.retrieve(
			'reanalysis-era5-single-levels',
			{
				'product_type': 'reanalysis',
				'variable':cdsGribcode[i],
				'year':str(cdsInitYear),
				## 'month':cdsInitMonth,
				'month': [
							'01', '02', '03',
							'04', '05', '06',
							'07', '08', '09',
							'10', '11', '12',
				],					
				'day': [
					'01', '02', '03',
					'04', '05', '06',
					'07', '08', '09',
					'10', '11', '12',
					'13', '14', '15',
					'16', '17', '18',
					'19', '20', '21',
					'22', '23', '24',
					'25', '26', '27',
					'28', '29', '30',
					'31',
				],
				'time': [
					'00:00', '01:00', '02:00',
					'03:00', '04:00', '05:00',
					'06:00', '07:00', '08:00',
					'09:00', '10:00', '11:00',
					'12:00', '13:00', '14:00',
					'15:00', '16:00', '17:00',
					'18:00', '19:00', '20:00',
					'21:00', '22:00', '23:00',
				],
				'area': [
					50, -10, 41.5,
					0,
				],
				'format': 'netcdf',
			},
			nameFile)
