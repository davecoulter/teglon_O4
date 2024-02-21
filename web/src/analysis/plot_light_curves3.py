import json
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np

gw170817_json_file = open("./web/models/kne/GW170817.json")
gw170817_json_data = json.load(gw170817_json_file)
merger_time = float(gw170817_json_data["GW170817"]["timeofmerger"][0]["value"])
sss17a_json_phot = gw170817_json_data["GW170817"]["photometry"]

model_table = Table.read('./web/models/kne/villar_2comp_0.0230_0.2560_0.5000_0.0230_0.2560_0.5000.dat', format='ascii.ecsv')

gw170817_r_band_phot = []
gw170817_r_band_phot_err = []
gw170817_r_band_phot_delta_mjd = []
all_data = {}
for p in sss17a_json_phot:
    if ('band' in p) and ('u_time' in p) and ('upperlimit' not in p) and ('system' in p) and ('e_magnitude' in p) and (p['u_time'] == "MJD") and ('model' not in p) and (p['system'] == 'AB'):
        filt = p['band']
        # if filt in ['K','H','J','y','z','i','r','V','g','B']:
        if filt in ['r']:
            if filt not in all_data.keys():
                all_data[filt]={'time':[], 'mag':[], 'magerr':[]}

            all_data[filt]['time'].append(float(p['time']) - merger_time)
            all_data[filt]['mag'].append(float(p['magnitude']))
            all_data[filt]['magerr'].append(float(p['e_magnitude']))

for key in all_data.keys():
    sorted_indices = np.argsort(all_data[key]['time'])
    all_data[key]['time'] = np.asarray(all_data[key]['time'])[sorted_indices]
    all_data[key]['mag'] = np.asarray(all_data[key]['mag'])[sorted_indices]
    all_data[key]['magerr'] = np.asarray(all_data[key]['magerr'])[sorted_indices]

# convert to app map
mu = 5*np.log10(43.2 * 1e6)-5

fig = plt.figure(figsize=(10, 10), dpi=600)
ax = fig.add_subplot(111)

filt_map = {
    'K':'ukirt_K',
            'H':'ukirt_H',
            'J':'ukirt_J',
            # 'y':'PS1_Y',
            'z':'PS1_z',
            'i':'PS1_i',
            'r':'PS1_r',
            'V':'johnson_V',
            'g':'PS1_g',
            'B':'johnson_B'}

# # Name, A_lambda
# SDSS u, 0.5126587253999999
# SDSS g, 0.3994601958
# SDSS r, 0.276344701
# SDSS i, 0.2053537428
# SDSS z, 0.1527454518
# Landolt B, 0.43852336359999994
# Landolt V, 0.3316136412
# Landolt R, 0.2623158234
# Landolt I, 0.18201259299999997
# UKIRT J, 0.0857454674
# UKIRT H, 0.054301431399999996
# UKIRT K, 0.0365234572
# Clear, 0.110054126
# ATLAS cyan, 0.35151045967199995
# ATLAS orange, 0.25633734464759994
# PS1 w, 0.28311726260000003

extinctions = {'K':0.0365234572,
                'H':0.054301431399999996,
                'J':0.0857454674,
                # 'y':'PS1_Y',
                'z':0.1527454518,
                'i':0.2053537428,
                'r':0.276344701,
                'V':0.3316136412,
                'g':0.3994601958,
                'B':0.43852336359999994}

colors = ['darkred','red','orangered','orange','goldenrod','gold','green','dodgerblue',
    'blue','magenta']

offsets=[-2.0,-1.5,-1.0,-0.5,-0.25,0.0,0.25,0.5,1.0,1.5,2.0]

# for i,filt in enumerate(['K','H','J','y','z','i','r','V','g','B']):
for i,filt in enumerate(['r']): # no y-band

    if filt not in all_data.keys():
        print(f'No data {filt}')

    model_filt = filt_map[filt]

    time = all_data[filt]['time']
    mag = all_data[filt]['mag']
    magerr = all_data[filt]['magerr']
    extinction = extinctions[filt]

    if offsets[i]==0.0:
        label=filt
    elif offsets[i] < 0:
        label=filt+str(offsets[i])
    else:
        label=filt+'+'+str(offsets[i])

    ax.errorbar(time, mag - mu + offsets[i] - extinction, yerr=magerr, color=colors[i], fmt='*',
        linestyle='None', label=label)
    ax.plot(model_table['time'], model_table[model_filt] + offsets[i], color=colors[i],
        linestyle='solid')

plt.ylim([-8,-18])
plt.xlim([0,22])
plt.legend()


output_file = "villar_models.png"
output_path = "./web/src/analysis/%s" % output_file
fig.savefig(output_path, bbox_inches='tight', format="png")
plt.close('all')
print("... Done.")