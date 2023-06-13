import os
import numpy as np
from astropy.table import Table


# abs_mag_array = np.linspace(-11.0, -23.0, 100)
abs_mag_array = np.linspace(-12.0, -21.0, 50)
dm_list_of_list = []
for am in abs_mag_array:
    # dm = []
    dm = list(np.logspace(np.log10(0.001), np.log10(2.0), 50))
    # if am > -19.5:
    #     # dm1 = [-2.0]
    #     # dm2 = list(np.linspace(-1.5, 1.5, 31))
    #     dm1 = [-2.0, -1.5]
    #     dm2 = list(np.linspace(-1.0, 1.0, 500))
    #     dm3 = [4.0, 5.0, 6.0, 7.0, 8.0, 10.0]
    #     dm += dm1
    #     dm += dm2
    #     dm += dm3
    # else:
    #     dm1 = [-2.0, -1.5]
    #     dm2 = list(np.linspace(-1.0, 10, 38))
    #     dm += dm1
    #     dm += dm2

    # print(len(dm))
    dm_list_of_list.append(dm)



# dm_array = np.linspace(-5.0, 10.0, 30)
time_array = np.linspace(0.00, 15.0, 200) # in days

# abs_mag_array = np.linspace(-16.0, -20.0, 40)
# dm_array = np.linspace(-1.0, 8.0, 40)
# time_array = np.linspace(0.00, 26.0, 251) # in days

abs_lcs_array = {}
for dm_index, abs_mag in enumerate(abs_mag_array):
    for dm in dm_list_of_list[dm_index]:
        key = (abs_mag, dm)
        lc = abs_mag + time_array * dm
        abs_lcs_array[key] = lc

print("Number of generated models: %s" % len(abs_lcs_array))


cols = ['time', 'sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z', 'Clear', 'johnson_B', 'johnson_V', 'johnson_R',
        'johnson_I', 'ukirt_J', 'ukirt_H', 'ukirt_K']
dtype = ['f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8']

j = 0
for i, (key, lc) in enumerate(abs_lcs_array.items()):

    current_sub_dir = "./web/models/linear/%s" % j
    if i % 25 == 0:
        j += 1
        current_sub_dir = "./web/models/linear/%s" % j
        os.mkdir(current_sub_dir)



    abs_mag = "%0.3f" % key[0]
    dm = "%0.3f" % key[1]

    meta = ["{key}={value}".format(key="M", value=key[0]), "{key}={value}".format(key="dM", value=key[1])]

    result_table = Table(dtype=dtype, names=cols)
    result_table.meta['comment'] = meta

    print("Writing linear model: (%s, %s)" % key)
    for i, epoch in enumerate(lc):
        result_table.add_row([time_array[i], epoch, epoch, epoch, epoch, epoch, epoch, epoch, epoch, epoch, epoch,
                              epoch, epoch, epoch])

    result_table.write("%s/%s_%s.dat" % (current_sub_dir, abs_mag.replace(".", "_"),
                                                      dm.replace(".", "_")), overwrite=True, format='ascii.ecsv')



