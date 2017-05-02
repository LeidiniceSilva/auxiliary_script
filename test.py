# -*- coding: utf-8 -*-

""" Verification of the precipitation ability of flow. """

import matplotlib as mpl; mpl.use('Agg')
import calendar
from datetime import date
from pylab import *
from hidropy.utils.hidropy_utils import create_path


__author__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = "07/04/2017"
__description__ = "Tests"

for k, yea in enumerate(range(2011, 2017)):
    aux = echam_corri[k::6]

    dat1 = date(yea, 8, 1)
    last_day_mon = calendar.monthrange(yea, 8)[1]
    dat2 = date(yea, 8, last_day_mon)
    new_start = dat1 + relativedelta(months=1)
    new_endd = dat2 + relativedelta(months=3)

    start_y = str(dat1)[0:4] + str(dat1)[5:7] + str(dat1)[8:10]
    new_y = str(new_start)[0:4] + str(new_start)[5:7] + str(new_start)[8:10]
    end_y = str(new_endd)[0:4] + str(new_endd)[5:7] + str(new_endd)[8:10]
    # print start_y, new_y, end_y
    # exit()

    if not os.path.exists(path_out2):
        create_path(path_out2)

    name_nc = write_thiessen(aux, new_y, end_y, 'monthly', 'pr', 'echam46_hind8110', 'fcst',
                             'corrigido_{0}'.format(basin_fullname), init_date=start_y, output_path=path_out2)


start_date = date(year1, 12, 1)
new_year = start_date + relativedelta(months=1)
new_year_y = str(new_year)[0:4] + str(new_year)[5:7] + str(new_year)[8:10]


last_day = calendar.monthrange(year1, 3)[1]
start_date = date(year1, 12, 1)
new_year = start_date + relativedelta(months=1)
new_start_date = date(year1, 12, last_day)
end_year = new_start_date + relativedelta(months=3)

new_year_y = str(new_year)[0:4] + str(new_year)[5:7] + str(new_year)[8:10]
end_year_y = str(end_year)[0:4] + str(end_year)[5:7] + str(end_year)[8:10]




