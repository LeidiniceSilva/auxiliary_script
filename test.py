# Write output thiessen in netCDF4 file
# for k, yea in enumerate(range(2011, 2017)):
#     aux = echam_corri[k::6]
#
#     dat1 = date(yea, 8, 1)
#     last_day_mon = calendar.monthrange(yea, 8)[1]
#     dat2 = date(yea, 8, last_day_mon)
#     new_start = dat1 + relativedelta(months=1)
#     new_endd = dat2 + relativedelta(months=3)
#
#     start_y = str(dat1)[0:4] + str(dat1)[5:7] + str(dat1)[8:10]
#     new_y = str(new_start)[0:4] + str(new_start)[5:7] + str(new_start)[8:10]
#     end_y = str(new_endd)[0:4] + str(new_endd)[5:7] + str(new_endd)[8:10]
#     # print start_y, new_y, end_y
#     # exit()
#
#     if not os.path.exists(path_out2):
#         create_path(path_out2)
#
#     name_nc = write_thiessen(aux, new_y, end_y, 'monthly', 'pr', 'echam46_hind8110', 'fcst',
#                              'corrigido_{0}'.format(basin_fullname), init_date=start_y, output_path=path_out2)















