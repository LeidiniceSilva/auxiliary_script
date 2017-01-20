# -*- coding: utf-8 -*-

import os
import numpy as np
import ftplib
import ftputil
from netCDF4 import Dataset
from PyFuncemeClimateTools import DefineGrid as Dg
from datetime import datetime
from dateutil.relativedelta import relativedelta


def write_3d_nc(ncname, var_array, time_array, lat_array, lon_array, var_units,
                var_shortname, var_longname, time_units, missing_value=-999.):
    """Write 3 dimensional (time, lat, lon) netCDF4

    :param ncname: Name of the output netCDF
    :type ncname: string object
    :param var_array: 3d prec array (time, lat, lon)
    :type var_array: numpy array
    :param time_array: 1d array with unique value which is the prediction hour
    :type time_array: numpy array
    :param lat_array: 1d array which contains the lat values
    :type lat_array: numpy array
    :param lon_array: 1d array which contains the lon values
    :type lon_array: numpy array
    :param var_units: variable units
    :type var_units: string object
    :param var_shortname: variable short name
    :type var_shortname: string object
    :param var_longname: variable long name
    :type var_longname: string object
    :param time_units: time units
    :type time_units: string object
    :param missing_value: missing value
    :type missing_value: float
    """

    foo = Dataset(ncname, 'w', format='NETCDF4_CLASSIC')

    foo.createDimension('time', time_array.shape[0])
    foo.createDimension('lat', lat_array.shape[0])
    foo.createDimension('lon', lon_array.shape[0])

    times = foo.createVariable('time', 'f8', 'time')
    times.units = time_units
    times.calendar = 'standard'
    times[:] = time_array

    laty = foo.createVariable('lat', 'f4', 'lat')
    laty.units = 'degrees_north'
    laty.long_name = 'latitude'
    laty[:] = lat_array

    lonx = foo.createVariable('lon', 'f4', 'lon')
    lonx.units = 'degrees_east'
    lonx.long_name = 'longitude'
    lonx[:] = lon_array

    var = foo.createVariable(var_shortname, 'f', ('time', 'lat', 'lon'))
    var.units = var_units
    var.long_name = var_longname
    var.missing_value = missing_value
    var[:] = var_array

    foo.close()


def cut2region(indata, inlat, inlon, minlat, minlon, maxlat, maxlon):
    """Cut input data to the region of interest

    :param indata: 3d data array that will be cut
    :type indata: numpy array
    :param inlat: 1d lat array that will be latitude
    :type inlat: numpy array
    :param inlon: 1d lon array that will be longitude
    :type inlon: numpy array
    :param minlat: value of the minimun latitude
    :type minlat: float
    :param minlon: value of the minimun longitude
    :type minlon: float
    :param maxlat: value of the maximum latitude
    :type maxlat: float
    :param maxlon: value of the maxium longitude
    :type maxlon: float
    :return: Return data array cut, lon array cut and lat array cut to
    the region of interest
    """

    min_lat, min_lon, min_lat_index, min_lon_index = Dg.gridpoint(inlat, inlon,
                                                                  minlat,
                                                                  minlon)
    max_lat, max_lon, max_lat_index, max_lon_index = Dg.gridpoint(inlat, inlon,
                                                                  maxlat,
                                                                  maxlon)
    cut_data = indata[..., min_lat_index:max_lat_index,
                      min_lon_index:max_lon_index]
    cut_lat = inlat[min_lat_index:max_lat_index]
    cut_lon = inlon[min_lon_index:max_lon_index]

    return cut_data, cut_lat, cut_lon


def basin_dict(basin_name='', macro=False, micro=False, all_basins=False):
    """
    Access a basin/sub-basin dictionary which contains:
    -basin/sub-basin path directory;
    -sub-basin's basin name;
    -
    :param basin_name: Name of the drainage basin
    :type basin_name: str
    :param macro: If true, returns keys of macro-basins dictionary
    :type macro: bool
    :param micro: If true, returns keys of micro-basins dictionary of input
    basin_name
    :type micro: bool
    :param all_basins: If true, returns keys for both macro-basins and
    micro-basins
    :type all_basins: bool
    :return: values of specific input basin or keys of dictionary
    """

    hidropy_path = get_hidropy_dirpath()

    shape_path = hidropy_path + 'hidropy/shapes/basins/'

    basins_list = os.listdir(shape_path)

    list_basins = []

    for i in basins_list:
        sub_basin = sorted(os.listdir(shape_path + i))
        for j in sub_basin:
            if i == j.split('.')[0]:
                list_basins.append("{0} {1}{0}/{0}.asc {0} {0} 0\n"
                                   .format(i, shape_path))
            else:
                mi = '_'.join(j.split('.')[0].split('_')[len(i.split('_')):])
                list_basins.append("{2} {1}{0}/{0}_{2}.asc {0} {0}_{2} 1\n"
                                   .format(i, shape_path, mi))

    basins = {}
    for line in list_basins:
        line = line.split()
        if not line:
            continue
        basins[line[0]] = line[1:]

    if macro:
        keys = [i for i in basins.keys() if basins[i][3] == "0"]
        return sorted(keys)
    elif micro:
        keys = [i for i in basins.keys()
                if basins[i][3] == "1" and basins[i][1] == basin_name or
                basins[i][3] == "1" and basin_name is '']
        return sorted(keys)
    elif all_basins:
        keys = basins.keys()
        return sorted(keys)
    else:
        basin_values = basins[basin_name.lower()]
        return basin_values


def date2index(dates, nctime):
    """
    Compute the index in time array for a given date.

    :param dates: input datetime
    :type dates: datetime
    :param nctime: input netCDF time variable
    :type nctime: netCDF variable
    :return: index for input date
    """

    if len(nctime.units.split(' ')) == 4:
        time_units = datetime.strptime(nctime.units.split(' ')[2] +
                                       nctime.units.split(' ')[3],
                                       '%Y-%m-%d%H:%M:%S')
        if time_units.hour == 12:
            time_units = time_units - relativedelta(hours=12)
    else:
        time_units = datetime.strptime(nctime.units.split(' ')[2],
                                       '%Y-%m-%d')

    if nctime.units.split(' ')[0] == 'hours':
        delta_time = relativedelta(dates, time_units).hours
        date_index = np.where(nctime[:] == delta_time)[0][0]
    elif nctime.units.split(' ')[0] == 'days':
        delta_time = (dates - time_units).days
        date_index = np.where(nctime[:] == delta_time)[0][0]
    elif nctime.units.split(' ')[0] == 'months':
        delta_time = 12*(dates.year - time_units.year) + \
                     (dates.month - time_units.month)
        date_index = np.where(nctime[:] == delta_time)[0][0]
    elif nctime.units.split(' ')[0] == 'years':
        delta_time = relativedelta(dates, time_units).years
        date_index = np.where(nctime[:] == delta_time)[0][0]
    else:
        date_index = 0
        exit()

    return date_index


def num2date(time_values, units):

    """
    Compute dates of a given time array.

    :param time_values: time array
    :type time_values: numpy array
    :param units: netcdf time units
    :type units: str
    :return: array of datetime
    """

    if len(units.split(' ')) == 4:
        time_units = datetime.strptime(units.split(' ')[2] +
                                       units.split(' ')[3],
                                       '%Y-%m-%d%H:%M:%S')
    else:
        time_units = datetime.strptime(units.split(' ')[2],
                                       '%Y-%m-%d')
    init_index = time_values[0]

    dates_array = []

    for j, i in enumerate(time_values):
        index_diff = int(i - init_index)
        if units.split(' ')[0] == 'hours':
            delta_time = relativedelta(time_units, hours=index_diff)
            cur_date = time_units + relativedelta(years=delta_time.years,
                                                  months=delta_time.months,
                                                  days=delta_time.days,
                                                  hours=delta_time.hours)
        elif units.split(' ')[0] == 'days':
            delta_time = relativedelta(time_units, days=index_diff)
            cur_date = time_units + relativedelta(years=delta_time.years,
                                                  months=delta_time.months,
                                                  days=delta_time.days)
        else:
            delta_time = relativedelta(time_units, months=index_diff)
            cur_date = time_units + relativedelta(years=delta_time.years,
                                                  months=delta_time.months)

        dates_array.append(cur_date)

    dates_array = np.asarray(dates_array)

    return dates_array


def get_hidropy_dirpath():

    """
    :return: Local hidropy's directory path.
    """

    hidropy_dir = '/'.join(os.path.dirname(os.path.realpath(__file__))
                           .split('/')[0:-2]) + '/'

    return hidropy_dir


def upload_obs_opendap(basin_name, var_name, time_freq, database, local_dir):

    """
    Upload all accumulated thiessen precipitation from local database directory
    to opendap's databese directory
    :param basin_name: Name of the basin where Thiessen will be computed.
    :type basin_name: str
    :param var_name: Name of the variable.
    :type var_name: str
    :param time_freq: Time frequency, must be 'daily' or 'monthly'
    :type time_freq: str
    :param database: Database used. Must be "inmet", "inmet_ana", "chirps".
    :type database: str
    :param local_dir: Local directory where data to be uploaded are".
    :type local_dir: str
    """

    basin = basin_dict(basin_name)[1]
    ftp = ftplib.FTP('opendap4.funceme.br')
    ftp.login('ftpuser', 'Funceme')
    for i in ['calibration', 'operation']:
        try:
            ftp.cwd('/io/{2}/{1}/{0}/'.format(time_freq, i, database))
            local_path = '{0}/{3}/{2}/{1}/'\
                .format(local_dir, time_freq, i, database)
        except Exception, e:
            print e, '--> No local directory named ' \
                     '{0}/{3}/{2}/{1}/, skipping ...'\
                .format(local_dir, time_freq, i, database)
            break
        try:
            ftp.cwd('/io/{3}/{2}/{0}/{1}_thiessen'.format(time_freq, var_name,
                                                          i, database))
        except Exception, e:
            print e, '--> No directory named: /io/{3}/{2}/{0}/{1}_thiessen,' \
                     ' creating it'.format(time_freq, var_name, i, database)
            ftp.mkd('/io/{3}/{2}/{0}/{1}_thiessen'.format(time_freq, var_name,
                                                          i, database))
            ftp.cwd('/io/{3}/{2}/{0}/{1}_thiessen'.format(time_freq, var_name,
                                                          i, database))
        m = var_name + '_thiessen/{0}'.format(basin)
        try:
            ftp.cwd('/io/{3}/{2}/{0}/{1}'.format(time_freq, m, i, database))
        except Exception, e:
            print e, '--> No directory named: /io/{3}/{2}/{0}/{1}, ' \
                     'creating it'.format(time_freq, m, i, database)
            ftp.mkd('/io/{3}/{2}/{0}/{1}'.format(time_freq, m, i, database))
            ftp.cwd('/io/{3}/{2}/{0}/{1}'.format(time_freq, m, i, database))

        opendap_pr_files = ftp.nlst()
        local_pr_files = sorted(os.listdir(local_path + '{0}/'.format(m)))
        upload_pr_files = [local_path + '{0}/'.format(m) + fi for fi in
                           list(set(local_pr_files) - set(opendap_pr_files))]
        for n in upload_pr_files:
            if n.endswith('.nc'):
                with open(n, 'rb') as f:
                    ftp.storbinary('STOR {0}'.format(os.path.basename(n)), f)
        print 'Files have been uploaded to opendap4'
        ftp.quit()


def create_path(path):

    """
    :param path: path to save data
    :type path: str
    """

    if not os.path.exists(path):
        os.system('mkdir -p ' + path)


def obs2opendap(upload_files, basin, data_usage, var_name, time_freq,
                database):

    """
    Upload thiessen .asc files from database's local directory to
    databases's opendap directory

    :param upload_files: name of the files to be uploaded
    :type upload_files: str
    :param basin: name of the basin or sub-basin
    :type basin: str
    :param data_usage: Data usage, must be calibration or operation
    :type data_usage: str
    :param var_name: name of the basin or sub-basin
    :type var_name: str
    :param time_freq: Time frequency, must be 'daily' or 'monthly'
    :type time_freq: str
    :param database: Database used. Must be "inmet", "inmet_ana", "chirps".
    :type database: str
    """

    url_musf = "opendap4.funceme.br"
    ftp_user = "ftpuser"
    ftp_pass = "Funceme"

    file_name = os.path.basename(upload_files)

    dir_target = '/io/{0}/{1}/{2}/{3}_thiessen/{4}/'\
        .format(database, data_usage, time_freq, var_name, basin_dict(basin)[1])

    ftp = ftputil.FTPHost(url_musf, ftp_user, ftp_pass)

    if ftp.path.isdir(dir_target):
        ftp.upload(upload_files, dir_target + file_name)
        print ""
        print ("Uploading {0} ...".format(file_name))
    else:
        ftp.makedirs(dir_target)
        ftp.upload(upload_files, dir_target + file_name)
        print ""
        print ("Uploading {0} ...".format(file_name))

    ftp.close()


def model2opendap(upload_files, basin, var_name, database, init_date):
    """
    Upload thiessen .asc files from database's local directory to
    databases's opendap directory

    :param upload_files: name of the files to be uploaded
    :type upload_files: str
    :param basin: name of the basin or sub-basin
    :type basin: str
    :param var_name: name of the basin or sub-basin
    :type var_name: str
    :param database: Database used. Must be "gfs05", "eta15", "eta40".
    :type database: str
    :param init_date: Start date of forecast.
    :type init_date: datetime
    """

    url_musf = "opendap4.funceme.br"
    ftp_user = "ftpuser"
    ftp_pass = "Funceme"

    file_name = os.path.basename(upload_files)

    dir_target = '/io/{0}/{1}/{1}{2:02d}/{3}/'\
        .format(database, init_date.year, init_date.month, var_name)
    if var_name == 'pr_thiessen' or var_name == 'pet_thiessen':
        dir_target = '/io/{0}/{1}/{1}{2:02d}/{3}/{4}/'\
            .format(database, init_date.year, init_date.month, var_name,
                    basin_dict(basin)[1])
    ftp = ftputil.FTPHost(url_musf, ftp_user, ftp_pass)

    if ftp.path.isdir(dir_target):
        ftp.upload(upload_files, dir_target + file_name)
        print ("-> Uploading {0}".format(file_name))
    else:
        ftp.makedirs(dir_target)
        ftp.upload(upload_files, dir_target + file_name)
        print ("-> Uploading {0}".format(file_name))

    ftp.close()


def obs2owncloud(upload_files, basin, data_usage, var_name, time_freq, database,
                 owncloud_path):

    """
    Upload thiessen .asc files from database's local directory to
    databases's owncloud directory

    :param upload_files: name of the files to be uploaded
    :type upload_files: str
    :param basin: name of the basin or sub-basin
    :type basin: str
    :param data_usage: Data usage, must be calibration or operation
    :type data_usage: str
    :param var_name: name of the basin or sub-basin
    :type var_name: str
    :param time_freq: Time frequency, must be 'daily' or 'monthly'
    :type time_freq: str
    :param database: Database used. Must be "inmet", "inmet_ana", "chirps".
    :type database: str
    :param owncloud_path: path to local ownCloud directory
    :type owncloud_path: str
    """

    cloud_dir = owncloud_path + 'FUNCEME-Petrobras/dados/data-thiessen/{0}/' \
                                '{1}/{2}/{3}_thiessen/{4}/'\
        .format(database, data_usage, time_freq, var_name, basin_dict(basin)[1])

    create_path(cloud_dir)

    os.system('mv {0} {1}'.format(upload_files, cloud_dir))


def thiessen_moving_sum(array, window):

    ms_array = np.zeros((array.shape[0] - window + 1))
    for i in range(array.shape[0] - window + 1):
        ms_array[i] = np.ma.sum(array[i:i + window])
    return ms_array


def thiessen_moving_avg(array, window):

    ms_array = np.zeros((array.shape[0] - window + 1))
    for i in range(array.shape[0] - window + 1):
        ms_array[i] = np.ma.mean(array[i:i + window])
    return ms_array
