# -*- coding: utf-8 -*-

__author__ = "Marcelo Rodrigues"
__email__ = "marcelorodriguesss@gmail.com"
__credits__ = ["Fco Vasconcelos", "Arthur Costa", "Aurélio Noronha"]

# TODO:
# Incluir função de probabilidade acumulada
# Incluir método estátisco de previsão objetiva


from mpl_toolkits.basemap import shiftgrid, interp
from PyFuncemeClimateTools import DefineGrid as dg
from netCDF4 import Dataset
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt


def _checkdims(model, obs):
    if not model.ndim == obs.ndim:
        print('Dims are not equals!')
        exit(1)

def compute_rmse(model, obs):

    """
    Os arrays de entrada devem ter a mesma dimensão
    :param model: numpy array com dados do modelo
    :param obs: numpy array com dados observados
    :return: Erro quadrático médio
    """

    _checkdims(model, obs)

    desv = (model - obs)**2

    myrmse = np.sqrt(np.nanmean(desv, axis=0))

    return myrmse


def compute_bias(model, obs):

    """
    Os arrays de entrada devem ter a mesma dimensão
    :param model: numpy array com os dados do modelo
    :param obs: numpy array com dados observados
    :return: Viés
    """

    _checkdims(model, obs)

    mybias = np.nanmean((model - obs), axis=0)

    return mybias


def compute_pearson(model, obs, **kwargs):

    """
    Retorna correlação de Pearson entre x e y.
    Todos os parametros são obrigatórios.

    model: numpy 3d array
    obs: : numpy 3d array

    - model e obs devem ter as mesma dimensões
    - a matriz deve ser: (tempo, lat, lon)
    """

    method = kwargs.pop('method', '3d')  # 1d ou 3d

    if method == '1d':

        corr, pvalue = st.pearsonr(model, obs)

        return corr

    elif method == '3d':

        _checkdims(model, obs)

        timelen = float(obs.shape[0])

        obs_mean = np.nanmean(obs, axis=0)
        obs_std = np.nanstd(obs, axis=0)

        model_mean = np.nanmean(model, axis=0)
        model_std = np.nanstd(model, axis=0)

        x1 = (obs - obs_mean)/obs_std
        y1 = (model - model_mean)/model_std

        xymult = x1 * y1

        xysum = np.nansum(xymult, axis=0)

        corr = xysum/timelen

        return corr

    else:

        print('--- Funcao ClimateStats.compute_pearson ---')
        print('--- Erro nos dados de entrada!---\nSaindo...')
        exit(1)


def compute_anomaly(model, obs):

    """
    Retorna anomalia e anomalia padronizada
    model = array 1d ou 3d (previsão)
    obs = array 1d ou 3d (climatologia)
    """

    # Média da climatologia
    y_mean = np.nanmean(obs, axis=0)

    # Desvio padrão
    y_std = np.nanstd(obs, axis=0)

    # Anomalia
    anom = model - y_mean

    # Anomalia padronizada
    anom_pad = (model - y_mean)/y_std

    return anom, anom_pad


def compute_righthind(forecast, hindcast, observation):

    """
    Esta função faz a correção do hincast usando
    a regressão linear (correção de bias). Mesmo
    método utilizado pelo Júnior.
    """

    # Calcula média e desvio padrão para todos os pontos da climatologia do modelo
    clim_mean = np.nanmean(hindcast, axis=0)
    clim_std = np.nanstd(hindcast, axis=0)

    # Calcula média e desvio padrão para todos os pontos da climatologia do
    # dados observado
    obs_mean = np.nanmean(observation, axis=0)
    obs_std = np.nanstd(observation, axis=0)

    # Calcula variável padronizada do modelo e observado
    clim_pad = (hindcast - clim_mean)/clim_std
    # obs_pad = (observation - obs_mean)/obs_std
    fcst_pad = (forecast - clim_mean)/clim_std

    # newhind é o hindcast corrigido
    newhind = clim_pad * obs_std + obs_mean

    newfcst = fcst_pad * obs_std + obs_mean

    return newhind, newfcst


def compute_setrighthind(hindcast, observation):

    """

    Esta função faz a correção do hincast usando
    a regressão linear (correção de bias). Mesmo
    método utilizado pelo Caio.

    - model e obs devem ter as mesma dimensões
    - a matriz deve ser: (tempo, lat, lon)

    """

    if hindcast.shape == observation.shape:

        lentime = float(observation.shape[0])

        # Calcula média e desvio padrão para todos os pontos da climatologia do modelo
        clim_mean = np.nanmean(hindcast, axis=0)
        clim_std = np.nanstd(hindcast, axis=0)

        # Calcula média e desvio padrão para todos os pontos da climatologia do
        # dados observado
        obs_mean = np.nanmean(observation, axis=0)
        obs_std = np.nanstd(observation, axis=0)

        # Calcula variável padronizada do modelo e observado
        clim_pad = (hindcast - clim_mean)/clim_std
        obs_pad = (observation - obs_mean)/obs_std

        # Calculando média e desvio padrão da variável
        # padronizada da climatologia do modelo e dado observado
        clim_pad_mean = int(np.nanmean(clim_pad))
        obs_pad_mean = int(np.nanmean(obs_pad))
        clim_pad_std = np.nanstd(clim_pad)
        obs_pad_std = np.nanstd(obs_pad)

        # Calcula correlação para todos
        # os pontos da grade usando a mesma metodologia
        # do Liqiang
        # corr = compute_pearson(observation, hindcast, lentime)
        corr = compute_pearson(observation, hindcast)

        # Júnior falou que não precisa de mascara
        # Ajuste para os pontos com significância estatística
        # if lentime >= 30 :
        #     corr = np.where(corr >= 0.3, corr, 0)
        # else:
        #     corr = np.where(corr >= 0.4, corr, 0)

        # Estima os parâmetros da reta de regressão
        # Método do mínimos quadrados
        b = corr * (obs_std / clim_std)
        a = obs_mean - b * clim_mean

        # Corrigindo o desvio padrão da previsão
        # fcst_std = obs_pad_std * (1 - corr**2)**0.5

        # Ver slides do Caio
        # newhind é o hindcast corrigido pela regressão
        newhind = a + b * (hindcast)

        return newhind

    else:

        print('Dados de entrada com dimensoes diferentes!\nSaindo...')
        exit()


def compute_probability(forecast, hindcast, observation):

    lentime = float(observation.shape[0])

    # Calcula média e desvio padrão para todos os pontos da climatologia do modelo
    clim_mean = np.nanmean(hindcast, axis=0)
    clim_std = np.nanstd(hindcast, axis=0)

    # Calcula média e desvio padrão para todos os pontos da climatologia do
    # dados observado
    obs_mean = np.nanmean(observation, axis=0)
    obs_std = np.nanstd(observation, axis=0)

    # Calcula variável padronizada do modelo e observado
    clim_pad = (hindcast - clim_mean)/clim_std
    obs_pad = (observation - obs_mean)/obs_std

    # Calculando média e desvio padrão da variável
    # padronizada da climatologia do modelo e dado observado
    clim_pad_mean = int(np.nanmean(clim_pad))  # igual a 0
    obs_pad_mean = int(np.nanmean(obs_pad))    # igual a 0
    clim_pad_std = np.nanstd(clim_pad)         # igual a 1
    obs_pad_std = np.nanstd(obs_pad)           # igual a 1

    # Calcula correlação para todos os pontos da grade
    # corr = compute_pearson(observation, hindcast, lentime)
    corr = compute_pearson(observation, hindcast)

    # Ajuste para os pontos com significância estatística
    if lentime >= 30 :
        corr = np.where(corr >= 0.3, corr, 0)
    else:
        corr = np.where(corr >= 0.4, corr, 0)

    # Estima os parâmetros da reta de regressão
    # Método do mínimos quadrados
    b = corr * (obs_pad_std / clim_pad_std)  # igual a corr
    a = obs_pad_mean - b * clim_pad_mean     # igual a 0

    # Corrigindo o desvio padrão da previsão
    fcst_std = obs_pad_std * (1 - corr**2)**0.5

    # Ver slides do Caio
    # fcst_signal_anom_pad -> Prev da anom padronizada (corrigida)
    fcst_signal_anom_pad = a + b * (forecast - clim_mean) / clim_std

    # usado para calcular anomalia corrigida
    bb = corr * (obs_std / clim_std)
    aa = obs_mean - bb * clim_mean
    # fcst_signal -> Prev corrigida pela reta de regressão
    fcst_signal = aa + bb * (forecast)

    # 0 e 1 pq a média e o desvio estão padronizados
    y_normal_obs = st.norm(loc=0, scale=1)
    y_normal_prev = st.norm(loc=fcst_signal_anom_pad, scale=fcst_std)

    r1 = 2.0/3
    r2 = 1.0/3

    xc_lower = y_normal_obs.isf(q=r1)
    xc_upper = y_normal_obs.isf(q=r2)

    p_lower = y_normal_prev.cdf(xc_lower)
    p_upper = y_normal_prev.cdf(xc_upper)

    # Probabilities categories
    # p_above = (1-p_upper)*100
    # p_normal = (p_upper - p_lower)*100
    # p_below = p_lower*100

    p_normal = np.around((p_upper - p_lower)*100)
    p_below = np.around(p_lower*100)
    p_above = 100 - (p_normal + p_below)

    # criando array para receber os tercis
    # result = np.ones(clim_mean.shape, dtype=int)*2

    # tercis: faixa abaixo(-1), faixa normal(0) e faixa acima(1)
    # result = np.where(obs < lim_lower, -1, result)
    # result = np.where((obs >= lim_lower) & (obs <= lim_upper), 0, result)
    # result = np.where(obs > lim_upper, 1, result)

    return p_below, p_normal, p_above, fcst_signal_anom_pad, \
    fcst_std, obs_pad, fcst_signal


def tercil_verification(obscli, obs):

    """
    :param obsclim: 3d numpy array climatology
    :param obs: 2d numpy array target
    :return tercilobs: 2d array tercis (-1: below) (0: normal) (1: above) (2: something is wrong)
    :return lim_lower: 2d array limite inferior
    :return lim_upper: 2d array limite superior
    """

    media = np.nanmean(obscli, axis=0)
    desvio = np.std(obscli, axis=0)

    # A curva da normal é a mesma usada na previsão
    y_normal_obs = st.norm(loc=media, scale=desvio)

    r1 = 2./3.
    r2 = 1./3.

    # limites da curva da normal
    lim_lower = y_normal_obs.isf(q=r1)
    lim_upper = y_normal_obs.isf(q=r2)

    # criando array para receber os tercis
    tercilobs = np.ones(media.shape, dtype=int)*2

    # tercis: faixa abaixo(-1), faixa normal(0) e faixa acima(1)
    tercilobs = np.where(obs < lim_lower, -1, tercilobs)
    tercilobs = np.where((obs >= lim_lower) & (obs <= lim_upper), 0, tercilobs)
    tercilobs = np.where(obs > lim_upper, 1, tercilobs)

    return tercilobs, lim_lower, lim_upper


def compute_lhss(fcst, hind, obs, climobs):

    """
    Calcula Likelihood Score (lhs) e Likelihood Skill Score (lhss).
    obsclim e hind devem ter as mesma dimensoes
    obs e clim devem ter as mesma dimensoes

    :param obsclim: - Climatologia observada. Ex. 30 anos do cmap. Array 3d (Tempo, Lat, Lon).
    :type obsclim: numpy array 3d
    :param obs:     - Dados observados. Array 3d (Tempo, Lat, Lon).
    :param hind:    - Hindcast do modelo. Ex: 30 anos de previsao do modelo. Array 3d (Tempo, Lat, Lon).
    :param fcst:    - Previsoes do modelo. Array 3d (Tempo, Lat, Lon)
    :returns: lhs, lhss

    """

    nyears = climobs[:, 0, 0].shape
    nyears = float(nyears[0])

    nlat = fcst[0, :, 0].shape
    nlat = int(nlat[0])

    nlon = fcst[0, 0, :].shape
    nlon = int(nlon[0])

    obstercis, tercilinf, tercilupp = tercil_verification(climobs, obs)

    below, normal, above, f_signal, f_std, o_pad, fcst_sig_anom = \
        compute_probability(fcst, hind, climobs)

    a = np.where(obstercis == 1,  above, 0)
    n = np.where(obstercis == 0,  normal, 0)
    b = np.where(obstercis == -1, below, 0)

    d = (a + n + b)/100.

    lhs = np.exp(np.nanmean(np.log(d), axis=0))

    lhsclim = np.ones((nlat, nlon)) * 1./3.

    lhss = (lhs - lhsclim) / (1 - lhsclim)

    return lhs, lhss


def compute_rpss(fcst, hind, obs, climobs):

    """
    Calcula Ranked Probability Score (RPS) e Ranked Probability Skill Score (RPSS)
    pbelow: Prob abaixo calculada pelo metódo do Caio. Array 2d (Lat, Lon)
    pnormal: Prob normal calculada pelo metódo do Caio. Array 2d (Lat, Lon)
    climobs: Climatologia observada. Ex. 30 anos do cmap. Array 3d (Tempo, Lat, Lon)
    obs: Dados observados. Array 2d (Tempo, Lat, Lon)
    return: rpss, rps (obs e model)
    """

    nlat = int(len(climobs[0, :, 0]))
    nlon = int(len(climobs[0, 0, :]))

    # A curva da normal é a mesma usada na previsão
    # Cálculo do percentis
    media = np.nanmean(climobs, axis=0)
    desvio = np.std(climobs, axis=0)
    y_normal_obs = st.norm(loc=media, scale=desvio)
    r1 = 2./3.
    r2 = 1./3.
    p33 = y_normal_obs.isf(q=r1)
    p66 = y_normal_obs.isf(q=r2)

    pbelow, pnormal, above, f_signal, f_std, o_pad, fcst_sig_anom = \
        compute_probability(fcst, hind, climobs)

    prob33fcst = pbelow / 100.
    prob66fcst = (pbelow + pnormal) / 100.

    prob33obs = (obs < p33)
    prob66obs = (obs < p66)

    a = ((prob33fcst[0, ...] - prob33obs)**2)/2.  # divisão 2 (número de categorias)
    b = ((prob66fcst[0, ...] - prob66obs)**2)/2.  # divisão 2 (número de categorias)

    # rps da previsão
    rps_prev = a + b

    prob33clim = np.zeros((nlat, nlon), dtype=np.float)+1/3.
    prob66clim = np.zeros((nlat, nlon), dtype=np.float)+2/3.

    a = ((prob33clim - prob33obs)**2)/2.  # divisão 2 (número de categorias)
    b = ((prob66clim - prob66obs)**2)/2.  # divisão 2 (número de categorias)

    # rps da climatologia
    rps_clim = a + b

    # rpss para cada ano da previsão
    rpss = 1 - (rps_prev/rps_clim)

    # rps médio dos anos
    rps_prev_mean = np.nanmean(rps_prev, axis=0)
    rps_clim_mean = np.nanmean(rps_clim, axis=0)
    rps_r = rps_prev_mean/rps_clim_mean
    rpss_mean = 1 - rps_r

    return rpss, rpss_mean, rps_clim, rps_prev


def ProbMemb(fcst_month, fcst_year, target_year, target_months,
             hind_period, nyears, shapef, nmemb=20):

    """
    Esta função calcula a curva de probabilida a previsão (sinal)
    para todos os membros.

    :param fcst_month: Mês previsão
    :type fcst_month: str

    :param fcst_year: Ano do mês previsão
    :type fcst_year: str

    :param target_year: Ano alvo da previsão
    :type target_year: str

    :param target_months: Trimestre alvo da previsão
    :type target_months: str

    :param hind_period: Período do hindcast
    :type hind_period: str

    :param nyears: Número de anos do hindcast
    :type nyears: int

    :param shapef: Arquivo de txt com pontos da região
    :type shapef: arquivo de txt

    :param nmemb: Número de membros da preevisão/hinscast
    :type nmemb: int

    :return sig_membs_ce: Sinal da previsão para todos os membros

    """

    #### OBSERVADO ####

    print('\n +++ LENDO A CLIMATOLOGIA OBSERVADA +++\n')

    if target_year > fcst_year:

        url = 'http://opendap2.funceme.br:8001/data/dados-obs/cmap/2.5dg/' \
              'cmap.precip.season.accum.standard.2.5dg.{0}.{1}.nc' \
              .format("1982-2011", target_months)

    else:

        url = 'http://opendap2.funceme.br:8001/data/dados-obs/cmap/2.5dg/' \
              'cmap.precip.season.accum.standard.2.5dg.{0}.{1}.nc' \
              .format("1981-2010", target_months)

    print(target_year, fcst_year)

    print(url)

    obs_data = Dataset(url, 'r')

    obs_aux = obs_data.variables['precip'][:]

    obs_lons_360 = obs_data.variables['lon'][:]

    obs_lats = obs_data.variables['lat'][:]

    obs_data.close()

    obs_aux, obs_lons = shiftgrid(180., obs_aux, obs_lons_360, start=False)

    print('\n +++ INTERPOLANDO CLIMATOLOGIA OBSERVADA +++ \n')

    # Interpolação para 1 grau

    # Nova grade de 1 grau
    newlats = np.linspace(-90, 90, 181)

    newlons = np.linspace(-180, 179, 360)

    obs = np.zeros((int(obs_aux.shape[0]), int(len(newlats)), int(len(newlons))))

    x, y = np.meshgrid(newlons, newlats)

    for i in range(0, int(obs_aux.shape[0])):

        obs[i, :, :] = interp(obs_aux[i, :, :], obs_lons, obs_lats, x, y, order=1)


    #### PREVISÃO ####

    print(' +++ LENDO OS MEMBROS PREVISÃO +++ \n')

    #~ Monta a url do ano da previsão com todos os membros
    url = 'http://opendap2.funceme.br:8001/data/dados-pos-processados-echam46' \
          '/forecasts/{0}/{2}/{1}/PCP/TRI/pcp-seasonacc-echam46-hind{0}-en%2.2d-{1}{2}_{3}{4}.ctl' \
          .format(hind_period, fcst_month, fcst_year, target_year, target_months)

    files = [url % d for d in range(51, 71)]

    # Inicializa array que irá receber os membros do ano da previsão
    fcst_t42 = np.zeros((20, 64, 128)) * 0

    for i, nc in enumerate(files):

        print(nc)

        fcst_data = Dataset(nc, 'r')

        fcst_aux = fcst_data.variables['pcp'][:]

        fcst_lons_360 = fcst_data.variables['longitude'][:]

        fcst_lats = fcst_data.variables['latitude'][:]

        fcst_data.close()

        fcst_aux, fcst_lons = shiftgrid(180., fcst_aux, fcst_lons_360, start=False)

        fcst_t42[i, :, :] = fcst_aux[0, :, :]

    print('\n    +++ INTERPOLANDO OS MEMBROS PREVISÃO +++ \n')

    # Interpolação para 1 grau

    # Nova grade de 1 grau

    fcst = np.zeros((20, int(len(newlats)), int(len(newlons)))) * 0

    for i in range(20):

        fcst[i, :, :] = interp(fcst_t42[i, :, :], fcst_lons, fcst_lats, x, y, order=1)


    #### HINDCAST ####

    print( '+++ LENDO OS MEMBROS DE CADA ANO DO HINDCAST +++\n')

    print( ' ==> AGUARDE...\n')

    # Inicializa array que irá receber todos os anos e todos os membros

    hind_t42 = np.zeros((30, 20, 64, 128)) * 0

    for i, hind_year in enumerate(range(1981, 2011)):

        if target_year > fcst_year:

            target_hind_year = int(hind_year) + 1

        else:

            target_hind_year = hind_year

        # Monta a url de cada ano com todos os membros
        url = 'http://opendap2.funceme.br:8001/data/dados-pos-processados-echam46' \
              '/hindcasts/{0}/{1}/PCP/TRI/pcp-seasonacc-echam46-hind{0}-en%2.2d-{1}{3}_{4}{2}.ctl' \
              .format(hind_period, fcst_month, target_months, hind_year, target_hind_year)

        files = [url % d for d in range(51, 71)]

        # Inicializa array que irá receber os membros para cada ano
        hindmemb = np.zeros((20, 64, 128)) * 0

        for j, nc in enumerate(files):

            print(nc)

            hind_data = Dataset(nc, 'r')

            hind_aux = hind_data.variables['pcp'][:]

            hind_lons_360 = hind_data.variables['longitude'][:]

            hind_lats = hind_data.variables['latitude'][:]

            hind_data.close()

            hind_aux, hind_lons = shiftgrid(180., hind_aux, hind_lons_360, start=False)

            hindmemb[j, :, :] = hind_aux[0, :, :]

        hind_t42[i, :, :, :] = hindmemb[:, :, :]

    print('\n +++ INTERPOLANDO OS ANOS DE CADA MEMBRO DO HINDCAST +++ \n')

    print(' ==> AGUARDE... \n')

    # Interpolação para 1 grau

    hind = np.zeros((30 , 20, int(len(newlats)), int(len(newlons)))) * 0

    for i in range(20):

        for j in range(30):

            hind[j, i, :, :] = interp(hind_t42[j, i, :, :], hind_lons, hind_lats, x, y, order=1)


    print(' +++ APLICANDO MÁSCARA DA REGIÃO +++ \n')

    # Retorna matriz com os pontos sobre o Ceará
    pointsgrid, lonlatgrid, mymatriz = dg.pointinside(newlats,
        newlons, shapefile=shapef)

    ce_fcst = np.ma.array(fcst, mask=np.tile(mymatriz, (fcst.shape[0], 1)))

    ce_hind = np.ma.array(hind, mask=np.tile(mymatriz, (hind.shape[0],
        hind.shape[1], 1)))

    ce_obs = np.ma.array(obs, mask=np.tile(mymatriz, (obs.shape[0], 1)))


    print(' +++ MÉDIA DOS PONTOS SOBRE A REGIÃO PARA TODOS OS MEMBROS +++ \n')

    ave_ce_fcst = np.zeros((int(nmemb), 1)) * 0

    ave_ce_hind = np.zeros((int(nyears), int(nmemb), 1)) * 0

    ave_ce_obs = np.zeros((int(nyears), 1)) * 0

    for i in range(int(nmemb)):
        ave_ce_fcst[i, 0] = ce_fcst[i, :, :].mean()
        for j in range(int(nyears)):
            ave_ce_hind[j, i, 0] = ce_hind[j, i, :, :].mean()

    for i in range(int(nyears)):
            ave_ce_obs[i, 0] = ce_obs[i, :, :].mean()


    print(' +++ CALCULANDO SINAL DA PREVISÃO PARA TODOS OS MEMBROS +++\n')

    sig_membs_ce = np.empty((int(nmemb)))

    sig_membs_ce[:] = np.nan

    for i in range(nmemb):
        below_ce, normal_ce, above_ce, sig_membs_ce[i], f_std_ce, \
        o_pad_ce, fcst_sig_anom = compute_probability(ave_ce_fcst[i, 0],
        ave_ce_hind[:, i, 0], ave_ce_obs[:, 0], nyears)

    return sig_membs_ce


def compute_exceed(fcst, hind, obs, **kwargs):
    '''
    Função para calcular a probabilidade de excedencia
    de um determinar limiar (thr). TAmbém faz o plot para um ponto.

    A função utiliza uma curva normal para calcular as probabilidades.

    A função funciona para matrizes 1d e 3d.

    Use thr ou perc.

    :param thr: limiar em mm - Ex: 300, 400...
    :param perc: percentual - Ex: 0,8, 0.75...
    :param fcst: dado da previsão - 2d (lat, lon)
    :param hind: dados do hindcast - 3d (time, lat, lon)
    :param obs: dados do observado - 3d (time, lat, lon)
    :param sv: Plotar e salvar as curvas cdf e pdf de uma série de dados.
               Usar somente para matrizes com uma dimensão (1d).
               Usar somente para thr.

    :return: probabilidades de excedencia de um determinado limiar
    '''

    thr = kwargs.pop('thr', False)
    perc = kwargs.pop('perc', False)
    sv = kwargs.pop('sv', False)  # Somente 1d
    fgnpdf = kwargs.pop('fgnpdf', False)  # se sv True
    fgncdf = kwargs.pop('fgncdf', False)  # se sv True
    climobs = kwargs.pop('climobs', False)
    pltclimobs = kwargs.pop('pltclimobs', False)
    pltpdftitle = kwargs.pop('pltpdftitle', '')
    pltcdftitle = kwargs.pop('pltcdftitle', '')

    if sv:
        corr = compute_pearson(hind, obs, method='1d')
    else:
        corr = compute_pearson(hind, obs)

    # reta de regressão
    b = corr * (np.nanstd(obs, axis=0) / np.nanstd(hind, axis=0))
    a = np.nanmean(obs, axis=0) - b * np.nanmean(hind, axis=0)
    # valor corrigido pela reta de regressão
    fsignal = a + b * fcst
    # desvio correlação
    fstd = np.nanstd(obs, axis=0) * (1 - corr**2)**0.5

    if thr:
        y_prev = st.norm.sf(thr, loc=fsignal, scale=fstd)

    if perc:
        y_prev = st.norm.isf(q=perc, loc=fsignal, scale=fstd)

    if pltclimobs:
        omean = np.nanmean(obs, axis=0)
        ostd = np.nanmean(obs, axis=0)
        yo = st.norm(loc=omean, scale=ostd)

    # Usado somente para um ponto
    # Somente para thr para verificação
    if sv:

        y = st.norm(loc=fsignal, scale=fstd)

        # plota curva pdf
        xO1 = np.linspace(-500, 1000)
        plt.plot(xO1, y.pdf(xO1), color='Blue', lw=2.0, label='Modelo')
        if pltclimobs:
            plt.plot(xO1, yo.pdf(xO1), color='Red', lw=2.0,
                     ls='--', label='Observado')
        plt.grid()
        plt.xlabel('(mm)')
        plt.ylabel('PROBABILIDADES')
        plt.title(pltpdftitle)
        # plt.xticks(np.arange(-600., 1000., 100))
        leg = plt.legend(loc='upper left')
        leg.get_frame().set_linewidth(0.0)
        plt.savefig(fgnpdf)
        plt.close()

        # plot curva cdf
        xO1 = np.linspace(-500, 1000)
        plt.plot(xO1, y.cdf(xO1), color='Blue', lw=2.0, label='Modelo')
        if pltclimobs:
            plt.plot(xO1, yo.cdf(xO1), color='Red', lw=2.0,
                     ls='--', label='Observado')
        plt.grid()
        plt.title(pltcdftitle)
        plt.xlabel('(mm)')
        plt.ylabel('PROBABILIDADES')
        # plt.ylim([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.])
        # plt.xticks(np.arange(-600., 1000., 100))
        plt.yticks(np.arange(0, 1.1, 0.1))
        leg = plt.legend(loc='upper left')
        leg.get_frame().set_linewidth(0.0)
        plt.savefig(fgncdf)
        plt.close()

    return y_prev


# def compute_rpss(climobs, fcst, obs):
#
#     '''
#     Calcula RPS e RPSS (Mesma metodologia usada pelo Dirceu Reis Jr (usando membros))
#     climobs = Climatologia observada. Ex. 30 anos do cmap. Array 3d (Tempo, Lat, Lon)
#     fcst = Membros das previsões a serem avaliadas. Array 4d (Membros, Tempo, Lat, Lon)
#     obs = Dados observados (Tempo, Lat, Lon)
#     fcst e obs devem ter a mesma dimensão!
#     :return: rps, rpss
#     '''
#
#     # A curva da normal é a mesma usada na previsão
#     # Cálculo do percentis
#     # nobs = float(len(obs[:, 0, 0]))  # Qtd de anos da climatologia do obs
#     media = np.nanmean(climobs, axis=0)
#     desvio = np.std(climobs, axis=0)
#     y_normal_obs = st.norm(loc=media, scale=desvio)
#     r1 = 2./3.
#     r2 = 1./3.
#     p33 = y_normal_obs.isf(q=r1)
#     p66 = y_normal_obs.isf(q=r2)
#
#     # RPS previsão
#     # Qtd de anos
#     nmemb = float(len(fcst[:, 0, 0, 0]))
#
#     # Metodo do Dirceu
#     prob33fcst = np.sum((fcst < p33), axis=0)/nmemb
#     prob66fcst = np.sum((fcst < p66), axis=0)/nmemb
#
#     prob33obs = (obs < p33)
#     prob66obs = (obs < p66)
#
#     a = ((prob33fcst - prob33obs)**2)/2.
#     b = ((prob66fcst - prob66obs)**2)/2.
#     rps_prev = a + b
#
#     prob33clim = np.zeros((2, 64, 128), dtype=np.float)+1/3.
#     prob66clim = np.zeros((2, 64, 128), dtype=np.float)+2/3.
#     a = (prob33clim - prob33obs)**2
#     b = (prob66clim - prob66obs)**2
#     rps_clim = a + b
#
#     rpsprev = np.nanmean(rps_prev, axis=0)
#
#     rpsclim = np.nanmean(rps_clim, axis=0)
#
#     rpss = 1 - (rpsprev/rpsclim)
#
#     return rpss
