# -*- coding: utf-8 -*-

"""
Dictionaries of the basins and the subbasins

> {basins}

    Itens format and locations:

        N:['macro basin name',
           'macro or micro basin name',
           'type macro or micro']

    where:

        N - number for ordenation, macro basin in 1:18, micro in 19:262.
        macro basin name - name of the macro basin, for macro and micro type.
        macro of micro basin name - name of the macro of micro basin in situ.
        smap parameter 1 - smap parameter 1
        smap parameter 2 - smap parameter 2
        smap parameter 3 - smap parameter 3
        smap parameter 4 - smap parameter 4
        type macro or micro - self explained

    Examples:

         1:['amazonas', 'amazonas', 'macro']
        19:['amazonas', 'balbina',  'micro'],

> {basin_n}

    Itens format and locations:

        'basin name':N

    where:

        basin name - if macro: 'macro name', if micro: 'macro_micro name'.

    Examples:

        'amazonas':1
        'amazonas_balbina':19

> {n_basin}

    Itens format and locations:

        N:'basin name'

    Exampĺes:

         1:'amazonas'
        19:'amazonas_balbina'


"""


basins_flow = {'amazonas_balbina': 269,
    'amazonas_belo_monte': 292,
    'amazonas_cachoeira_caldeirao': 204,
    'amazonas_coaracy_nunes': 280,
    'amazonas_colider': 228,
    'amazonas_curua_una': 277,
    'amazonas_dardanelos': 291,
    'amazonas_ferreira_gomes': 297,
    'amazonas_guapore': 296,
    'amazonas_jirau': 285,
    'amazonas_rondon_ii': 145,
    'amazonas_samuel': 279,
    'amazonas_santo_antonio': 287,
    'amazonas_santo_antonio_do_jari': 290,
    'amazonas_sao_manuel': 230,
    'amazonas_sinop': 227,
    'amazonas_teles_pires': 229,

    'atlantico_leste_irape': 255,
    'atlantico_leste_itapebi': 188,
    'atlantico_leste_pedra_do_cavalo': 25,
    'atlantico_leste_santa_clara': 283,

    'atlantico_sudeste_henry_borden': 318,
    'atlantico_sudeste_lajes_fontes_nova_pereira_passos': 202,
    'atlantico_sudeste_pedras': 116,
    'atlantico_sudeste_rosal': 196,

    'atlantico_sul_capivari_cachoeira': 115,
    'atlantico_sul_salto_pilao': 101,

    'doce_aimores': 148,
    'doce_baguari': 141,
    'doce_cadonga': 149,
    'doce_guilman': 262,
    'doce_mascarenhas': 144,
    'doce_porto_estrela': 263,
    'doce_sa_carvalho': 183,
    'doce_salto_grande': 134,

    'grande_agua_vermelha': 18,
    'grande_as_oliveira': 16,
    'grande_caconde': 14,
    'grande_camargos': 1,
    'grande_euclides_da_cunha': 15,
    'grande_funil_grande': 211,
    'grande_furnas': 6,
    'grande_igarapava': 10,

    'grande_itutinga': 2,
    'grande_jaguara': 9,
    'grande_lc_barreto': 8,
    'grande_marimbondo': 17,
    'grande_mascarenhas_de_moraes': 7,
    'grande_porto_colombia': 12,
    'grande_volta_grande': 11,

    'iguacu_baixo_iguacu': 81,
    'iguacu_foz_do_areia': 74,
    'iguacu_fundao': 72,
    'iguacu_gov_jose_richa': 222,
    'iguacu_jordao': 73,
    'iguacu_salto_osorio': 78,
    'iguacu_salto_santiago': 77,
    'iguacu_santa_clara': 71,
    'iguacu_segredo': 76,

    'jacui_14_de_julho': 284,
    'jacui_castro_alves': 98,
    'jacui_dona_francisca': 114,
    'jacui_ernestina': 110,
    'jacui_itauba': 113,
    'jacui_jacui': 112,
    'jacui_jacui_inc': 120,
    'jacui_monte_claro': 97,
    'jacui_passo_real': 111,

    'paraguai_itiquira': 259,
    'paraguai_jauru': 295,
    'paraguai_manso': 278,
    'paraguai_ponte_de_pedra': 281,

    'paraiba_do_sul_anta': 129,
    'paraiba_do_sul_funil': 123,
    'paraiba_do_sul_ilha_dos_pombos': 130,
    'paraiba_do_sul_jaguari': 120,
    'paraiba_do_sul_itaocara_i': 199,
    'paraiba_do_sul_paraibuna': 121,
    'paraiba_do_sul_picada': 197,
    'paraiba_do_sul_santa_branca': 54,
    'paraiba_do_sul_santa_cecilia': 125,
    'paraiba_do_sul_santana': 203,
    'paraiba_do_sul_sobragi': 198,
    'paraiba_do_sul_tocos': 201,

    'parana_ilha_solteira': 34,
    'parana_ilha_solteira_equivalente': 244,
    'parana_itaipu': 266,
    'parana_jupia': 245,
    'parana_porto_primavera': 246,

    'paranaiba_barra_dos_coqueiros': 248,
    'paranaiba_batalha': 22,
    'paranaiba_cachoeira_dourada': 32,
    'paranaiba_cacu': 247,
    'paranaiba_capim_branco_i': 207,
    'paranaiba_capim_branco_ii': 28,
    'paranaiba_corumba_i': 209,
    'paranaiba_corumba_iii': 23,
    'paranaiba_corumba_iv': 205,
    'paranaiba_emborcacao': 24,
    'paranaiba_espora': 99,
    'paranaiba_foz_do_rio_claro': 261,
    'paranaiba_itumbiara': 31,
    'paranaiba_miranda': 206,
    'paranaiba_nova_ponte': 25,
    'paranaiba_salto': 294,
    'paranaiba_salto_rio_verdinho': 241,
    'paranaiba_sao_simao': 33,
    'paranaiba_serra_do_facao': 251,

    'paranapanema_canoas_i': 52,
    'paranapanema_canoas_ii': 51,
    'paranapanema_capivara': 61,
    'paranapanema_chavantes': 49,
    'paranapanema_jurumirim': 47,
    'paranapanema_maua': 57,
    'paranapanema_ourinhos': 249,
    'paranapanema_piraju': 48,
    'paranapanema_rosana': 63,
    'paranapanema_salto_grande_l_n_garcez': 50,
    'paranapanema_taquarucu_escola_politecnica': 62,

    'parnaiba_boa_esperanca': 190,

    'sao_francisco_apolonio_sales': 173,
    'sao_francisco_complexo_paulo_afonso': 176,
    'sao_francisco_itaparica': 172,
    'sao_francisco_paulo_afonso': 175,
    'sao_francisco_queimado': 158,
    'sao_francisco_retiro_baixo': 155,
    'sao_francisco_sobradinho': 169,
    'sao_francisco_tres_marias': 156,
    'sao_francisco_xingo': 178,

    'tiete_bariri': 238,
    'tiete_barra_bonita': 273,
    'tiete_billings': 118,
    'tiete_billings_mais_pedras': 119,
    'tiete_edgard_de_souza': 161,
    'tiete_guarapiranga': 117,
    'tiete_ibitinga': 239,
    'tiete_nova_avanhandava': 242,
    'tiete_ponte_nova': 160,
    'tiete_promissao': 240,
    'tiete_traicao': 104,
    'tiete_tres_irmaos': 243,

    'tocantins_cana_brava': 191,
    'tocantins_estreito_tocantins': 271,
    'tocantins_lajeado': 273,
    'tocantins_peixe_angical': 257,
    'tocantins_sao_salvador': 253,
    'tocantins_serra_da_mesa': 270,
    'tocantins_tucurui': 275,

    'uruguai_barra_grande': 215,
    'uruguai_campos_novos': 216,
    'uruguai_foz_do_chapeco': 94,
    'uruguai_garibaldi': 89,
    'uruguai_ita': 92,
    'uruguai_machadinho': 217,
    'uruguai_monjolinho': 220,
    'uruguai_passo_fundo': 93,
    'uruguai_passo_sao_joao': 103,
    'uruguai_quebra_queixo': 286,
    'uruguai_sao_jose': 102,
    'uruguai_sao_roque': 88}
