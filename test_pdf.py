# -*- coding: utf-8 -*-

"""
This script creates PDF of accumulated and anomaly
    precipitation simulated by GFS0p25 model.
"""


import re
import os
import argparse
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib.units import inch
import datetime
import glob

today = datetime.datetime.today().strftime('%d-%m-%Y')

__author__ = 'Leidinice Silva'
__email__ = 'leidinice.silva@funceme.br'
__date__ = '25/01/2018'
__description__ = "This script creates PDF of accum and anom pr by GFS0p25 model"


def finish_page(pdf, page):
    w, h = landscape(letter)    
    pdf.setFont("Helvetica", 22)
    pdf.drawString(w - 35, 10, str(page))
    pdf.showPage()


def cover(pdf):
    w, h = landscape(letter)
    pdf.setFont("Helvetica", 40)
    pdf.drawString(150, h / 2 + 100, "Monitoramento/Previsão da")
    pdf.drawString(80, h / 2 + 60, "Anomalia e Acumulado de PR")
    pdf.drawString(100, h / 2 + 20, "para o NEB do Brail")
    pdf.setFont("Helvetica-Bold", 40)
    pdf.drawImage("images/funceme.jpg", inch, inch, width=100, height=100)
    pdf.drawImage("images/ceara.jpg", h - 50, inch - 30, width=160, height=160)
    pdf.setFont("Helvetica", 25)
    pdf.drawString(w/2-80, 150, 'Emitido em:')
    pdf.setFont("Helvetica-Bold", 30)
    pdf.drawString(w/2-95, 120, today)
    pdf.showPage()


def pages(pdf):
    w, h = landscape(letter)

    pdf.setFont("Helvetica", 22)
    pdf.drawString(30, h - 150, "À direita: Evolução temporal (eixo y) da")
    pdf.drawString(30, h - 175, "da anomalia de TSM no Pacífico")
    pdf.drawString(30, h - 200, "equatorial entre as longitudes")
    pdf.drawString(30, h - 225, "120° leste e 80° oeste (eixo x).")
    pdf.drawString(30, h - 270, "Abaixo: Área em que foi calculada")
    pdf.drawString(30, h - 295, "a média meridional (5N-5S).")
    pdf.drawImage(images[2], h - 230, 90, width=400, height=430)
    pdf.drawImage("images/{0}_EVOL_AREA.png".format(args.region), 30, 150, width=349, height=97.5)
    finish_page(pdf, 2)

    pdf.setFont("Helvetica", 18)
    #Abaixo: Média da anomalia de TSM das ultimas 4 (quatro) semanas.
    #A direita: Condições da anomalia de TSM das ultimas 4 (quatro) semanas.
    pdf.drawString(30, h - 50, "Abaixo: Média da anomalia de TSM das últimas 4 (quatro)")
    pdf.drawString(30, h - 70, "semanas.")
    pdf.drawString(30, h - 100, "À direita: Condições da anomalia de TSM das últimas 4")
    pdf.drawString(30, h - 120, "(quatro) semanas.")
    pdf.drawImage(images[0], h - 120, 0, width=302.94, height=604.12)
    pdf.drawImage(images[4], 30, inch, width=452, height=333)
    finish_page(pdf, 4)


def arguments():
    global args

    parser = argparse.ArgumentParser(prog="Script PDF-TSM", description=__description__)
    parser.add_argument("data1", help="Initialize the server in Debug mode")
    parser.add_argument("data2", help="Initialize the server in Debug mode")
    parser.add_argument("data3", help="Initialize the server in Debug mode")
    parser.add_argument("data4", help="Initialize the server in Debug mode")
    parser.add_argument("-f", "--filename", default='my_file.pdf', help="Define o nome do arquivo.")
    parser.add_argument("-v", "--version", action='version',
                        version='%(prog)s Version {0}'.format(__version__),
                        help="Shows application version.")

    args = parser.parse_args()

    args.region = args.region.upper()


if __name__ == "__main__":
    arguments()

    obs_dir = os.environ["OBS_DIR"]
    ssta_dir = os.path.join(obs_dir, "ssta")
    subta_dir = os.path.join(obs_dir, "subta")

    sst_years = [os.path.join(ssta_dir, f) for f in os.listdir(ssta_dir) if f.isdigit()]
    sst_years = sorted(sst_years, reverse=True)

    subta_years = [os.path.join(subta_dir, f) for f in os.listdir(subta_dir) if f.isdigit()]
    subta_years = sorted(subta_years)
    subta_dir = subta_years[-1]

    files = []
    for year in sst_years:
        files += [os.path.join(year, im) for im in os.listdir(year)]
        if len(files) > 50:
            break
 
    r_img = re.compile(
        "^ssta_+(diff|evol|index|mean|tend|4week|5week)+_oisst_pcf_1dg_\d+-+\d+-+\d+_1w.(png)$")
    r_txt = re.compile("^ssta_+mnindex+_oisst_pcf_1dg_\d+-+\d+-+\d+_1w.(txt)$")

    images_aux = [im for im in files if r_img.match(im.split(os.sep)[-1])]
    #images_aux = [os.path.join(data_dir, f) for f in os.listdir(data_dir) if r_img.match(f)]
    images_aux = sorted(images_aux)

    last_file = (images_aux[-1].split(os.sep))[-1]
    last_date = last_file.split("_")[5]
    images = [img for img in images_aux if last_date in img]
    avg_list = sorted(glob.glob(subta_dir + "/*_4pentmean_*.png"))
    evol_list = sorted(glob.glob(subta_dir + '/*_evol_*.png'))
    images.append(avg_list[-1])
    images.append(evol_list[-1])
    images.append(os.path.join(ssta_dir, 'figure3.png'))
    images.append(os.path.join(ssta_dir, 'figure6.png'))

    aux_dir = os.path.join(ssta_dir, "2018")
    for f in files:
        if r_txt.match(f.split(os.sep)[-1]):
            txt = f
    pdf = canvas.Canvas(args.filename, pagesize=landscape(letter))
    cover(pdf)
    pages(pdf)
    pdf.save()
