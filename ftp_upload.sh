#!/bin/bash

HOST='ftp.funceme.br'
USER='cptec2'
PASSWD='20140805'

cd /home/junior/funceme/tsub/data/godas/fig/atsub_mean

data1=$(ls *AVG4PENT*.png | tail -1)

cd /home/junior/funceme/tsub/data/godas/fig/hovmoller

data2=$(ls *EVOL*.png | tail -1)

ftp -n $HOST <<EOF
user $USER $PASSWD
bin
prom
pass
cd /marcelorodrigues/enzo/mon_semanal/tsub
mdelete *AVG4PENT*.png
mdelete *ATSUB_EVOL*.png

lcd /home/junior/funceme/tsub/data/godas/fig/atsub_mean
mput ${data1}

lcd /home/junior/funceme/tsub/data/godas/fig/hovmoller
mput ${data2}

bye

EOF

cd /home/junior/funceme/tsub/data/godas/fig/atsub

data3=$(ls *ANIM15*.mp4 | tail -1)

HOST2='baco2.funceme.br'
USER2='ftpguest'
PASSWD2='Meef8ol{'

ftp -n $HOST2 <<EOF
user $USER2 $PASSWD2
bin
prom
pass
cd /modelagem/sst_monitoring
ls
pwd

mdelete *_PT_*.mp4

lcd /home/junior/funceme/tsub/data/godas/fig/atsub
mput ${data3}
ls
bye
EOF

exit 0












































