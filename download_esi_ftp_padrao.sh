#!/bin/bash

ANO=`date +%Y`
PATH="/home/leidinice/Documentos/esi/results/2016/"
echo
cd ${PATH}

if [ ! -d "${ANO}" ]; then
    echo "Diretorio ${ANO} nao existe"
    mkdir ${ANO}
fi

cd ${ANO}
pwd

/usr/bin/wget -N ftp://ftp.star.nesdis.noaa.gov/pub/smcd/emb/chain/bbarker/DFPPM_*WK_SMN_${ANO}???.dat
/usr/bin/wget -N ftp://ftp.star.nesdis.noaa.gov/pub/smcd/emb/chain/bbarker/DFPPM_*WK_${ANO}???.dat



skip(){
HOST2='ftp.star.nesdis.noaa.gov'
USER2='anonymous'
PASSWD2='@Funceme'
for AAA in $(/usr/bin/seq -w 008 7 365) ; do

    data1=DFPPM_12WK_${ANO}${AAA}.dat
    data2=DFPPM_12WK_SMN_${ANO}${AAA}.dat
    data3=DFPPM_4WK_${ANO}${AAA}.dat
    data4=DFPPM_4WK_SMN_${ANO}${AAA}.dat
    
    echo ${data1}
    echo ${data2}
    echo ${data3}
    echo ${data4}    

    /usr/bin/ftp -n $HOST2 << EOF

user $USER2 $PASSWD2
ls
bin
prompt
pass
cd /pub/smcd/emb/chain/bbarker/ 
pwd
get $data1
get $data2
get $data3
get $data4
quit
EOF
done
exit
}







