
'reinit'

'open anual.nc.ctl'

********
** VARIAVEIS
********

varu=ua
varv=va
varq=qas
varp=ps
varn=flux_uvq


********

'c'

'set mpdset brmap_hires'
'set map 1 1 4'
'set grads off'
'set grid off'
'set lat -50 15'
'set lon -90 -30'

**********


'define 'varq'z=vint(('varp'),('varu'*'varq'),100)'
'define 'varq'm=vint(('varp'),('varv'*'varq'),100)'

'set gxout shaded'
'color.gs -400 400 10 -kind navy->mediumblue->deepskyblue->white->goldenrod->chocolate->red'

'd (2592000*(hdivg(('varq'z),('varq'm))))/12'

'xbar.gs -direction vertical -400 400 -fskip 5'

*'set gxout stream'

'd skip(('varq'z),6,6);('varq'm)'
*'d 'varq'z;'varq'm)'

'printim 'varn'.png x800 y600 white'

