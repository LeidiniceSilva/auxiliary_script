#!/bin/bash
#Gerar figura no ferret
#Leidinice Silva / funceme

week=4
ano=2004

echo '
use esi_($1)WK_($2)_SA.nc

fill ESI[l=22]

FRAME/FILE=esi_($1)WK_($2)_SA.gif

quit

'> plot_fig.jnl

ferret -script plot_fig.jnl $week $ano

rm plot_fig.jnl
