;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
							   ;;;;
;;;  libs
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/taylor_diagram.ncl"
;load "./taylor_diagram_newOpts.ncl"
load "./taylor_diagram_second_quadrant_2.ncl"
							   ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;  LEITURA DE DADOS E DEFIÇÃO COMO VARIAVEL  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
												 ;;;;

	NOME="PR_CORREL_DESVIO-PADRAO_ESPACIAL_DA_MEDIA_NEB"


ARQUIVO	= systemfunc("sed 's/  / /g' "+NOME+".txt")

FLDCOR	= stringtofloat(str_get_field(ARQUIVO,  1, " "))
FLDSTD	= stringtofloat(str_get_field(ARQUIVO,  2, " "))
NOMES   =               str_get_field(ARQUIVO,  3, " ")

;;;;;;  RATIOS  ;;;;;;;;;;;;;	razão entre o desvio padrão do modelo e da observação (NORMALIZANDO)
RATIOS     	= new((/1,38/),float)			; 2 --> MGR e Had ; 1 --> um valor de ratio para cada, em relação ao CRU
RATIOS(0,0:36)	= FLDSTD(7:43)/FLDSTD(3)
RATIOS(0,  37)	= FLDSTD(  45)/FLDSTD(3)
print(max(RATIOS))
;;;;;;  CORREL  ;;;;;;;;;;;;;
CORREL     	= new((/1,38/),float)			; 2 --> MGR e Had ; 1 --> um valor de correl para cada, em relação ao CRU
CORREL(0,0:36)	= FLDCOR(7:43)
CORREL(0,  37)	= FLDCOR(  45)
												 ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;  DEFININDO O "WKS" E "NOME" PARA O GRAFICO  ;;;;;;;;;;;;;;;;;;;;;;;
									;;;
wks   = gsn_open_wks("ps","TAYLOR_DIAGRAM_"+NOME) 
									;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;  PERSONALIZAÇÂO DO PLOT  ;;;;;;;;;;;;;;;;;;;;;;;;;;
							;;;
res				= True
;
;res@tiMainString		= "Diagrama de Taylor"		; not using title makes plot bigger 
;res@tiMainFontHeightF		= 0.0225			; tiMainString size
;
res@txFontHeightF		= 0.009;75			; Size of text in axis x, y, Correlation, and of "CorLabel"
;
atilde				= "a~H-13V2F35~D~FV-2H3~"
cedil				= "c~H-13F35~K~FH2~"
res@drawCorLabel		= "Correla"+cedil+atilde+"o"	; Title Label of Correlation arc 
res@ccRays			= (/-0.95,-0.8,-0.6,-0.3,0.3,0.6,0.8,0.95/)		; Default correlation rays lines are desable
res@ccRays_color		= "LightGray"			; Default correlation rays lines color 
res@OneX			= "REF"				; 1.00, REF, OBS, DATA, etc.
;
res@drawStnDevLabel		= "Desvio Padr"+atilde+"o (Normalizado)"	; Stn. Dev. axis Label
res@stnRad			= (/0.50,1.0,1.5,2.0,2.5,3/)		; locations of Standard Dev. Radius
;
res@centerDiffRMS		= True				; Concentric arcs of RMS difference ; True to draw  
res@centerDiffRMS_color		= "green"			; Concentric arcs of RMS difference color  
;
res@Colors			= (/"gray"/)
res@Markers			= (/16/)			; Marker Indices
res@gsMarkerThicknessF		= 1.0				; external boarder of markers
res@gsMarkerSizeF		= 0.013				; marker size
res@gsMarkerOpacityF		= 0.95				; Opacity of marker ; range of 0.0 (completely transparent) to 1.0 (fully opaque)
res@markerTxFontHeightF		= 0.0067;25 			; Size of markers number text
res@markerTxYOffset		= -0.025			; offset for text above marker
res@markerTxFontColor		= "black" 			; Color of markers number text ; Default is the same color of marker
; 
;res@caseLabels			= (/"Case A","Case B","..."/)	; Default not draw it
;res@caseLabelsFontHeightF	= 0.05				; Font size for case legend
;res@legendWidth		= 0.15				; width  of legend (NDC)
;res@legendHeight		= 0.030				; height of legend (NDC) ; space for each name in legend
; 
;res@varLabels			= (/ "Var 1", "Var 2","..."/)	; Default not draw it
;res@varLabelsXloc		= 0.125				; X location of Variables legend
;res@varLabelsYloc		= 1.0				; Y location of Variables legend ; default is proportional to number of variables starting from 0.3
;res@varLabelsFontHeightF	= 0.013				; Size of text in legend list of variables
;
res@taylorDraw			= False
res@taylorFrame			= False

; Criando o PLOT
PLOT				= taylor_diagram(wks, RATIOS, CORREL, res)

; Maximizando o PLOT
pres 				= True				; "True" = attributes to be applied, "False" = don't have any resources to set, or you don't want the resources applied.
maximize_output(wks,pres)					; calls draw and frame for you / Maximize plots in frame
							;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;  CONVERTENDO PS PARA JPG  ;;;;;;;;;;;;;;;;;;;;;;;;;
							;;;

system("convert -density 300 TAYLOR_DIAGRAM_"+NOME+".ps -resize 1024x1024 -rotate -90 -crop +45+135 -crop -50-110 TAYLOR_DIAGRAM_"+NOME+".jpg &") 

							;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end
