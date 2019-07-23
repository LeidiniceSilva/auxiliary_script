;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
							   ;;;;
;;;  libs
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
							   ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;  LEITURA DE DADOS E DEFIÇÃO COMO VARIAVEL  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
												 ;;;;

NOME1	="PR_INDICE_DE_CONFIANCA_MENSAL_NEB"
ARQUIVO1= systemfunc("sed 's/  / /g' "+NOME1+".txt")
VALOR1	= stringtofloat(str_get_field(ARQUIVO1,  1, " "))

NOME2	="PR_INDICE_DE_CONFIANCA_SAZONAL_NEB"
ARQUIVO2= systemfunc("sed 's/  / /g' "+NOME2+".txt")
VALOR2	= stringtofloat(str_get_field(ARQUIVO2,  1, " "))

NOME3	="PR_INDICE_DE_CONFIANCA_ANUAL_NEB"
ARQUIVO3= systemfunc("sed 's/  / /g' "+NOME3+".txt")
VALOR3	= stringtofloat(str_get_field(ARQUIVO3,  1, " "))

;;;
;;; definindo os valores que serão plotados
;;;
Cindex     	= new((/3,dimsizes(VALOR1)/),float)			
Cindex(0,0:36)	= VALOR1(7:43)
Cindex(0,  37)	= VALOR1(  45)
Cindex(1,0:36)	= VALOR2(7:43)
Cindex(1,  37)	= VALOR2(  45)
Cindex(2,0:36)	= VALOR3(7:43)
Cindex(2,  37)	= VALOR3(  45)

;;;
;;; ordenando os valores
;;;
SortCindex     	= new((/3,(38)/),float)	
	
SortCindex1	= Cindex(0,0:37)
qsort(SortCindex1)
SortCindex(0,:)	= SortCindex1

SortCindex2	= Cindex(1,0:37)
qsort(SortCindex2)
SortCindex(1,:)	= SortCindex2

SortCindex3	= Cindex(2,0:37)
qsort(SortCindex3)
SortCindex(2,:)	= SortCindex3

;;;
;;; localizando os valores ordenados
;;;
SortInd     	= new((/3,(38)/),integer)			

do caso=0,2
	do i=0,37
		conta=0
		do n=0,i
			if ( SortCindex(caso,i).eq.SortCindex(caso,n) ) then
			conta=conta+1
			end if
		end do
		SInd = ind(Cindex(caso,:).eq.SortCindex(caso,i))
		SortInd(caso,i) = SInd(conta-1)
		delete(SInd)
	end do
end do

;;;
;;; Criando a CDF
;;;
CDFperc   	= new((/3,(38)/),float)
CDF   		= new((/1/),float)					

do caso=0,2
	CDF=0.
	do i=0,37
		CDF	   	  = (1./38.) + CDF
		CDFperc(caso,i) = CDF
	end do
end do
												 ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;  PLOTAGEM DO GRAFICO  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
												 ;;;;


;;;;;;  DEFININDO OS EIXOS PARA O GRAFICO  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
									;;;
x 	   = ispan(1,38,1)		; Criando a variavel "x" que guiará o tamanho do eixo horizontal
									;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;  DEFININDO O "WKS" E "NOME" PARA O GRAFICO  ;;;;;;;;;;;;;;;;;;;;;;;
									;;;
wks   = gsn_open_wks("ps","DIST.ACUM_DO_INDICE_DE_CONFIANCA") 
									;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;  PERSONALIZAÇÂO DO PLOT  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
									;;;


;;;;;;  PLOT_1 >  LINHAS  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
							;;;
; Config. iniciais
res_1					= True					; Ativando o customizador res_1 [ se "False" o plot usaria as caracteristicas padrões ]
res_1@gsnDraw				= False					; don't draw yet		> não execute o res_1 ainda
res_1@gsnFrame				= False					; don't advance frame yet	> não salve/plot o res_1 ainda
res_1@gsnMaximize			= True					; Maximize plot
res_1@gsnPaperOrientation		= "auto"				; auto portrait or landscape
res_1@vpHeightF				= 0.40
res_1@vpWidthF				= 0.60
; Config. Titles
; X
cedil					= "c~H-13F35~K~FH2~"			; ç
atilde					= "a~H-13V2F35~D~FV-2H3~"		; ã
Iacute					= "I~H-08V6F35~B~FV-6H3~"		; Í
res_1@tiXAxisString			= Iacute+"ndice de Confian"+cedil+"a"
res_1@tiXAxisFontHeightF		= 0.009
res_1@tiXAxisFontColor			= "gray50"
res_1@tiXAxisOffsetXF			= 0.27
res_1@tiXAxisOffsetYF			= 0.065
; Y
res_1@tiYAxisString			= "Distribui"+cedil+atilde+"o Acumulada"
res_1@tiYAxisFontHeightF		= 0.009
res_1@tiYAxisFontColor			= "gray50"
res_1@tiYAxisOffsetXF			= 0.09
res_1@tiYAxisOffsetYF			= 0.1
; Config. X axis
res_1@tmXBMode				= "Explicit" 
res_1@tmXBValues			= (/-0.1,  0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8/)	; major tm
res_1@tmXBLabels			= (/"-0,1","0","0,1","0,2","0,3","0,4","0,5","0,6","0,7","0,8"/)
res_1@tmXBMajorLengthF			= 0.005					; default=0.02 for a vpHeightF=0.6
res_1@tmXBLabelFontThicknessF		= 5.5
res_1@tmXBLabelFontHeightF		= 0.008
res_1@tmXBMajorOutwardLengthF		= -0.00
res_1@tmXBLabelDeltaF			= -0.30					; push the XB labels down from plot
res_1@tmXBMinorValues			= (/-0.05,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75/)			; 
res_1@tmXBMinorLengthF			= 0.003					; Comprimento dos ticks
res_1@tmXBMinorOutwardLengthF		= -0.00					; Posisão dos ticks, colocando o sinal de menos os ticks ficaram do lado de dentro
res_1@tmXMajorGrid      		= True					; "Linhas de grade" > habilitando
res_1@tmXMajorGridLineColor		= "black"				; "Linhas de grade" > cor
res_1@tmXMajorGridLineDashPattern	= 2					; "Linhas de grade" > tipo de linha
res_1@tmXMajorGridThicknessF		= 0.05					; "Linhas de grade" > largura da linha
res_1@trXMinF				= -0.1					; x minimum value
res_1@trXMaxF				= 0.8					; x maximum value
; Config. Y axis
res_1@tmYLMode				= "Explicit" 
res_1@tmYLValues			= (/0.00,  0.2,  0.4,  0.6,  0.8,  1.0  /) ; major tm
res_1@tmYLLabels			= (/ "0%","20%","40%","60%","80%","100%" /)
res_1@tmYLMajorLengthF			= res_1@tmXBMajorLengthF
res_1@tmYLMajorOutwardLengthF		= -0.00
res_1@tmYLLabelFontHeightF		= res_1@tmXBLabelFontHeightF
res_1@tmYLLabelDeltaF			= -0.30					; push the XB labels down from plot
res_1@tmYLMinorValues			= (/0.1,0.3,0.5,0.7,0.9/)		 
res_1@tmYLMinorLengthF			= 0.003					; Comprimento dos ticks
res_1@tmYLMinorOutwardLengthF		= -0.00					; Posisão dos ticks, colocando o sinal de menos os ticks ficaram do lado de dentro
res_1@tmYMajorGrid     			= True					; "Linhas de grade" > habilitando
res_1@tmYMajorGridLineColor		= "black"				; "Linhas de grade" > cor
res_1@tmYMajorGridLineDashPattern	= 2					; "Linhas de grade" > tipo de linha
res_1@tmYMajorGridThicknessF		= 0.05					; "Linhas de grade" > largura da linha
res_1@trYMinF				= 0.					; y minimum value
res_1@trYMaxF				= 1.1					; y maximum value
; Config. Linhas a serem plotadas
res_1@xyLineColors			= (/"black","red","blue"/)		; cor das linhas
res_1@xyLineThicknesses			= (/2,2,2/)				; largura das linhas
res_1@xyDashPatterns	 		= (/0,0,0/)				; make all lines solid
res_1@xyMonoLineColor			= False					; want colored lines	; "True" todas linhas ficam em preto  ; "False" habilita usar linhas em cores diferentes
; Config. Legenda
res_1@xyExplicitLegendLabels		= (/" Mensal"," Sazonal "," Anual"/)			; texto da legenda ; os espaços extra foram para ajeitar a legenda
res_1@pmLegendDisplayMode		= "Always"				; habilita a legenda
res_1@pmLegendSide			= "Bottom"				; lugar da legenda (top,Buttom,Right,Left)
res_1@pmLegendParallelPosF		= 0.9					; move a legenda para direita
res_1@pmLegendOrthogonalPosF		= -0.27					; move a legenda para baixo
res_1@pmLegendWidthF			= 0.04					; largura da legenda
res_1@pmLegendHeightF			= 0.05					; altura da legenda
res_1@lgLabelFontHeightF		= 0.006					; tamanho do texto da legenda
res_1@lgOrientation			= "Vertical"				; Faz a legenda ficar disposta na horizonal, mas as linhas ficam verticais (bug) ; defalt é "Vertical"
res_1@lgPerimOn				= True					; habilita a caixa de legenda

; Criando o plot
PLOT_1					= gsn_csm_xy(wks,SortCindex(0:2,:),CDFperc(0:2,:),res_1)

							;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;  PLOT_2 e PLOT_3 > MARKER AND NUMBER  ;;;;;;;;;;;;;
							;;;

res_2					= True
res_2@gsMarkerThicknessF		= 1.					; default=1.0
res_2@gsMarkerSizeF			= 0.014					; Default: 0.007 
res_2@gsMarkerOpacityF			= 0.85					; Opacity of marker
res_2@gsMarkerIndex			= 16					; marker style (+)
	
res_3					= True					; text options for points
res_3@txJust				= "BottomCenter"			; Default="CenterCenter".
res_3@txFontThicknessF			= 1.5					; default=1.00
res_3@txFontHeightF			= 0.0068				; default=0.05
res_3@txFontColor			= "black"				; default=0.05
      markerTxYOffset			= -0.0064				; offset for text above marker

PLOT_2 = new((3*38),graphic)
PLOT_3 = new((3*38),graphic)

N_maker	= (/ "1", "2", "3", "4", "5", "7","10","11", "6","12", \
	    "13","14","15","16","17","18","19","20","21","22", \
	    "23","24","25","27","28","26","31","32", "8", "9", \
	    "29","30","33","34","35","36","37","M" /)

N_num	= (/ 1, 2, 3, 4, 5, 7,10,11, 6,12, \
	    13,14,15,16,17,18,19,20,21,22, \
	    23,24,25,27,28,26,31,32, 8, 9, \
	    29,30,33,34,35,36,37,38 /)

conta=0
do caso=0,2
  do i=0,37

	NmakerTx  = N_maker(SortInd(caso,i))
	Nnum	  =   N_num(SortInd(caso,i))

	if (Nnum.le.32) then
	res_2@gsMarkerColor	= "gray"					
	res_3@txFontColor	= "black"				
	end if
	if ((Nnum.ge.33).and.(Nnum.le.35)) then
	res_2@gsMarkerColor	= "goldenrod"				
	res_3@txFontColor	= "black"				
	end if
	if (Nnum.eq.36) then
	res_2@gsMarkerColor	= "deepskyblue3"			
	res_3@txFontColor	= "black"				
	end if
	if (Nnum.eq.37) then
	res_2@gsMarkerColor	= "chartreuse4"				
	res_3@txFontColor	= "black"				
	end if
	if (Nnum.eq.38) then
	res_2@gsMarkerColor	= "black"				
	res_3@txFontColor	= "white"				
	end if

	if (caso.eq.0) then
	PLOT_2(conta) = gsn_add_polymarker(wks,PLOT_1,SortCindex(caso,i)-0.00,CDFperc(caso,i),res_2) 				; Plot creation
	PLOT_3(conta) = gsn_add_text      (wks,PLOT_1,NmakerTx,SortCindex(caso,i)-0.00,CDFperc(caso,i)+markerTxYOffset,res_3)	; Plot creation
	end if

	if (caso.eq.1) then
	PLOT_2(conta) = gsn_add_polymarker(wks,PLOT_1,SortCindex(caso,i)+0.00,CDFperc(caso,i),res_2) 				; Plot creation
	PLOT_3(conta) = gsn_add_text      (wks,PLOT_1,NmakerTx,SortCindex(caso,i)+0.00,CDFperc(caso,i)+markerTxYOffset,res_3)	; Plot creation
	end if

	if (caso.eq.2) then
	PLOT_2(conta) = gsn_add_polymarker(wks,PLOT_1,SortCindex(caso,i),CDFperc(caso,i),res_2) 				; Plot creation
	PLOT_3(conta) = gsn_add_text      (wks,PLOT_1,NmakerTx,SortCindex(caso,i),CDFperc(caso,i)+markerTxYOffset,res_3)	; Plot creation
	end if

	conta=conta+1

  end do
end do
							;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;  PLOT_4 e PLOT_5 > CLASSIFICACAO  ;;;;;;;;;;;;;;;;;
							;;;

rsrRes					= True
rsrRes@gsLineThicknessF			= 4.0					; line thickness
rsrRes@gsLineDashPattern		= 0					; solid line pattern
rsrRes@gsLineColor			= "gray50";"black";"gray"		; solid line pattern

res_5					= True					; text options for points
res_5@txJust				= "BottomCenter"			; Default="CenterCenter".
res_5@txFontThicknessF			= 1.5					; default=1.00
res_5@txFontHeightF			= 0.009					; default=0.05
res_5@txFontColor			= rsrRes@gsLineColor			; default=0.05

PLOT_4 = new((10),graphic)
PLOT_5 = new((10),graphic)

; Plot creation
PLOT_4(1)			= gsn_add_polyline(wks,PLOT_1,(/-0.095,0.795/),(/1.025,1.025/),rsrRes) 		; linha horizontal
rsrRes@gsLineThicknessF		= 6.0					; line thickness
PLOT_4(2)			= gsn_add_polyline(wks,PLOT_1,(/0.4,0.4/),(/1.025+0.01,1.025-0.01/),rsrRes)	; marca
PLOT_4(3)			= gsn_add_polyline(wks,PLOT_1,(/0.5,0.5/),(/1.025+0.01,1.025-0.01/),rsrRes)	; marca
PLOT_4(4)			= gsn_add_polyline(wks,PLOT_1,(/0.6,0.6/),(/1.025+0.01,1.025-0.01/),rsrRes)	; marca
PLOT_4(5)			= gsn_add_polyline(wks,PLOT_1,(/0.65,0.65/),(/1.025+0.01,1.025-0.01/),rsrRes)	; marca
PLOT_4(6)			= gsn_add_polyline(wks,PLOT_1,(/0.75,0.75/),(/1.025+0.01,1.025-0.01/),rsrRes)	; marca

aacute  = "a~H-13V2F35~B~FV-2H3~"
eacute  = "e~H-13V2F35~B~FV-2H3~"

PLOT_5(1)			= gsn_add_text	  (wks,PLOT_1,"P"+eacute+"ssimo",  0.2  ,1.04,res_5)	; tx
PLOT_5(2)			= gsn_add_text	  (wks,PLOT_1,"Ruim",              0.45 ,1.04,res_5)	; tx 
PLOT_5(3)			= gsn_add_text	  (wks,PLOT_1,"Toler"+aacute+"vel",0.55 ,1.04,res_5)	; tx 
PLOT_5(4)			= gsn_add_text	  (wks,PLOT_1,"Mediano",           0.625,1.04,res_5)	; tx 
PLOT_5(5)			= gsn_add_text	  (wks,PLOT_1,"Bom",               0.7  ,1.04,res_5)	; tx 
PLOT_5(6)			= gsn_add_text	  (wks,PLOT_1,"Muito",             0.775,1.065,res_5)	; tx 
PLOT_5(7)			= gsn_add_text	  (wks,PLOT_1,"Bom",               0.775,1.04,res_5)	; tx 

							;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;  FAZENDO OS PLOTS  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
							;;;
;draw(PLOT_1)
;frame(wks)
;pres = True							
maximize_output(wks,True)					; "True" = attributes to be applied, "False" = don't have any resources to set, or you don't want the resources applied.
							;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;  CONVERTENDO PS PARA JPG  ;;;;;;;;;;;;;;;;;;;;;;;;;
							;;;

system("convert -density 300 DIST.ACUM_DO_INDICE_DE_CONFIANCA.ps -resize 1600x1600 -rotate -90 -crop +65+130 -crop -60-130 DIST.ACUM_DO_INDICE_DE_CONFIANCA.jpg &") 

							;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

									;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


end
