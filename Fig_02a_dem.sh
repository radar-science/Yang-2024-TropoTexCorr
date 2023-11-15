#!/bin/sh
# GMT, version 6.4

gmt set FONT_LABEL 12p
gmt set FONT_ANNOT_PRIMARY 12p
gmt set FORMAT_GEO_MAP ddd.xF
gmt set MAP_FRAME_TYPE inside
gmt set MAP_FRAME_PEN thick,black
gmt set MAP_TICK_PEN thin,black

cd ./DEM
FILE='dem_km.grd'
OUTFILE='dem.ps'

################### DEM ###################
#gmt makecpt -Cgray -T-3000/3000/1 -V > topo.cpt
gmt makecpt -Cwiki-2.0.cpt -T-0.650/1.850/0.01 -V > topo.cpt
gmt grdgradient $FILE -Nt -A0 -fg -Gtopo_i.nc

#gmt grdimage $FILE -R130.72/130.99/31.84/32.06 -JM5i -BwEsN -B0.1 \
#gmt grdimage $FILE -R130.7/131.0/31.8/32.1 -JM3i -BwEsN -B0.1 \
gmt grdimage $FILE -R130.82/130.92/31.88/31.98 -JM3i -BwEsN -B0.1 \
    -Ctopo.cpt -Itopo_i.nc -P -V -K > $OUTFILE
gmt pscoast -R -J -LjBR+c31.9+o0.8i/0.3i+w2k+l+f -W -V -O -K >> $OUTFILE        #scale-bar
#gmt pscoast -R -J -LjBR+c31.9+v+o1.8i/0.3i+w3k+l+f -W -V -O -K >> $OUTFILE     #scale-bar

gmt set MAP_FRAME_TYPE plain
gmt psscale -R -J -DjBL+w0.8i/0.1i+v+o0.2i/1.8i -Ctopo.cpt -By1+lkm \
    -Bx0.5 -G0.3/1.85 -V -O >> $OUTFILE                                         #color-map


################### Convert to TIFF ###################
#gmt psconvert -E600 -Tt $OUTFILE
gmt psconvert -A0.2c -E600 -Tt $OUTFILE
