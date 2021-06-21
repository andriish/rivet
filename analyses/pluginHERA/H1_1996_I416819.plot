# ... add more histograms as you need them ...
BEGIN PLOT /H1_1996_I416819/d*
LogX=1
LogY=0
XMin=1E-5
XMax=1.
LegendYPos=0.95
LegendXPos=0.4

XLabel= $x$
YLabel= $\sigma_{red}$
#XMin=0.0005
# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /H1_1996_I416819/Q2
LogX=0
LogY=0
XMax=100
GofLegend=1
GofType=chi2
Taylor
END PLOT


BEGIN PLOT /H1_1996_I416819/d01-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 1.5$ GeV$^2$
CustomLegend= $Q^2=[1.0,2.3] $  GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d03-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 2.5 $ GeV$^2$
CustomLegend= $Q^2=[2.3,2.85] $  GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d04-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 3.5  $ GeV$^2$
CustomLegend= $Q^2=[2.85,4.3] $  GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d05-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 5 $ GeV$^2$
CustomLegend= $Q^2=[4.3,5.9] $  GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d06-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 6.5  $ GeV$^2$
CustomLegend= $Q^2=[5.9,7.2] $  GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d07-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 8.5 $ GeV$^2$
CustomLegend= $Q^2=[7.2,10.1] $  GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d08-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 12 $ GeV$^2$
CustomLegend= $Q^2=[11.5,12.5]$  GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d09-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 15 $ GeV$^2$
CustomLegend= $Q^2=[12.5,18.3]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d10-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 20 $ GeV$^2$
CustomLegend= $Q^2 = [18.3,22]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d11-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 25 $ GeV$^2$
CustomLegend= $Q^2 = [22,29]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d12-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 35 $ GeV$^2$
CustomLegend= $Q^2=[29,42]$  GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d13-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 45 $ GeV$^2$
CustomLegend= $Q^2=[42,50]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d14-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 60 $ GeV$^2$
CustomLegend= $Q^2=[50,73]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d15-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 90 $ GeV$^2$
CustomLegend= $Q^2 = [73,115]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d16-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 120 $ GeV$^2$
CustomLegend= $Q^2=[115,130]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d17-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 150 $ GeV$^2$
CustomLegend= $Q^2=[130,180]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d18-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 200 $ GeV$^2$
CustomLegend= $Q^2=[180,233]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d19-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 250 $ GeV$^2$
CustomLegend= $Q^2=[233,280]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d20-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 350 $ GeV$^2$
CustomLegend= $Q^2=[280,455]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d21-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 500 $ GeV$^2$
CustomLegend= $Q^2=[455,558]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d22-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 650 $ GeV$^2$
CustomLegend= $Q^2=[558,780]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d23-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 800 $ GeV$^2$
CustomLegend= $Q^2=[780,830]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d24-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 1000 $ GeV$^2$
CustomLegend= $Q^2=[830,1290]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d25-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 2000 $ GeV$^2$
CustomLegend= $Q^2=[1290,4000]$ GeV$^2$
END PLOT

BEGIN PLOT /H1_1996_I416819/d26-x01-y01
Title= H1 1996 $F_2(x,Q^2)$  $Q^2 = 5000 $ GeV$^2$
CustomLegend= $Q^2=[4000,6700]$ GeV$^2$
END PLOT


