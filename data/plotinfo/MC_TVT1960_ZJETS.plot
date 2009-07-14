# BEGIN PLOT /MC_TVT1960_ZJETS/Z_mass
Title=Z mass
XLabel=$m_{\text{Z}}$ [GeV]
YLabel=$1/\sigma \text{d}\sigma/\text{d}m_{\text{Z}}$
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/dR_jet2_jet3
Title=
XLabel=$|\Delta{R}(\text{2nd jet, 3rd jet})|$
YLabel=$1/\sigma \text{d}\sigma/\text{d}|\Delta{R}(\text{2nd jet, 3rd jet})|$
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/deta_Z_jet2
Title=
XLabel=$|\Delta{\eta}(\text{Z, 1st jet})|$
YLabel=$1/\sigma \text{d}\sigma/\text{d}|\Delta{\eta}(\text{Z, 1st jet})|$
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/jet10_multi_exclusive
Title=Exclusive jet multiplicity
XLabel=$N_{\text{jet(}\geq 10 \text{GeV)}}$
YLabel=$\sigma(N_{\text{jet}})/\sigma(N_{\text{jet}}=0)$
XMajorTickMarks=1
XMinorTickMarks=0
ErrorBands=1
LegendXPos=1.15
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/jet10_multi_inclusive
Title=Inclusive jet multiplicity
XLabel=$N_{\text{jet(}\geq 10 \text{GeV)}}$
YLabel=$\sigma(\geq N_{\text{jet}})/\sigma(\text{inclusive})$
XMajorTickMarks=1
XMinorTickMarks=0
ErrorBands=1
LegendXPos=1.15
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/jet10_multi_ratio
Title=Ratio of jet multiplicity
XLabel=$N_{\text{jet(}\geq 10 \text{GeV)}}$
YLabel=$\sigma(\geq N_{\text{jet}})/\sigma(\geq N_{\text{jet}}-1)$
XMajorTickMarks=10
XMinorTickMarks=0
LogY=0
ErrorBands=1
LegendXPos=0.5
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/jet1_pT
Title=pT of 1st jet
XLabel=$p_{\perp}^{\text{1st jet}}$ [GeV]
YLabel=$1/\sigma \text{d}\sigma/\text{d}p_{\perp}^{\text{1st jet}}$
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/jet20_multi_exclusive
Title=Exclusive jet multiplicity
XLabel=$N_{\text{jet(}\geq 20 \text{GeV)}}$
YLabel=$\sigma(N_{\text{jet}})/\sigma(N_{\text{jet}}=0)$
XMajorTickMarks=1
XMinorTickMarks=0
ErrorBands=1
LegendXPos=1.15
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/jet20_multi_inclusive
Title=Inclusive jet multiplicity
XLabel=$N_{\text{jet(}\geq 20 \text{GeV)}}$
YLabel=$\sigma(\geq N_{\text{jet}})/\sigma(\text{inclusive})$
XMajorTickMarks=1
XMinorTickMarks=0
ErrorBands=1
LegendXPos=1.15
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/jet20_multi_ratio
Title=Ratio of jet multiplicity
XLabel=$N_{\text{jet(}\geq 20 \text{GeV)}}$
YLabel=$\sigma(\geq N_{\text{jet}})/\sigma(\geq N_{\text{jet}}-1)$
XMajorTickMarks=10
XMinorTickMarks=0
LogY=0
ErrorBands=1
LegendXPos=0.5
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/jet2_pT
Title=pT of 2nd jet
XLabel=$p_{\perp}^{\text{2nd jet}}$ [GeV]
YLabel=$1/\sigma \text{d}\sigma/\text{d}p_{\perp}^{\text{2nd jet}}$
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/jet3_pT
Title=pT of 3rd jet
XLabel=$p_{\perp}^{\text{3rd jet}}$ [GeV]
YLabel=$1/\sigma \text{d}\sigma/\text{d}p_{\perp}^{\text{3rd jet}}$
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/jet4_pT
Title=pT of 4th jet
XLabel=$p_{\perp}^{\text{4th jet}}$ [GeV]
YLabel=$1/\sigma \text{d}\sigma/\text{d}p_{\perp}^{\text{4th jet}}$
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/log10_R_0
Title=$\log_{10}$(Integrated $0$ jet rate in $k_\perp$ [GeV])
XLabel=$\log_{10}(d_{\text{cut}}/\text{GeV})$
YLabel=$R_{0}$
Rebin=2
YMin=9e-4
YMax=1.5
LegendYPos=0.8
LegendXPos=1.2
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/log10_R_1
Title=$\log_{10}$(Integrated $1$ jet rate in $k_\perp$ [GeV])
XLabel=$\log_{10}(d_{\text{cut}}/\text{GeV})$
YLabel=$R_{1}$
YMin=9e-6
YMax=0.5
Rebin=2
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/log10_R_2
Title=$\log_{10}$(Integrated $2$ jet rate in $k_\perp$ [GeV])
XLabel=$\log_{10}(d_{\text{cut}}/\text{GeV})$
YLabel=$R_{2}$
YMin=9e-6
YMax=0.5
Rebin=2
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/log10_R_3
Title=$\log_{10}$(Integrated $3$ jet rate in $k_\perp$ [GeV])
XLabel=$\log_{10}(d_{\text{cut}}/\text{GeV})$
YLabel=$R_{3}$
YMin=9e-6
YMax=0.5
Rebin=2
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/log10_R_4
Title=$\log_{10}$(Integrated $4$ jet rate in $k_\perp$ [GeV])
XLabel=$\log_{10}(d_{\text{cut}}/\text{GeV})$
YLabel=$R_{\geq4}$
YMin=9e-6
YMax=0.5
Rebin=2
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/log10_d_01
Title=$\log_{10}$($k_\perp$ jet resolution $0 \to 1$ [GeV])
XLabel=$\log_{10}(d_{01}/\text{GeV})$
YLabel=$\text{d}\sigma/\text{d}\log_{10}(d_{01})$
LegendXPos=0.5
LegendYPos=0.5
YMin=9e-6
YMax=3.0
Rebin=2
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/log10_d_12
Title=$\log_{10}$($k_\perp$ jet resolution $1 \to 2$ [GeV])
XLabel=$\log_{10}(d_{12}/\text{GeV})$
YLabel=$\text{d}\sigma/\text{d}\log_{10}(d_{12})$
LegendXPos=0.5
LegendYPos=0.5
YMin=9e-6
YMax=3.0
Rebin=2
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/log10_d_23
Title=$\log_{10}$($k_\perp$ jet resolution $2 \to 3$ [GeV])
XLabel=$\log_{10}(d_{23}/\text{GeV})$
YLabel=$\text{d}\sigma/\text{d}\log_{10}(d_{23})$
LegendXPos=0.5
LegendYPos=0.5
YMin=9e-6
YMax=3.0
Rebin=2
# END PLOT

# BEGIN PLOT /MC_TVT1960_ZJETS/log10_d_34
Title=$\log_{10}$($k_\perp$ jet resolution $3 \to 4$ [GeV])
XLabel=$\log_{10}(d_{34}/\text{GeV})$
YLabel=$\text{d}\sigma/\text{d}\log_{10}(d_{34})$
LegendXPos=0.5
LegendYPos=0.5
YMin=9e-6
YMax=3.0
Rebin=2
# END PLOT

