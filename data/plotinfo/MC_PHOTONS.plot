# BEGIN PLOT /MC_PHOTONS/Ptgamma
Title=Photon $p_\perp$
XLabel=$p_\perp^\gamma$ [GeV]
#YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}p_\perp^\gamma$ [GeV$^{-1}$]
# END PLOT

# BEGIN PLOT /MC_PHOTONS/Egamma
Title=Photon energy
XLabel=$E_\gamma$ [GeV]
#YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}E_\gamma$ [GeV$^{-1}$]
# END PLOT

# BEGIN PLOT /MC_PHOTONS/sumPtgamma
Title=Scalar sum of photon $p_\perp$
XLabel=$\sum{p_\perp^\gamma}$ [GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}\sum{p_\perp^\gamma}$ [GeV$^{-1}$]
# END PLOT

# BEGIN PLOT /MC_PHOTONS/sumEgamma
Title=Sum of photon energies
XLabel=$\sum{E_\gamma}$ [GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}\sum{E_\gamma}$ [GeV$^{-1}$]
# END PLOT

# BEGIN PLOT /MC_PHOTONS/deltaR
Title=$\Delta{R}$ between photons and associated leptons
XLabel=$\Delta{R}$
#YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}\Delta{R}$
# END PLOT

# BEGIN PLOT /MC_PHOTONS/deltaR_weighted
Title=$p_\perp$-weighted $\Delta{R}$ between photons and associated leptons
XLabel=$\Delta{R}$
#YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}\Delta{R}$
# END PLOT

# BEGIN PLOT /MC_PHOTONS/deltaR_R
Title=Phase-space weighted $\Delta{R}$ between photons and associated leptons
XLabel=$\Delta{R}$
#YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}\Delta{R}$
# END PLOT

# BEGIN PLOT /MC_PHOTONS/deltaR_R_weighted
Title=Phase-space \& $p_\perp$-weighted $\Delta{R}$ between photons and associated leptons
XLabel=$\Delta{R}$
#YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}\Delta{R}$
# END PLOT

# BEGIN PLOT /MC_PHOTONS/.*_vs_pTlep
XLabel=Lepton $p_\perp$ [GeV]
# END PLOT

# BEGIN PLOT /MC_PHOTONS/DeltaR_vs_pTlep
Title=Mean $\Delta{R}$ between photons and associated leptons vs. lepton $p_\perp$
YLabel=$\angle\Delta{R}\rangle$
# END PLOT

# BEGIN PLOT /MC_PHOTONS/DeltaR_weighted_vs_pTlep
Title=Mean $p_\perp$-weighted $\Delta{R}$ between photons and associated leptons vs. lepton $p_\perp$
YLabel=$\angle\text{Weighted}\quad\Delta{R}\rangle$
# END PLOT

# BEGIN PLOT /MC_PHOTONS/sumPtGamma_vs_pTlep
Title=Mean scalar $\sum{p_\perp}$ of photons on associated leptons vs. lepton $p_\perp$
YLabel=$\angle\sum{p_\perp}\rangle$
# END PLOT
