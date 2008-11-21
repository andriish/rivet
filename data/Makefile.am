dist_pkgdata_DATA = \
  refdata/ALEPH_1991_S2435284.aida       \
  refdata/ALEPH_1996_S3486095.aida       \
  refdata/DELPHI_1994_S3021912.aida      \
  refdata/DELPHI_1995_S3137023.aida      \
  refdata/DELPHI_1996_S3430090.aida      \
  refdata/DELPHI_2002_069_CONF_603.aida  \
  refdata/DELPHI_2003_WUD_03_11.aida     \
  refdata/OPAL_2004_S6132243.aida        \
  refdata/OPAL_1998_S3780481.aida        \
  refdata/CLEO_1997_S3486664.aida        \
  refdata/CLEO_2001_S4557530.aida        \
  refdata/CLEO_2004_S5809304.aida        \
  refdata/BABAR_2006_S6511112.aida       \
  refdata/BELLE_2006_S6265367.aida       \
  refdata/H1_1994_S2919893.aida          \
  refdata/H1_1995_S3167097.aida          \
  refdata/H1_2000_S4129130.aida          \
  refdata/ZEUS_2001_S4815815.aida        \
  refdata/PHENIX_2003_S5538505.aida      \
  refdata/BRAHMS_2007_S7052448.aida      \
  refdata/UA5_1989_S1926373.aida         \
  refdata/UA5_1988_S1867512.aida         \
  refdata/CDF_1994_S2952106.aida         \
  refdata/CDF_2000_S4155203.aida         \
  refdata/CDF_2001_S4751469.aida         \
  refdata/CDF_2002_S4796047.aida         \
  refdata/CDF_2004_S5839831.aida         \
  refdata/CDF_2005_S6217184.aida         \
  refdata/CDF_2006_S6653332.aida         \
  refdata/CDF_2007_S7057202.aida         \
  refdata/CDF_2008_S7541902.aida         \
  refdata/CDF_2008_S7554427.aida         \
  refdata/CDF_2008_S7540469.aida         \
  refdata/CDF_2008_S7782535.aida         \
  refdata/CDF_2008_S7828950.aida         \
  refdata/CDF_2008_NOTE_9337.aida        \
  refdata/CDF_2008_NOTE_9351.aida        \
  refdata/CDF_2008_LEADINGJETS.aida      \
  refdata/D0_2001_S4674421.aida          \
  refdata/D0_2004_S5992206.aida          \
  refdata/D0_2006_S6438750.aida          \
  refdata/D0_2007_S7075677.aida          \
  refdata/D0_2008_S6879055.aida          \
  refdata/D0_2008_S7554427.aida          \
  refdata/D0_2008_S7662670.aida          \
  refdata/D0_2008_S7719523.aida          \
  refdata/D0_2008_S7863608.aida          \
  refdata/PDG_HADRON_MULTIPLICITIES.aida \
  refdata/PDG_HADRON_MULTIPLICITIES_RATIOS.aida


EXTRA_DIST =

## bash completion
if ENABLE_PYRIVET

dist_pkgdata_DATA = rivet-completion
bashcomp_dir = $(prefix)/etc/bash_completion.d
install-data-local:
	if [[ -d "$(bashcomp_dir)" && -w "$(bashcomp_dir)" ]]; then \
      install --mode 644 rivet-completion $(bashcomp_dir)/; fi
uninstall-local:
	rm -f $(bashcomp_dir)/rivet-completion

else

EXTRA_DIST += rivet-completion

endif
