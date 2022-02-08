# Generation of the PSNs (then analyzed with psntools in ../psntools_analysis)

# The analyses performed in this folder have been run with PyInteraph2
# (version 1861c4861a419a7a5be69155a7ac76bea627a20f, 13 Jan 2022)
# on the following systems (path in this folder -> name in the OSF repository
# associated with the publication (https://osf.io/9v6wz)

- free/2lpc_1-169/replicate1/CHARMM22star_TIP3P -> free_replicate1_charmm22star_tip3p
- free/2lpc_1-169/replicate3/CHARMM22star_TIP3P -> free_replicate3_charmm22star_tip3p 
- beclin1/2PON_1-156_21-43/replicate1/CHARMM22star_TIP3P -> bclxl-beclin1_replicate1_charmm22star_tip3p
- beclin1/2PON_1-156_21-43/replicate2/CHARMM22star_TIP3P -> bclxl-beclin1_replicate2_charmm22star_tip3p

# In acPSN analyses, the occurrence cut-off is given as a frequency and not
# as a percentage because, up to the PyInteraph2 version used here, frequencies
# were used as cut-offs for acPSNs. If you use newer versions of the software,
# please multiply the cut-off by 100 to obtain the correct cut-off.
