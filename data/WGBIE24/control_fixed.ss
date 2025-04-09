#V3.30.xx.yy;_safe;_compile_date:_Feb  1 2022;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_12.3
#_Stock_Synthesis_is_a_work_of_the_U.S._Government_and_is_not_subject_to_copyright_protection_in_the_United_States.
#_Foreign_copyrights_may_apply._See_copyright.txt_for_more_information.
#_User_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_User_info_available_at:https://vlab.noaa.gov/group/stock-synthesis
#_Source_code_at:_https://github.com/nmfs-stock-synthesis/stock-synthesis

#C file created using the SS_writectl function in the R package r4ss
#C file write time: 2021-01-27 10:37:57
#_data_and_control_files: shake_data.ss // control_fixed_growth.ss
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns (Growth Patterns, Morphs, Bio Patterns, GP are terms used interchangeably in SS3)
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Platoon_within/between_stdev_ratio (no read if N_platoons=1)
#_Cond  1 #vector_platoon_dist_(-1_in_first_val_gives_normal_approx)
#
3 # recr_dist_method for parameters:  2=main effects for GP, Area, Settle timing; 3=each Settle entity; 4=none (only when N_GP*Nsettle*pop==1)
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
2 # 3 #  number of recruitment settlement assignments 
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
 1 7 1 0
# 1 4 1 0
#
#_Cond 0 # N_movement_definitions goes here if Nareas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
1 #_Nblock_Patterns
1 #_blocks_per_pattern 
# begin and end years of blocks
 1960 1993
# 2004 2021
#
# controls for all timevary parameters 
1 #_time-vary parm bound check (1=warn relative to base parm bounds; 3=no bound check); Also see env (3) and dev (5) options to constrain with base bounds
#
# AUTOGEN
 1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen time-varying parms of this category; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
#_Available timevary codes
#_Block types: 0: P_block=P_base*exp(TVP); 1: P_block=P_base+TVP; 2: P_block=TVP; 3: P_block=P_block(-1) + TVP
#_Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
#_EnvLinks:  1: P(y)=P_base*exp(TVP*env(y));  2: P(y)=P_base+TVP*env(y);  3: P(y)=f(TVP,env_Zscore) w/ logit to stay in min-max;  4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
#_DevLinks:  1: P(y)*=exp(dev(y)*dev_se;  2: P(y)+=dev(y)*dev_se;  3: random walk;  4: zero-reverting random walk with rho;  5: like 4 with logit transform to stay in base min-max
#_DevLinks(more):  21-25 keep last dev for rest of years
#
#_Prior_codes:  0=none; 6=normal; 1=symmetric beta; 2=CASAL's beta; 3=lognormal; 4=lognormal with biascorr; 5=gamma
#
# setup for M, growth, wt-len, maturity, fecundity, (hermaphro), recr_distr, cohort_grow, (movement), (age error), (catch_mult), sex ratio 
#_NATMORT
1 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=BETA:_Maunder_link_to_maturity
4 #_N_breakpoints
 0 1 5 15 # age(real) at M breakpoints
#
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0.5 #_Age(post-settlement)_for_L1;linear growth below this
999 #_Growth_Age_for_L2 (999 to use as Linf)
-998 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0  #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
#
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
2 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
2 #_parameter_offset_approach for M, G, CV_G:  1- direct, no offset**; 2- male=fem_parm*exp(male_parm); 3: male=female*exp(parm) then old=young*exp(parm)
#_** in option 1, any male parameter with value = 0.0 and phase <0 is set equal to female parameter
#
#_growth_parms
#_ LO HI INIT PRIOR PR_SD PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# Sex: 1  BioPattern: 1  NatMort
 0.15 3 1.19 0 0 0 -3 0 0 0 0 0 0 0 # NatM_break_1_Fem_GP_1
 0.15 1.4 0.64 0 0 0 -3 0 0 0 0 0 0 0 # NatM_break_2_Fem_GP_1
 0.15 0.4 0.34 0 0 0 -3 0 0 0 0 0 0 0 # NatM_break_3_Fem_GP_1
 0.15 0.4 0.2 0 0 0 -3 0 0 0 0 0 0 0 # NatM_break_4_Fem_GP_1
# Sex: 1  BioPattern: 1  Growth
 16 23 20.5095 16.165 0 0 2 0 0 0 0 0 0 0 # L_at_Amin_Fem_GP_1
 60 140 110 0 0 0 -3 0 0 0 0 0 0 0 # L_at_Amax_Fem_GP_1
 0.05 0.4 0.14 0 0 0 -2 0 0 0 0 0 0 0 # VonBert_K_Fem_GP_1
 0.005 0.5 0.15 0 0 0 -6 0 0 0 0 0 0 0 # CV_young_Fem_GP_1
 0.005 0.5 0.15 0 0 0 -6 0 0 0 0 0 0 0 # CV_old_Fem_GP_1
# Sex: 1  BioPattern: 1  WtLen
 -1 1 0.00000377 0 0 0 -3 0 0 0 0 0 0 0 # Wtlen_1_Fem_GP_1
 2 4 3.16826 0 0 0 -3 0 0 0 0 0 0 0 # Wtlen_2_Fem_GP_1
# Sex: 1  BioPattern: 1  Maturity&Fecundity
 30 55 42.36 0 0 0 -3 0 0 0 0 0 0 0 # Mat50%_Fem_GP_1
 -1 1 -0.265 0 0 0 -3 0 0 0 0 0 0 0 # Mat_slope_Fem_GP_1
 -3 3 1 0 0 0 -3 0 0 0 0 0 0 0 # Eggs/kg_inter_Fem_GP_1
 -3 3 0 0 0 0 -3 0 0 0 0 0 0 0 # Eggs/kg_slope_wt_Fem_GP_1
# Sex: 2  BioPattern: 1  NatMort
 -0.7 0.7 0 0 0 0 -3 0 0 0 0 0 0 0 # NatM_break_1_Mal_GP_1
 -0.7 0.7 0 0 0 0 -3 0 0 0 0 0 0 0 # NatM_break_2_Mal_GP_1
 -0.7 0.7 0.163 0 0 0 -3 0 0 0 0 0 0 0 # NatM_break_3_Mal_GP_1
 -0.7 0.7 0.336 0 0 0 -3 0 0 0 0 0 0 0 # NatM_break_4_Mal_GP_1
# Sex: 2  BioPattern: 1  Growth
 -0.7 0.7 0 0 0 0 -2 0 0 0 0 0 0 0 # L_at_Amin_Mal_GP_1
 -0.6 0 -0.4 0 0 0 -4 0 0 0 0 0 0 0 # L_at_Amax_Mal_GP_1
 -0.6 0.6 0 0 0 0 -3 0 0 0 0 0 0 0 # VonBert_K_Mal_GP_1
 -0.7 0.7 0 0 0 0 -3 0 0 0 0 0 0 0 # CV_young_Mal_GP_1
 -0.7 0.7 0 0 0 0 -3 0 0 0 0 0 0 0 # CV_old_Mal_GP_1
# Sex: 2  BioPattern: 1  WtLen
 -1 1 0.00000377 0 0 0 -3 0 0 0 0 0 0 0 # Wtlen_1_Mal_GP_1
 2 4 3.16826 0 0 0 -3 0 0 0 0 0 0 0 # Wtlen_2_Mal_GP_1
# Hermaphroditism
#  Recruitment Distribution  
 -22 35 20.5897 0 0 0 -3 0 0 0 0 0 0 0 # RecrDist_GP_1_area_1_month_1
 -22 35 21.5486 -0.56 0 0 3 0 23 1998 2023 4 0 0  # RecrDist_GP_1_area_1_month_7
#  Cohort growth dev base
 0 2 1 1 0 0 -3 0 0 0 0 0 0 0 # CohortGrowDev
#  Movement
#  Age Error from parameters
#  catch multiplier
#  fraction female, by GP
 1e-06 0.999999 0.5 0.5 0 0 -99 0 0 0 0 0 0 0 # FracFemale_GP_1
#  M2 parameter for each predator fleet
#
# timevary MG parameters 
#_ LO HI INIT PRIOR PR_SD PR_type  PHASE
 0.0001 2 1.5 0.5 0.5 0 -5 # RecrDist_GP_1_area_1_month_1_dev_se
 -0.99 0.99 0 0 0.5 0 -6 # RecrDist_GP_1_area_1_month_1_dev_autocorr
# info on dev vectors created for MGparms are reported with other devs after tag parameter section 
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; Options: 1=NA; 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn #  parm_name
           3.5            20       		14         11.75             0             0          1          0          0          0          0          0          0          0 # SR_LN(R0)
           0.2         0.999          0.88           0.9             0             0         -4          0          0          0          0          0          0          0 # SR_BH_steep
           0.1             2           0.6           0.4             0             0         -1          0          0          0          0          0          0          0 # SR_sigmaR
            -5             6             0             0             0             0         -1          0          0          0          0          0          0          0 # SR_regime
             0             0             0             0             0             0         -1          0          0          0          0          0          0          0 # SR_autocorr
# timevary SR parameters
# -5 6 0.0419613 0 0 0 3 # SR_regime_BLK2repl_2004
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1983 # first year of main recr_devs; early devs can preceed this era
2023 # last year of main recr_devs; forecast devs start in following year
3 #_recdev phase 
1 # (0/1) to read 13 advanced options
 -10 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 4 #_recdev_early_phase
 -99 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1965 #_last_yr_nobias_adj_in_MPD; begin of ramp
 1989 #_first_yr_fullbias_adj_in_MPD; begin of plateau
 2019 #_last_yr_fullbias_adj_in_MPD
 2024 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS3 sets bias_adj to 0.0 for fcast yrs)
 0.972 #_max_bias_adj_in_MPD (typical ~0.8; -3 sets all years to 0.0; -2 sets all non-forecast yrs w/ estimated recdevs to 1.0; -1 sets biasadj=1.0 for all yrs w/ recdevs)
 0 #_period of cycles in recruitment (N parms read below)
 -5 #min rec_dev
 5 #max rec_dev
 0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
# all recruitment deviations
#  1973E 1974E 1975E 1976E 1977E 1978E 1979E 1980E 1981E 1982E 1983R 1984R 1985R 1986R 1987R 1988R 1989R 1990R 1991R 1992R 1993R 1994R 1995R 1996R 1997R 1998R 1999R 2000R 2001R 2002R 2003R 2004R 2005R 2006R 2007R 2008R 2009R 2010R 2011R 2012R 2013R 2014R 2015R 2016R 2017R 2018R 2019R 2020R 2021F 2022F 2023F
#  -0.00528082 0.0070678 0.0104296 0.0458488 0.0844483 0.113904 0.110933 -0.0217014 0.179204 -0.00722445 -0.135493 0.0611366 0.182234 0.429866 -0.318378 0.460211 0.0881247 0.233169 -0.0109562 -0.587502 -0.375425 0.346508 -0.945245 0.150641 0.154065 0.0820621 0.0462207 -0.118004 -0.366446 -0.0650642 -0.109136 0.454289 0.72512 0.342662 0.530547 0.121746 0.166914 -0.330609 0.0799752 -0.279088 -0.159681 -0.271126 -0.0839391 -0.267002 0.0781038 -0.0523761 0.164211 -0.422337 0 0 0
#
#Fishing Mortality info 
0.3 # F ballpark value in units of annual_F
-2001 # F ballpark year (neg value to disable)
4 # F_Method:  1=Pope midseason rate; 2=F as parameter; 3=F as hybrid; 4=fleet-specific parm/hybrid (#4 is superset of #2 and #3 and is recommended)
10 # max F (methods 2-4) or harvest fraction (method 1)
# read list of fleets that do F as parameter; unlisted fleets stay hybrid, bycatch fleets must be included with start_PH=1, high F fleets should switch early
# (A) fleet, (B) F_starting_value (used if start_PH=1), (C) start_PH for parms (99 to stay in hybrid, <0 to stay at starting value)
# (A) (B) (C)  (terminate list with -9999 for fleet)
 1 0.05 1 # trawlers
 2 0.05 1 # volpal
 3 0.05 99 # artisanal
 4 0.05 1 # cdTrw
-9999 1 1 # end of list
4 #_number of loops for hybrid tuning; 4 good; 3 faster; 2 enough if switching to parms is enabled
#
#_initial_F_parms; for each fleet x season that has init_catch; nest season in fleet; count = 8
#_for unconstrained init_F, use an arbitrary initial catch and set lambda=0 for its logL
#_ LO HI      INIT PRIOR PR_SD  PR_type  PHASE
 1e-05 4 0.0621701   0.1   0.1        0      1 # InitF_seas_1_flt_1trawlers
 1e-05 4 0.104834    0.1   0.1        0      1 # InitF_seas_1_flt_2volpal
 1e-05 4 0.0737925   0.1   0.1        0      1 # InitF_seas_2_flt_1trawlers
 1e-05 4 0.110133    0.1   0.1        0      1 # InitF_seas_2_flt_2volpal
 1e-05 4 0.074167    0.1   0.1        0      1 # InitF_seas_3_flt_1trawlers
 1e-05 4 0.0694638   0.1   0.1        0      1 # InitF_seas_3_flt_2volpal
 1e-05 4 0.0558943   0.1   0.1        0      1 # InitF_seas_4_flt_1trawlers
 1e-05 4 0.0470007   0.1   0.1        0      1 # InitF_seas_4_flt_2volpal
#
# F rates by fleet x season
# Yr:  1960 1960 1960 1960 1961 1961 1961 1961 1962 1962 1962 1962 1963 1963 1963 1963 1964 1964 1964 1964 1965 1965 1965 1965 1966 1966 1966 1966 1967 1967 1967 1967 1968 1968 1968 1968 1969 1969 1969 1969 1970 1970 1970 1970 1971 1971 1971 1971 1972 1972 1972 1972 1973 1973 1973 1973 1974 1974 1974 1974 1975 1975 1975 1975 1976 1976 1976 1976 1977 1977 1977 1977 1978 1978 1978 1978 1979 1979 1979 1979 1980 1980 1980 1980 1981 1981 1981 1981 1982 1982 1982 1982 1983 1983 1983 1983 1984 1984 1984 1984 1985 1985 1985 1985 1986 1986 1986 1986 1987 1987 1987 1987 1988 1988 1988 1988 1989 1989 1989 1989 1990 1990 1990 1990 1991 1991 1991 1991 1992 1992 1992 1992 1993 1993 1993 1993 1994 1994 1994 1994 1995 1995 1995 1995 1996 1996 1996 1996 1997 1997 1997 1997 1998 1998 1998 1998 1999 1999 1999 1999 2000 2000 2000 2000 2001 2001 2001 2001 2002 2002 2002 2002 2003 2003 2003 2003 2004 2004 2004 2004 2005 2005 2005 2005 2006 2006 2006 2006 2007 2007 2007 2007 2008 2008 2008 2008 2009 2009 2009 2009 2010 2010 2010 2010 2011 2011 2011 2011 2012 2012 2012 2012 2013 2013 2013 2013 2014 2014 2014 2014 2015 2015 2015 2015 2016 2016 2016 2016 2017 2017 2017 2017 2018 2018 2018 2018 2019 2019 2019 2019 2020 2020 2020 2020 2021 2021 2021 2021 2022 2022 2022 2022 2023 2023 2023 2023
# seas:  1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4
# trawlers 0.105164 0.129437 0.130391 0.0995844 0.121655 0.150086 0.150686 0.116002 0.128057 0.157688 0.158425 0.121298 0.146296 0.180798 0.182711 0.138695 0.199428 0.250391 0.255813 0.195993 0.255933 0.323034 0.335197 0.258857 0.233129 0.29198 0.296256 0.224514 0.222332 0.273642 0.275378 0.206777 0.215542 0.264076 0.264267 0.19759 0.198801 0.242646 0.241816 0.18073 0.124725 0.148424 0.144344 0.106117 0.109691 0.13063 0.12722 0.0942724 0.241572 0.301918 0.310101 0.239839 0.361466 0.472097 0.508784 0.401384 0.298817 0.373964 0.384343 0.2903 0.423696 0.547442 0.587175 0.452042 0.452563 0.586252 0.6307 0.476928 0.304536 0.366384 0.36141 0.25673 0.234251 0.278033 0.270668 0.192656 0.261149 0.318119 0.316397 0.231281 0.288658 0.35862 0.366415 0.276691 0.25587 0.313799 0.312304 0.227573 0.18029 0.222392 0.226017 0.169992 0.239426 0.308718 0.330514 0.25782 0.253619 0.327801 0.350796 0.263874 0.301688 0.382273 0.402048 0.295373 0.286598 0.351981 0.358528 0.252544 0.190384 0.237602 0.246593 0.190709 0.254395 0.320232 0.337514 0.240318 0.25315 0.310937 0.317659 0.231836 0.228796 0.279795 0.285762 0.206171 0.20287 0.249317 0.255307 0.189899 0.216553 0.279806 0.300156 0.23903 0.189955 0.243738 0.257654 0.197294 0.3088 0.277869 0.278164 0.293506 0.321973 0.548149 0.556975 0.46573 0.404453 0.555703 0.489307 0.338054 0.433937 0.478994 0.547967 0.34231 0.417123 0.428358 0.415031 0.347122 0.277808 0.285702 0.461881 0.306628 0.348476 0.428605 0.409977 0.281036 0.243559 0.396407 0.361743 0.26119 0.207408 0.420225 0.353956 0.186315 0.157521 0.399078 0.388882 0.221925 0.187844 0.334473 0.348686 0.21338 0.234619 0.363003 0.293233 0.182396 0.167919 0.294148 0.353384 0.24088 0.178048 0.366338 0.34589 0.223836 0.215998 0.260321 0.347052 0.20637 0.271782 0.36736 0.35087 0.218294 0.24087 0.439967 0.272673 0.131862 0.308193 0.350378 0.380361 0.246676 0.260375 0.326199 0.261103 0.424978 0.185214 0.172692 0.18611 0.297255 0.22883 0.228746 0.197036 0.1561 0.133029 0.222754 0.316456 0.233368 0.186477 0.312837 0.206985 0.246518 0.137204 0.162308 0.147186 0.135221 0.112526 0.169163 0.200803 0.193994 0.148009 0.175031 0.253611 0.191021 0.12508 0.10753 0.0950657 0.0984805 0.178012 0.20853 0.253658 0.223197 0.178012 0.20853 0.253658 0.223197 0.178012 0.20853 0.253658 0.223197
# volpal 0.17921 0.191637 0.12328 0.0837297 0.208392 0.224195 0.144578 0.0988309 0.223019 0.242257 0.155079 0.106242 0.26391 0.286286 0.185809 0.127217 0.374936 0.41491 0.273889 0.188424 0.506963 0.569178 0.382634 0.267677 0.504531 0.559877 0.369945 0.254195 0.526288 0.576512 0.37555 0.254813 0.547478 0.59267 0.380742 0.254858 0.525095 0.559438 0.354739 0.236192 0.329245 0.337976 0.207192 0.13415 0.277079 0.282978 0.173391 0.11151 0.573773 0.626781 0.409005 0.278202 0.857834 0.997298 0.696642 0.500216 0.780343 0.869936 0.58506 0.406329 1.25301 1.50846 1.08241 0.779831 1.68686 2.06075 1.49034 1.08288 1.45675 1.59755 1.01393 0.652986 1.19839 1.2223 0.732487 0.449039 1.19485 1.24852 0.763051 0.476184 1.17878 1.25555 0.78785 0.515665 0.942531 0.991855 0.623066 0.403198 0.603785 0.639183 0.408352 0.270745 0.974945 1.10066 0.749186 0.521115 1.36653 1.64964 1.16904 0.837268 1.25076 1.43856 1.01323 0.708327 1.45895 1.65437 1.10483 0.737162 1.8755 2.07004 1.34066 0.848934 1.3652 1.44163 0.913262 0.595666 0.990878 1.05609 0.687216 0.451414 0.985942 1.02775 0.633323 0.397202 0.772564 0.793068 0.491621 0.316373 0.668774 0.716554 0.45894 0.306878 0.76882 0.839312 0.555485 0.380192 0.87953 0.757695 0.516847 0.353657 0.615126 0.813392 0.553759 0.358831 0.558823 0.819011 0.446777 0.346597 0.657843 0.671696 0.458282 0.25889 0.465197 0.535926 0.435217 0.328522 0.265732 0.153127 0.146821 0.0888954 0.357192 0.273422 0.189634 0.0824555 0.150431 0.210808 0.11339 0.0543728 0.21277 0.133406 0.119901 0.0660011 0.0850096 0.169079 0.123151 0.0584045 0.123467 0.109881 0.0932807 0.0500356 0.164192 0.106172 0.0894757 0.074738 0.115817 0.16185 0.124115 0.0629021 0.190043 0.211638 0.219454 0.166655 0.232013 0.230858 0.305245 0.166592 0.33183 0.242119 0.277958 0.164425 0.278975 0.286097 0.0865264 0.0389789 0.249595 0.247931 0.244643 0.187913 0.284624 0.170387 0.0780653 0.11681 0.143636 0.159838 0.185332 0.143121 0.133263 0.236583 0.23794 0.242551 0.15761 0.224754 0.183525 0.196426 0.244714 0.230014 0.177505 0.163624 0.341455 0.193119 0.104841 0.180849 0.144212 0.223881 0.153032 0.161017 0.217845 0.239238 0.108169 0.109655 0.130331 0.146706 0.0900913 0.101122 0.227303 0.281515 0.162168 0.171632 0.227303 0.281515 0.162168 0.171632 0.227303 0.281515 0.162168 0.171632
# artisanal 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.144566 0.265146 0.223871 0.153861 0.154775 0.297468 0.265641 0.191941 0.157891 0.31169 0.282786 0.204776 0.192601 0.367426 0.322079 0.227054 0.186436 0.346343 0.29334 0.199814 0.179282 0.318604 0.262988 0.177929 0.212004 0.41566 0.377675 0.272744 0.164335 0.292115 0.238321 0.159114 0.125549 0.228899 0.193189 0.132492 0.133718 0.242943 0.203802 0.139287 0.173901 0.325811 0.286954 0.206309 0.165664 0.321674 0.289834 0.210148 0.165026 0.250723 0.286264 0.201325 0.16898 0.330008 0.250069 0.164236 0.130242 0.318522 0.247462 0.183814 0.174616 0.297138 0.226475 0.1335 0.14047 0.247241 0.223463 0.176132 0.104237 0.2636 0.318841 0.160153 0.149567 0.249364 0.253139 0.144097 0.119933 0.257679 0.276208 0.142591 0.134191 0.228057 0.24104 0.122459 0.0682976 0.160634 0.156266 0.115423 0.0876656 0.156432 0.172847 0.0943414 0.0726999 0.121741 0.133055 0.0628478 0.048034 0.0649132 0.0829929 0.0398057 0.0545989 0.0658453 0.0730232 0.0573135 0.0767719 0.0704988 0.0528094 0.0486138 0.0863766 0.0767365 0.0660419 0.0598393 0.0540165 0.0906869 0.0670102 0.049093 0.0809796 0.0803169 0.0824073 0.0659375 0.0961251 0.0729691 0.101994 0.0757335 0.0653026 0.085925 0.126694 0.0994083 0.0987259 0.10169 0.11431 0.066893 0.0630168 0.0870665 0.109918 0.0769116 0.0713715 0.0830116 0.125844 0.0844839 0.073057 0.06453 0.059469 0.0523112 0.0447544 0.0541551 0.0593958 0.0462201 0.0483679 0.0535025 0.0792896 0.0526956 0.0574895 0.0554518 0.0557422 0.0375891 0.0695274 0.0752967 0.0897542 0.0630151 0.0695274 0.0752967 0.0897542 0.0630151 0.0695274 0.0752967 0.0897542 0.0630151
# cdTrw 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0295788 0.0468949 0.0393973 0.0186095 0.0374268 0.061784 0.0542091 0.0263466 0.0549828 0.0882227 0.0759129 0.0342985 0.0647474 0.102867 0.0885547 0.0396335 0.0808256 0.124287 0.104046 0.0444191 0.0691508 0.114237 0.100155 0.0511861 0.0952335 0.148029 0.124969 0.0515292 0.0736476 0.116818 0.0994166 0.0466386 0.104711 0.1649 0.139162 0.0620984 0.102723 0.163855 0.13931 0.065119 0.0885025 0.150246 0.135227 0.0698156 0.0689295 0.111532 0.0960506 0.0443377 0.0524741 0.0536885 0.0486018 0.0196223 0.0302082 0.0739532 0.0711784 0.0489166 0.14233 0.263589 0.229052 0.0822809 0.126626 0.203137 0.151663 0.0535786 0.0803236 0.0833155 0.0814422 0.04447 0.0631892 0.0884742 0.0989939 0.0568065 0.0757153 0.112426 0.0843497 0.0509621 0.15733 0.278449 0.191515 0.0936641 0.144354 0.180096 0.145252 0.0887355 0.19204 0.241028 0.185062 0.0654329 0.140247 0.173238 0.134648 0.034469 0.0727383 0.0657713 0.0836184 0.0292439 0.0347616 0.048235 0.0412802 0.00822384 0.0234792 0.0354714 0.0347343 0.0107523 0.0234807 0.0229015 0.0458251 0.0129672 0.0169241 0.047838 0.048328 0.0108308 0.0266197 0.061294 0.0628225 0.0217848 0.04152 0.0531824 0.0368462 0.016827 0.0168438 0.0362183 0.0397727 0.0196241 0.0411254 0.0639803 0.0723961 0.0359649 0.0705898 0.0717293 0.0748221 0.0475879 0.0551879 0.0535153 0.0612421 0.0413134 0.0603177 0.0628066 0.0650494 0.0469369 0.0396802 0.0379845 0.0451753 0.0255047 0.0465941 0.0347172 0.0514211 0.0312401 0.0382412 0.0364411 0.0313435 0.0175346 0.0231381 0.032339 0.0283209 0.0180565 0.049844 0.0477778 0.0512807 0.0308515 0.049844 0.0477778 0.0512807 0.0308515 0.049844 0.0477778 0.0512807 0.0308515
#
#_Q_setup for fleets with cpue or survey data
#_1:  fleet number
#_2:  link type: (1=simple q, 1 parm; 2=mirror simple q, 1 mirrored parm; 3=q and power, 2 parm; 4=mirror with offset, 2 parm)
#_3:  extra input for link, i.e. mirror fleet# or dev index number
#_4:  0/1 to select extra sd parameter
#_5:  0/1 for biasadj or not
#_6:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
         5         1         0         1         0         1  #  SpSurv
         6         1         0         1         0         1  #  PtSurv
         7         1         0         1         0         1  #  CdSurv
         8         1         0         0         0         1  #  SpCPUE_trawlers
         9         1         0         0         0         1  #  SpCPUE_volpal
-9999 0 0 0 0 0
#
#_Q_parms(if_any);Qunits_are_ln(q)
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
           -15            -1      -6.54593           -10             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_SpSurv(5)
             0             1           0.2           0.1           0.1             0         -4          0          0          0          0          0          0          0  #  Q_extraSD_SpSurv(5)
           -15            -1      -6.49808           -10             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_PtSurv(6)
             0             1           0.2           0.1           0.1             0         -4          0          0          0          0          0          0          0  #  Q_extraSD_PtSurv(6)
           -15          -0.1      -7.18498           -10             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_CdSurv(7)
             0             1           0.2           0.1           0.1             0         -4          0          0          0          0          0          0          0  #  Q_extraSD_CdSurv(7)
           -15            -1      -10.5205           -10             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_SpCPUE_trawlers(8)
           -16            -1      -9.86751           -10             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_SpCPUE_volpal(9)
#_no timevary Q parameters
#
#_size_selex_patterns
#Pattern:_0;  parm=0; selex=1.0 for all sizes
#Pattern:_1;  parm=2; logistic; with 95% width specification
#Pattern:_2;  parm=6; modification of pattern 24 with improved sex-specific offset
#Pattern:_5;  parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
#Pattern:_11; parm=2; selex=1.0  for specified min-max population length bin range
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_6;  parm=2+special; non-parm len selex
#Pattern:_43; parm=2+special+2;  like 6, with 2 additional param for scaling (average over bin range)
#Pattern:_8;  parm=8; double_logistic with smooth transitions and constant above Linf option
#Pattern:_9;  parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
#Pattern:_21; parm=2+special; non-parm len selex, read as pairs of size, then selex
#Pattern:_22; parm=4; double_normal as in CASAL
#Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
#Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL), using joiners
#Pattern:_25; parm=3; exponential-logistic in length
#Pattern:_27; parm=special+3; cubic spline in length; parm1==1 resets knots; parm1==2 resets all 
#Pattern:_42; parm=special+3+2; cubic spline; like 27, with 2 additional param for scaling (average over bin range)
#_discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
 24 1 0 0 # 1 trawlers
  1 0 0 4 # 2 volpal
 24 0 0 0 # 3 artisanal
 24 0 0 0 # 4 cdTrw
 24 0 0 0 # 5 SpSurv
 24 0 0 0 # 6 PtSurv
 24 0 0 0 # 7 CdSurv
 15 0 0 1 # 8 SpCPUE_trawlers
 15 0 0 2 # 9 SpCPUE_volpal
#
#_age_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
#Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
#Pattern:_11; parm=2; selex=1.0  for specified min-max age
#Pattern:_12; parm=2; age logistic
#Pattern:_13; parm=8; age double logistic
#Pattern:_14; parm=nages+1; age empirical
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_16; parm=2; Coleraine - Gaussian
#Pattern:_17; parm=nages+1; empirical as random walk  N parameters to read can be overridden by setting special to non-zero
#Pattern:_41; parm=2+nages+1; // like 17, with 2 additional param for scaling (average over bin range)
#Pattern:_18; parm=8; double logistic - smooth transition
#Pattern:_19; parm=6; simple 4-parm double logistic with starting age
#Pattern:_20; parm=6; double_normal,using joiners
#Pattern:_26; parm=3; exponential-logistic in age
#Pattern:_27; parm=3+special; cubic spline in age; parm1==1 resets knots; parm1==2 resets all 
#Pattern:_42; parm=2+special+3; // cubic spline; with 2 additional param for scaling (average over bin range)
#Age patterns entered with value >100 create Min_selage from first digit and pattern from remainder
#_Pattern Discard Male Special
 0 0 0 0 # 1 trawlers
 0 0 0 0 # 2 volpal
 0 0 0 0 # 3 artisanal
 0 0 0 0 # 4 cdTrw
 0 0 0 0 # 5 SpSurv
 0 0 0 0 # 6 PtSurv
 0 0 0 0 # 7 CdSurv
 0 0 0 0 # 8 SpCPUE_trawlers
 0 0 0 0 # 9 SpCPUE_volpal
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
# 1   trawlers LenSelex
             6            60       12.4617            15          0.01             0          5          0          0          0          0          0          0          0  #  Size_DblN_peak_trawlers(1)
           -16             5      -1.78604            -2          0.01             0          6          0          0          0          0          0          0          0  #  Size_DblN_top_logit_trawlers(1)
           -16            14       1.63788             4          0.01             0          5          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_trawlers(1)
           -16            40       6.72094            10          0.01             0          5          0          0          0          0          0          0          0  #  Size_DblN_descend_se_trawlers(1)
          -999          -999          -999          -999             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_start_logit_trawlers(1)
          -999          -999          -999          -999             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_end_logit_trawlers(1)
             4            30        27.711            27          0.01             0          5          0          3       1994       1997          5          1          2  #  Retain_L_infl_trawlers(1)
          0.00001         10       1.65486           0.8          0.01             0          5          0          0          0         0            0         1          2  #  Retain_L_width_trawlers(1)
           -10           999           999            10          0.01             0         -6          0          0          0          0          0          0          0  #  Retain_L_asymptote_logit_trawlers(1)
             0             0             0             0          0.01             0         -6          0          0          0          0          0          0          0  #  Retain_L_maleoffset_trawlers(1)
# 2   volpal LenSelex

19	   70	45	50	99	0	5	0	0	0	0	0.5	 0	0		#infl_for_logistic
0.01   60	20	15	99	0	5	0	0	0	0	0.5	 0	0		#95%width_for_logistic

# 3   artisanal LenSelex
             6            70       34.4646            15          0.01             0          5          0          0          0          0          0          0          0  #  Size_DblN_peak_artisanal(3)
           -16             5      -13.1426           -10            10             0          6          0          0          0          0          0          0          0  #  Size_DblN_top_logit_artisanal(3)
           -16            14       4.09915             4          0.01             0          5          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_artisanal(3)
           -16            40       6.61423            10          0.01             0          5          0          0          0          0          0          0          0  #  Size_DblN_descend_se_artisanal(3)
          -999          -999          -999          -999             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_start_logit_artisanal(3)
          -999          -999          -999          -999             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_end_logit_artisanal(3)
# 4   cdTrw LenSelex
             6            60       9.64525            15          0.01             0          5          0          0          0          0          0          0          0  #  Size_DblN_peak_cdTrw(4)
           -16             5      -3.62299            -10         10               0          6          0          0          0          0          0          0          0  #  Size_DblN_top_logit_cdTrw(4)
           -16            14      -6.52996             4          0.01             0          5          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_cdTrw(4)
           -16            40       6.25427            10          0.01             0          5          0          0          0          0          0          0          0  #  Size_DblN_descend_se_cdTrw(4)
          -999          -999          -999          -999             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_start_logit_cdTrw(4)
          -999          -999          -999          -999             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_end_logit_cdTrw(4)
# 5   SpSurv LenSelex
             6            60       6.00004            10            40             0          5          0          0          0          0          0          0          0  #  Size_DblN_peak_SpSurv(5)
           -17             5      -14.8573           -10            10             0          6          0          0          0          0          0          0          0  #  Size_DblN_top_logit_SpSurv(5)
           -16            14       11.6792            10           100             0          5          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_SpSurv(5)
           -16            40        5.6154            10          0.01             0          5          0          0          0          0          0          0          0  #  Size_DblN_descend_se_SpSurv(5)
          -999          -999          -999          -999             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_start_logit_SpSurv(5)
          -999          -999          -999          -999             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_end_logit_SpSurv(5)
# 6   PtSurv LenSelex
             6            65       18.0326            15          0.01             0          5          0          0          0          0          0          0          0  #  Size_DblN_peak_PtSurv(6)
           -16             5      -14.4234            -2          0.01             0          6          0          0          0          0          0          0          0  #  Size_DblN_top_logit_PtSurv(6)
           -11            14        3.6772             4          0.01             0          5          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_PtSurv(6)
           -16            40       5.88544            10          0.01             0          5          0          0          0          0          0          0          0  #  Size_DblN_descend_se_PtSurv(6)
          -999          -999          -999          -999             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_start_logit_PtSurv(6)
          -999          -999          -999          -999             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_start_logit_PtSurv(6)
# 7   CdSurv LenSelex
             6            60       6.00004            10            40             0          6          0          0          0          0          0          0          0  #  Size_DblN_peak_CdSurv(7)
           -17             5           -15            -2          0.01             0         -6          0          0          0          0          0          0          0  #  Size_DblN_top_logit_CdSurv(7)
           -15            14       7.90058            10          100              0          6          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_CdSurv(7)
           -16            40       6.26776            10          0.01             0          6          0          0          0          0          0          0          0  #  Size_DblN_descend_se_CdSurv(7)
          -999          -999          -999          -999             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_start_logit_CdSurv(7)
          -999          -999          -999          -999             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_end_logit_CdSurv(7)
# 8   SpCPUE_trawlers LenSelex
# 9   SpCPUE_volpal LenSelex
# 10   PtCPUE LenSelex
# 1   trawlers AgeSelex
# 2   volpal AgeSelex
# 3   artisanal AgeSelex
# 4   cdTrw AgeSelex
# 5   SpSurv AgeSelex
# 6   PtSurv AgeSelex
# 7   CdSurv AgeSelex
# 8   SpCPUE_trawlers AgeSelex
# 9   SpCPUE_volpal AgeSelex
# 10   PtCPUE AgeSelex
#_No_Dirichlet parameters
# timevary selex parameters 
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type    PHASE  #  parm_name
             4            30             4            27          0.01             0      -2  # Retain_L_infl_trawlers(1)_BLK1repl_1960
             0.5          30             1            27          0.01             0      -2  # Retain_L_widt_trawlers(1)_BLK1repl_1960
        0.0001             2           0.5           0.5           0.5             0      -5  # Retain_L_infl_trawlers(1)_dev_se
         -0.99          0.99             0             0           0.5             0      -6  # Retain_L_infl_trawlers(1)_dev_autocorr
# info on dev vectors created for selex parms are reported with other devs after tag parameter section 
#
0   #  use 2D_AR1 selectivity(0/1)
#_no 2D_AR1 selex offset used
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read and autogen if tag data exist; 1=read
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# deviation vectors for timevary parameters
#  base   base first block   block  env  env   dev   dev   dev   dev   dev
#  type  index  parm trend pattern link  var  vectr link _mnyr  mxyr phase  dev_vector
#      2     4     1     2     2     0     0     0     0     0     0     0
#      5     7     2     1     2     0     0     1     3  1994  1997     5 -4.7594 0.161304 0.657099 0.419124
#      5    14     5    -1     0     0     0     0     0     0     0     0
#
# Input variance adjustments factors: 
 #_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
#_Factor  Fleet  Value
 -9999   1    0  # terminator
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 0 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method

4	1	1	0.3	 1	#_length
4	2	1	0.3	 1	#_length
4	3	1	0.3	 1	#_length
4	4	1	0.3	 1	#_length
4	5	1	0.3	 1	#_length
4	6	1	0.3	 1	#_length
4	7	1	0.3	 1	#_length
4	8	1	0.3	 1	#_length
4	9	1	0.3	 1	#_length

6	1	1	0.3	 1	#_sizefreq
6	2	1	0.3	 1	#_sizefreq
-9999  1  1   1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  0 #_CPUE/survey:_1
#  0 #_CPUE/survey:_2
#  0 #_CPUE/survey:_3
#  0 #_CPUE/survey:_4
#  1 #_CPUE/survey:_5
#  1 #_CPUE/survey:_6
#  1 #_CPUE/survey:_7
#  1 #_CPUE/survey:_8
#  1 #_CPUE/survey:_9
#  1 #_CPUE/survey:_10
#  1 #_discard:_1
#  0 #_discard:_2
#  0 #_discard:_3
#  0 #_discard:_4
#  0 #_discard:_5
#  0 #_discard:_6
#  0 #_discard:_7
#  0 #_discard:_8
#  0 #_discard:_9
#  0 #_discard:_10
#  1 #_lencomp:_1
#  1 #_lencomp:_2
#  1 #_lencomp:_3
#  1 #_lencomp:_4
#  1 #_lencomp:_5
#  1 #_lencomp:_6
#  1 #_lencomp:_7
#  0 #_lencomp:_8
#  0 #_lencomp:_9
#  0 #_lencomp:_10
#  1 #_sizefreq:_1
#  1 #_sizefreq:_2
#  1 #_init_equ_catch1
#  1 #_init_equ_catch2
#  1 #_init_equ_catch3
#  1 #_init_equ_catch4
#  1 #_init_equ_catch5
#  1 #_init_equ_catch6
#  1 #_init_equ_catch7
#  1 #_init_equ_catch8
#  1 #_init_equ_catch9
#  1 #_init_equ_catch10
#  1 #_recruitments
#  1 #_parameter-priors
#  1 #_parameter-dev-vectors
#  1 #_crashPenLambda
#  0 # F_ballpark_lambda
0 # (0/1/2) read specs for more stddev reporting: 0 = skip, 1 = read specs for reporting stdev for selectivity, size, and numbers, 2 = add options for M,Dyn. Bzero, SmryBio
 # 0 2 0 0 # Selectivity: (1) fleet, (2) 1=len/2=age/3=both, (3) year, (4) N selex bins
 # 0 0 # Growth: (1) growth pattern, (2) growth ages
 # 0 0 0 # Numbers-at-age: (1) area(-1 for all), (2) year, (3) N ages
 # -1 # list of bin #'s for selex std (-1 in first bin to self-generate)
 # -1 # list of ages for growth std (-1 in first bin to self-generate)
 # -1 # list of ages for NatAge std (-1 in first bin to self-generate)
999

