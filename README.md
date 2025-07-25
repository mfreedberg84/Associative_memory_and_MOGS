# Investigating the effect of task engagement during inter-trial rest periods on micro-offline gains

Nafiz Ishtiaque Ahmed, Tharan Suresh, Sara Jeanne Hussain, Michael Freedberg

Explicit motor sequence learning (ESL) consolidation manifests as skill improvements during rest periods that last minutes or longer. However, recent evidence suggests that ESL improvements also occur during inter-trial rest periods that last seconds (micro-offline gains;  MOGS). While indirect evidence supports that MOGS is a phenomenon tied to brief periods of wakeful rest, this hypothesis has never been directly tested. Prior studies suggest that MOGS and sequence learning in general relies on associative memory processes that link sequence elements across time and space. However, evidence supporting this hypothesis in healthy humans is lacking. We reasoned that if wakeful rest during inter-trial ESL periods is necessary for MOGS, replacing these periods with an engaging task should degrade MOGS. Additionally, we explored whether associative processes support MOGS, predicting that associative memory encoding during inter-trial periods would further degrade MOGS. To this end, we compared the performance of three groups that differed only in whether inter-trial periods included no task (REST), an associative encoding task (ENC), or a semantic judgment task (SEM). We identified no significant group differences in ESL or MOGS, which were confirmed by Bayesian analyses. These results are consistent with the notion that MOGS capture inter-trial performance changes rather than a rapid form of consolidation that requires rest or depends on associative memory processes.

Instructions: Download all files and set the relevant paths at the top of all scripts. Open matlab and navigate to the folder with all the scripts. The scripts are organized to display the data as shown in the paper. The scripts are meant to be run in this order in MATLAB:

1) **a_MOGS_calc.m**.This script will calculate
    a) the number of correct sequences on each trial,
    b) the keypresses per second (KPS) of all correct sequences,
    c) associative memory encoding success (ENC group only),
    d) microoffline gains, and
    e) microonline gains for all participants.

   The .mat files with these data gets saved to the out_dir folder. 

3) **b_reproducing_bonstrup.m** 
4) **d_break_point_analysis_MOGS.m**
5) **e_group_differences_MOnGS.m**    
6) **f_MOGS_and_MOnGS_corr.m**
7) **g_break_point_analysis_MOnGS.m**
8) **h1_acquisition_retention.m**
    a) This script has four subscripts (only h1 needs to be run on the command line).
       i) h2_curvefit_mogs.m
       ii) h3_exponentialfit.m
       iii) h4_exponential_ls.m
       iv) h5_exponential.m
9) **i_post_hoc_sample_size_analysis.m**
10) **j_associative_memory_and_acq.m**
11) **k_associative_memory_and_mogs.m**
    
Additionally, there are several downloadable functions that need to go in the folder with your scripts:
a) Violinplot: https://www.mathworks.com/matlabcentral/fileexchange/170126-violinplot-matlab
b) shadedErrorBar: https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar

The scripts will generate several pop-up figures that are also shown in the paper. The result of statistical tests and descriptive statistics are outputted on the command line.

For the Bayes analyses, results are already stored in the "Group_data_ANOVA.jasp" file, which reads "Group_data_ANOVA.csv." This file contains the results contrasting groups on acquisition asymptote, acquisition magnitude, acquisition rate, retention, and MOGS. 
