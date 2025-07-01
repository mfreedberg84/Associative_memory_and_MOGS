# Examining the effect of task engagement during inter-trial rest periods on “Micro-offline gains”

Nafiz Ahmed, Tharan Suresh, Sara Jeanne Hussain, Michael Freedberg

Explicit motor sequence learning (ESL) consolidation manifests as skill improvements during rest periods that last minutes or longer. More recent evidence, however, suggest that ESL improvements also occur during inter-trial rest periods that last seconds (micro-offline gains, or MOGS). While indirect evidence supports that MOGS is a phenomenon tied to brief periods of wakeful rest, no study to our knowledge has tested whether wakeful rest is necessary for MOGS. Motor and procedural learning, in general, are also thought to rely on medial temporal lobe activity, possibly through associative memory processes supported by the hippocampus, suggesting that MOGS also relies on associative processes. However, this hypothesis remains to be tested. In this study, we reasoned that if wakeful rest is necessary for MOGS, replacing inter-trial rest periods with an engaging task should degrade MOGS. Further, if associative processes support MOGS during these inter-trial periods, encoding associative memories during them should further degrade MOGS. To test these hypotheses, we compared the performance of three groups that differed only in whether the 10 second inter-trial periods included no task (REST), an associative encoding task (ENC), or a semantic judgment task (SEM). Our results revealed no significant group differences in motor skill acquisition, retention, or MOGS. These results are consistent with the theory that MOGS capture inter-trial performance changes in ESL  rather than a rapid form of consolidation supported by associative processes. 

Instructions: Download all files and set the relevant paths at the top of the **a_MOGS_calc.m**, **b_MOGS_group_analysis.m**, and **c_acquisition_retention.m** scripts. Open matlab and navigate to the folder with all the scripts. The scripts are meant to be run in this order in MATLAB:

1) **a_MOGS_calc.m**. To initialize this script, you need to set the path to your scripts (scripts_dir) and output (out_dir) folders. This script will calculate
    a) the number of correct sequences on each trial,
    b) the keypresses per second (KPS) of all correct sequences,
    c) associative memory encoding success (ENC group only),
    d) microoffline gains, and
    e) microonline gains for all participants.

   The .mat files with this information get saved to the out_dir folder. 

3) **b_MOGS_group_analysis.m**. To initiailze this script, you need to set the path to where "a_MOGS_calc.m" outputs its files. You also have the option to define the early learning break point at the top of the script. 
    
4) **c1_acquisition_retention.m**.  To initiailze this script, you need to set the path to where "a_MOGS_calc.m" outputted its files. 

5) **d_associative_memory_and_acq.m**  To initiailze this script, you need to set the path to where "a_MOGS_calc.m" outputs its files. You also have the option to define the early learning break point at the top of the script.
   
6) **e_associative_memory_and_mogs**.m  To initiailze this script, you need to set the path to where "a_MOGS_calc.m" outputs its files. You also have the option to define the early learning break point at the top of the script. 

In addition, there are several downloadable functions that need to go in the folder with your scripts:
a) Violinplot: https://www.mathworks.com/matlabcentral/fileexchange/170126-violinplot-matlab
b) shadedErrorBar: https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar

The scripts will generate several pop-up figures that are also shown in the paper. The result of statistical tests and descriptive statistics are outputted on the command line.

For the Bayes analyses, results are already stored in the "Group_data_ANOVA.jasp" file, which reads "Group_data_ANOVA.csv." This file contains the results contrasting groups on acquisition asymptote, acquisition magnitude, acquisition rate, retention, and MOGS. 
