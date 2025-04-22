# Do “Micro-offline gains” rely on rest or hippocampal-mediated processes?

Nafiz Ahmed, Tharan Suresh, Sara Hussain, Michael Freedberg

Motor learning consolidation manifests as skill improvements during rest periods that last minutes or longer. More recent evidence, however, suggest that skill improvements also occur during inter-trial rest periods that last seconds (micro-offline gains, or MOGS). While indirect evidence supports that MOGS is a phenomenon tied to brief periods of wakeful rest, no study to our knowledge has tested whether wakeful rest is necessary for MOGS. MOGS is also thought to rely on activity in the hippocampus, a region classically associated with episodic rather than motor memory. However, experimental evidence showing that hippocampal activity is necessary for MOGS in healthy humans is also lacking. In this study, we reasoned that if wakeful rest is necessary for MOGS, replacing inter-trial rest periods with an engaging task should degrade MOGS. Further, if hippocampal-mediated processes support MOGS during these inter-trial periods, performing a hippocampal-dependent task like episodic memory encoding during them should further degrade MOGS. To test these hypotheses, we recruited forty-five neurotypical adult participants to perform an explicit sequence learning (ESL) task. Participants were evenly divided into three groups that differed only in whether the 10 second inter-trial periods included no task (REST), a word pair encoding task (ENC), or a semantic judgment task (SEM). Contrary to our hypotheses, our results revealed no significant differences in acquisition rate, acquisition magnitude, retention, or MOGS between groups. Because ENC or SEM task performance did not affect either MOGS or conventional motor learning metrics, our results are consistent with the notion that MOGS capture inter-trial performance changes rather than a hippocampal-dependent rapid form of consolidation. 

Instructions: Download all files and set the relevant paths at the top of the **a_MOGS_calc.m**, **b_MOGS_group_analysis.m**, and **c_acquisition_retention.m** scripts. Open matlab and navigate to the folder with all the scripts. The scripts are meant to be run in this order in MATLAB:

1) **a_MOGS_calc.m**. To initialize this script, you need to set the path to your scripts (scripts_dir) and output (out_dir) folders. This script will calculate
    a) the number of correct sequences on each trial,
    b) the keypresses per second (KPS) of all correct sequences,
    c) episodic memory encoding success (ENC group only),
    d) microoffline gains, and
    e) microonline gains for all participants.

   The .mat files with this information get saved to the out_dir folder. 

3) **b_MOGS_group_analysis.m**. To initiailze this script, you need to set the path to where "a_MOGS_calc.m" outputs its files. You also have the option to define the early learning break point at the top of the script. 
    
4) **c1_acquisition_retention.m**.  To initiailze this script, you need to set the path to where "a_MOGS_calc.m" outputted its files. 

In addition, there are several downloadable functions that need to go in the folder with your scripts:
1) Violinplot: https://www.mathworks.com/matlabcentral/fileexchange/170126-violinplot-matlab
2) shadedErrorBar: https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar

The scripts will generate several pop-up figures that are also shown in the paper. The result of statistical tests and descriptive statistics are outputted on the command line.
