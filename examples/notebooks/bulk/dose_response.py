import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import os
from matplotlib import gridspec


def get_pheno(i):
    if i==1: return 3,"HY-2"
    if i==2: return 2,"Trudy"
    if i==3: return 4,"NE"
    if i==4: return 1,"ML"
    raise(ValueError("Bad phenotype value"))

directory = "/data/dwooten/School/qlab/Projects/SCLC/Teicher"
dose_response = pd.read_csv(os.path.join(directory,"sclc_export_experiment_dose_resp_element_19FEB2016.csv"))
sensitivites_parent = pd.read_csv(os.path.join(directory,"sensitivity_all.csv"),index_col="experiment_id")





#for drug in ['Topotecan','Etoposide','Irinotecan','Paclitaxel']:
#for drug in ['Topotecan','AZD-1152','Doxorubicin','VS-507','AZD-0530','AZD-8330','PF-04929113','Etoposide','Irinotecan','Paclitaxel','Tamoxifen']:
#for drug in ['MLN-8237 ']:
print len(sensitivites_parent['drug'].unique())
dlist = list(sensitivites_parent['drug'].unique())
dlist.sort()
for drug in dlist:


    sensitivites = sensitivites_parent.loc[sensitivites_parent['drug']==drug]
    clusters = dict()
    for i in range(1,5): clusters[i] = sensitivites.loc[sensitivites['cluster_2017']==i]

    if False:
        fig = plt.figure()
        for i in range(1,5):
            ax = fig.add_subplot(2,2,i)
            for experiment in clusters[i].index:
                dr = dose_response.loc[dose_response['experiment_id']==experiment]
                ax.plot(np.log10(dr['concentration']),dr['mean_pct_ctrl'],c='k',alpha=0.1)
            ax.set_xlim(np.log10(1.52420000e-09),-5)
            ax.set_ylim(0,120)
            ax.grid()
        plt.show()


    if True:
        fig = plt.figure()
        
        ax = plt.subplot(111)
        ax.set_xlabel("\nlog10 drug concentration")
        ax.set_ylabel("Percent relative to CTRL\n")
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        for j in range(1,5):
            i,cell_type = get_pheno(j)

            
            ax = fig.add_subplot(2,2,i)
            dr = dose_response.loc[dose_response['experiment_id'].isin(clusters[j].index)]
            print i, dr.shape, dr.isnull().sum()['mean_pct_ctrl'],drug
            dr_curve_x = []
            dr_curve_y = []
            dr_curve_lower = []
            dr_curve_upper = []
            for con in dr['concentration'].unique():
                dr_con = dr.loc[dr['concentration']==con]
    #            print con, dr_con.isnull().sum()['mean_pct_ctrl']
                dr_con = dr_con.dropna()
#                if dr_con.shape[0] < 6: continue
                dr_curve_x.append(con)
                dr_curve_y.append(dr_con.median()['mean_pct_ctrl'])
                dr_curve_lower.append(dr_con.quantile(q=0.25)['mean_pct_ctrl'])
                dr_curve_upper.append(dr_con.quantile(q=0.75)['mean_pct_ctrl'])
            ax.fill_between(np.log10(dr_curve_x),dr_curve_lower,dr_curve_upper,color='r',alpha=0.3)
            ax.plot(np.log10(dr_curve_x),dr_curve_y,c='k')
            ax.set_xlim(np.log10(1.52420000e-09),-5)
            ax.plot([np.log10(1.52420000e-09),-5],[50,50])
            ax.set_ylim(0,110)
            ax.grid()
            if i == 1 or i == 2: ax.set_xticklabels([])
            if i == 2 or i == 4: ax.set_yticklabels([])
            ax.set_title("%s"%(cell_type))
    #    axy = plt.subplot2grid((gs,gs),(0,0),rowspan=gs)
    #    axx = plt.subplot2grid((gs,gs),(gs-1,0),colspan=gs)

#        plt.text(-10.8,-20,"log10 drug concentration")
#        plt.text(-14.3,170,"Percent relative to CTRL",rotation='vertical')
        plt.suptitle("%s dose response"%drug)
        plt.savefig("all_data/All_Drugs/%s.png"%repr(drug).replace("/","|"))
        plt.clf()
        plt.cla()
        plt.close()
