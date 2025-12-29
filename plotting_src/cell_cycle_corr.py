from ggp_notebooks.plotting_src.header_settings import *

from ggp_notebooks.plotting_src.read_ggp_run import * # you might need to change these two line, such that you can use read ggp_run
from ggp_notebooks.plotting_src.cells_utils import * # you might need to change these two line, such that you can use read ggp_run
from ggp_notebooks.plotting_src.read_utils import * # you might need to change these two line, such that you can use read ggp_run



def split_array(arr, sorting_criterium, nbins):
    bin_edges = [np.min(sorting_criterium)]
    splitted_arr = []
    for i in range(nbins):
        bin_edges.append(np.percentile(sorting_criterium, (i+1)/nbins * 100))
        
        if i < nbins-1:
            arr_in_bin = arr[np.where( np.logical_and(sorting_criterium >= bin_edges[-2], sorting_criterium < bin_edges[-1]) )]

        else:
            arr_in_bin = arr[np.where(sorting_criterium >= bin_edges[-2] )]
        splitted_arr.append(arr_in_bin)
    return splitted_arr, bin_edges



def cell_cycle_dependent_cross_corr(filename,  
                                    xlim=[None, None], 
                                    ylim=[None, None],
                                    outfile=None, 
                                    show=True, 
                                    xaxis_parameter="volume",
                                    alignment="norm", 
                                    label="both",
                                    title=None,
                                    color_code="condition", 
                                    bins=10):    
    
    data_set_name = re.split('_f\d', filename.split('/')[-1] )[0]
    condition = data_set_name.split("_")[0]
    promoter = data_set_name.split("_")[1]

    if label=="both":
        labeli = condition + " " +promoter
    elif label=="condition":
        labeli = condition
    elif label=="promoter":
        labeli = promoter

    if bins>1:
        fig, axes = plt.subplots( 1,bins+1, figsize=(20,3.2))
        colors=["tab:red", "tab:orange" , "tab:olive", "tab:green", "tab:cyan", "tab:blue", "tab:purple"]
        
    else:
        fig, axes = plt.subplots( 1,bins, figsize=(5,3.2))
        if color_code=="promoter":
            colors = [color_by_promoter[promoter]]

        elif color_code=="condition":
            colors = [color_by_condition[condition]]
        elif color_code==None:
            colors = ["tab:blue"]*10
        
    fig.suptitle(labeli)
        
    ax = np.ravel(axes)
    
    for i in range(1,bins):
        ax[i].sharex(ax[0])
        ax[i].sharey(ax[0])

    for a in ax:
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)

        # turn off the right spine/ticks
        a.spines['right'].set_color('none')
        a.yaxis.tick_left()

        # turn off the top spine/ticks
        a.spines['top'].set_color('none')
        a.xaxis.tick_bottom()
        # remove first tick
#         a.xaxis.get_major_ticks()[0].set_visible(False)


    data = pd.read_csv(filename, skiprows=header_lines(filename))
    cells_data = df2ggp_cells(data)

    collect_mean_l = np.array([])
    collect_var_l = np.array([])

    collect_mean_q = np.array([])
    collect_var_q = np.array([])
    
    collect_cov_lq = np.array([])
    

    xaxis = np.array([])

    for i, cell in enumerate(cells_data):
        if len(cell.time)<2:
            continue

        if xaxis_parameter == "volume":
            single_cell_xaxis = np.exp(cell.mean_x)
        elif xaxis_parameter == "time":
            single_cell_xaxis = np.array(cell.time)

        if alignment=="norm":                    
            dx = (single_cell_xaxis -single_cell_xaxis[0])  \
                / (single_cell_xaxis[-1] - single_cell_xaxis[0]) 

        elif alignment== "birth":
            dx = single_cell_xaxis - single_cell_xaxis[0]
        else:
            return

        xaxis = np.append(xaxis, dx)

        collect_mean_l = np.append(collect_mean_l, cell.mean_l)
        collect_var_l = np.append(collect_var_l, cell.cov_ll)

        collect_mean_q = np.append(collect_mean_q, cell.mean_q)
        collect_var_q = np.append(collect_var_q, cell.cov_qq)
        
        collect_cov_lq = np.append(collect_cov_lq, cell.cov_lq)

    xaxis_split, _ = split_array(xaxis, xaxis, bins)

    mean_l_split, bin_edges = split_array(collect_mean_l, xaxis, bins)
    cov_ll_split, _ = split_array(collect_var_l, xaxis, bins)
    
    mean_q_split, _ = split_array(collect_mean_q, xaxis, bins)
    cov_qq_split, _ = split_array(collect_var_q, xaxis, bins)
    
    cov_lq_split, _ = split_array(collect_cov_lq, xaxis, bins)

    bin_centers = []
    corrs = []
    for i, _ in enumerate(xaxis_split):        
        cov_lq = np.mean(np.array(mean_l_split[i]) * np.array(mean_q_split[i]) +  np.array(cov_lq_split[i])) \
                - np.mean(np.array(mean_l_split[i]))*np.mean(np.array(mean_q_split[i]))
        
        cov_ll = np.mean(np.array(mean_l_split[i]) * np.array(mean_l_split[i]) +  np.array(cov_ll_split[i])) \
                - np.mean(np.array(mean_l_split[i]))*np.mean(np.array(mean_l_split[i]))
        
        cov_qq = np.mean(np.array(mean_q_split[i]) * np.array(mean_q_split[i]) +  np.array(cov_qq_split[i])) \
                - np.mean(np.array(mean_q_split[i]))*np.mean(np.array(mean_q_split[i]))
        
#         cov_lq = np.mean(np.array(mean_l_split[i]) * np.array(mean_q_split[i])) \
#                 - np.mean(np.array(mean_l_split[i]))*np.mean(np.array(mean_q_split[i]))
        
#         cov_ll = np.mean(np.array(mean_l_split[i]) * np.array(mean_l_split[i])) \
#                 - np.mean(np.array(mean_l_split[i]))*np.mean(np.array(mean_l_split[i]))
        
#         cov_qq = np.mean(np.array(mean_q_split[i]) * np.array(mean_q_split[i])) \
#                 - np.mean(np.array(mean_q_split[i]))*np.mean(np.array(mean_q_split[i]))
        
        corr = cov_lq/np.sqrt((cov_ll * cov_qq))
                
        label = "Corr={:.2f}".format(corr)
        ax[i].set_title ("Corr={:.2f}".format(corr))
        ax[i].errorbar(mean_l_split[i]*60, mean_q_split[i]*60, 
#                        xerr=np.sqrt(cov_ll_split[i])*60, 
#                        yerr=np.sqrt(cov_qq_split[i])*60, 
                       marker='.', 
                       markersize=3, 
                       alpha=0.2, 
                       color=colors[i], 
                       fmt='o')       
#         ax[i].legend(loc='lower center', bbox_to_anchor=(0.5, 1.02))
        
        bin_center = (np.max(xaxis_split[i])+np.min(xaxis_split[i]))/2
        ax[i].set_xlabel(r"$\lambda$ (1/h)")
        
        if bins>1:
            ax[-1].scatter(bin_center , corr, color=colors[i])
        

    ax[0].set_ylabel("q (1/h)")
    
    if bins>1:
        ax[-1].set_xlabel("norm. added length")
        ax[-1].set_ylabel("Corr")
        
    fig.tight_layout()
    

    if outfile != None:
        fig.savefig(outfile.format(xaxis_parameter), dpi=600, facecolor="white", bbox_inches="tight")
        print("Saved:", outfile)
        plt.close('all')
    if show:
        plt.show()
    else:
        return ax
# ==================================================== #