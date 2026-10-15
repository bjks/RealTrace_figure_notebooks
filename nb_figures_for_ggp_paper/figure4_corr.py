#!/usr/bin/env python
# coding: utf-8

# # Figure 4: Plot correlation function

# In[1]:


import sys

import builtins
builtins.NAMING = "new"


sys.path.insert(0, "../../")
from RealTrace_paper_analysis.plotting_src.header_settings import *

from RealTrace_paper_analysis.plotting_src.plot_correlations_from_npz import *



# In[2]:


# np.exp(1,  dtype="float96")


# ### Reading files

# In[3]:


def corr_to_csv(correlations, output_file):
    
    ["x(t+dt)", "g(t+dt)", "l(t+dt)", "q(t+dt)", "x(t)", "g(t)", "l(t)", "q(t)"]
    
    
    columns = ["dt", 
               "cov_l(t+dt)l(t)", "cov_l(t+dt)l(t)_err", 
               "cov_l(t+dt)q(t)", "cov_l(t+dt)q(t)_err",
               "cov_q(t+dt)l(t)", "cov_q(t+dt)l(t)_err", 
               "cov_q(t+dt)q(t)", "cov_q(t+dt)q(t)_err", 
               "cov_c(t+dt)c(t)", "cov_c(t+dt)c(t)_err",
               ###
               "corr_l(t+dt)l(t)", "corr_l(t+dt)l(t)_err",
               "corr_l(t+dt)q(t)", "corr_l(t+dt)q(t)_err",
               "corr_q(t+dt)l(t)", "corr_q(t+dt)l(t)_err", 
               "corr_q(t+dt)q(t)", "corr_q(t+dt)q(t)_err", 
               "corr_c(t+dt)c(t)", "corr_c(t+dt)c(t)_err", 
               ###
               "corr_naive_l(t+dt)l(t)",
               "corr_naive_l(t+dt)q(t)",
               "corr_naive_q(t+dt)l(t)", 
               "corr_naive_q(t+dt)q(t)",
               "corr_naive_c(t+dt)c(t)"]
    
    with open(output_file, 'w') as f:
        for i,c in enumerate(columns):
            f.write(c)
            if i<len(columns)-1:
                f.write(",")
        f.write("\n")
        
        for corr in correlations:
            f.write(str(corr.dt)+",")
            
            f.write(str(corr.cov_mle[2, 6])+",")
            f.write(str(corr.cov_mle_err[2, 6])+",")
            
            f.write(str(corr.cov_mle[2, 7])+",")
            f.write(str(corr.cov_mle_err[2, 7])+",")
            
            f.write(str(corr.cov_mle[3, 6])+",")
            f.write(str(corr.cov_mle_err[3, 6])+",")
            
            f.write(str(corr.cov_mle[3, 7])+",")
            f.write(str(corr.cov_mle_err[3, 7])+",")
            
            f.write(str(corr.cov_concentration_mle[0, 1])+",")
            f.write(str(corr.cov_concentration_mle_err[0, 1])+",")
            
            #######
            f.write(str(corr.corr_mle[2, 6])+",")
            f.write(str(corr.corr_mle_err[2, 6])+",")
            
            f.write(str(corr.corr_mle[2, 7])+",")
            f.write(str(corr.corr_mle_err[2, 7])+",")
            
            f.write(str(corr.corr_mle[3, 6])+",")
            f.write(str(corr.corr_mle_err[3, 6])+",")
            
            f.write(str(corr.corr_mle[3, 7])+",")
            f.write(str(corr.corr_mle_err[3, 7])+",")
            
            f.write(str(corr.corr_concentration_mle[0, 1])+",")
            f.write(str(corr.corr_concentration_mle_err[0, 1])+",")
            
            #######
            f.write(str(corr.corr_naive[2, 6])+",")
            f.write(str(corr.corr_naive[2, 7])+",")
            f.write(str(corr.corr_naive[3, 6])+",")
            f.write(str(corr.corr_naive[3, 7])+",")
            f.write(str(corr.corr_concentration_naive[0, 1]))
            
            f.write("\n")
    


input_dir = "../../data_realtrace_paper/ecoli_data/"
# input_dir = "../../RealTrace_test_data/test_new_npz/"



corralations_filenames = get_input_files(input_dir, keyword="correlations", ext=".npz")
promoters_to_plot = ['hi1', 'hi3', 'rplN', 'rrnB']
corralations_filenames = [cf for cf in corralations_filenames if cf.split("/")[-1].split("_")[1]]
# Filter if wanted
corralations_filenames = [cf for cf in corralations_filenames if cf.split("/")[-1].split("_")[1] in promoters_to_plot]
# corralations_filenames


# In[4]:


def get_mean_growth_rate(cells):
    mean_lambda = np.sum(np.concatenate([np.array(cell.mean_l)/np.array(cell.cov_ll) for cell in cells])) \
                /np.sum(np.concatenate([1./np.array(cell.cov_ll) for cell in cells]))
    return mean_lambda


# In[5]:


corr_dict_all = {}
gamma_lambdas = {}
mean_lambdas = {}
mean_qs = {}
gamma_qs = {}
mean_sq_conc = {}

for cf in corralations_filenames:
    # print(cf)
    data = np.load(cf, allow_pickle=True)

    data_set_name = re.split('_f\d', cf.split('/')[-1] )[0]
    # print(data_set_name, end="\t") 
    # print(list(data.keys()))

    corr_dict_all[data_set_name] = Correlation_function(data['correlations'])

    corr_to_csv(data['correlations'], cf.replace(".npz", ".csv"))
    
    # Collect relevant parameters 
    gamma_lambdas[data_set_name] = data['gamma_lambda']
    mean_lambdas[data_set_name] = data['mean_lambda']
    
    gamma_qs[data_set_name] = data['gamma_q']
    mean_qs[data_set_name] = data['mean_q']
    
    prediction_file = cf.replace("correlations.npz", "prediction.csv")
    cells = df2ggp_cells(pd.read_csv(prediction_file, skiprows=header_lines(prediction_file), 
                                 dtype={"parent_id":str, "cell_id":str}))
    
    mean_sq_conc[data_set_name] = get_mean_growth_rate(cells)**2



# In[ ]:


def sort_dict(indict, sort_key):
    out = {}
    for sk in sort_key:
        for k in indict.keys():
            if sk in k:
                out[k] = indict[k]
    return out


# In[ ]:


conditions_to_plot =  [conditions[0],conditions[1], conditions[2]]
corr_dict = sort_dict(corr_dict_all,conditions_to_plot)

gamma_lambdas = sort_dict(gamma_lambdas, conditions_to_plot)
gamma_qs = sort_dict(gamma_qs, conditions_to_plot)
mean_lambdas = sort_dict(mean_lambdas, conditions_to_plot)



# In[ ]:


altcolor_by_condition = {cond: col for cond,col  in zip(conditions, ["#101518", "#5e798f", "#95a9ba", "#adbdc9" ])}


# # By condition

# In[ ]:


def plot_xy_correlation_list(correlation_function_list, x, y, 
                                plot_file=None, 
                                label=[], 
                                gamma=[], 
                                mean_lambda = [],
                                scale_t=False, 
                                normalize=[],
                                xlim=[None, None],
                                ylim=[None,None], 
                                ylabel=None,
                                err_style="bar", 
                                log=False, 
                                min_joint_number=0, 
                                fit=False, 
                                highlight_x0=False,
                                highlight_y0=False,
                                cov=False, 
                                 color_by="condition",
                                legend=True,
                                title=None):

    if cov:
        get_func=Correlation_function.get_cov_mle
        get_err_func=Correlation_function.get_cov_mle_err
    else:
        get_func=Correlation_function.get_corr_mle
        get_err_func=Correlation_function.get_corr_mle_err
        
        
    labels = copy.deepcopy(label)
    gammas = copy.deepcopy(gamma)
    mean_lambdas = copy.deepcopy(mean_lambda)

    # =========== figure =========== #
    fig, ax0 = plt.subplots(1, 1, figsize=(FULL_SIZE/2, FULL_HEIGHT))


    axs = [ax0]

    ax0.set_title(label_by_condition[title])
    for ax in axs:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        # set the x-spine
#         ax.spines['left'].set_position('zero')

        # turn off the right spine/ticks
        ax.spines['right'].set_color('none')
        ax.yaxis.tick_left()

    # "pad" the following lists with Nones
    labels += [None] * len(correlation_function_list)
    gammas = list(gammas)+ [None] * len(correlation_function_list)
    mean_lambdas = list(mean_lambdas) + [None] * len(correlation_function_list)
    normalize = list(normalize) + [None] * len(correlation_function_list)

    new_ylim = [np.inf,0]
    # =========== plot =========== #
    for i, cf in enumerate(correlation_function_list):
        condition = labels[i].split('_')[0]
        promoter = labels[i].split('_')[1]
      
        color = altcolor_by_condition[condition]

        correlation_function = copy.deepcopy(cf)
        correlation_function.filter_by_n(min_joint_number)

        gamma = gammas[i]
        mean_lambda = mean_lambdas[i]
        norm = normalize[i]

        # =========== correlation =========== #
        dts = correlation_function.dt
        rs = get_func(correlation_function, x,y)
        errs = get_err_func(correlation_function, x,y)
        if norm != None:
            rs /= norm
            errs /= norm

        ### Scale dt by mean doubling time 
        if scale_t:
            dts /= (np.log(2.)/mean_lambda)
            if gamma!=None:
                gamma*=(np.log(2)/mean_lambda)
        else:
            if mean_lambda!=None:
                ax0.axvline(np.log(2)/mean_lambda, color=color, alpha=0.5)

        if gamma!=None:
            r0 = rs[dts==0][0]
            ax0.plot(dts, r0*np.exp(-dts*gamma), ls='--', color=color, alpha=1)

        if xlim!=[None,None]:
            dt_filter = (dts<=xlim[-1])*(dts>=xlim[0])
            dts=dts[dt_filter]
            rs=rs[dt_filter]
            errs=errs[dt_filter]
        ### correlation plot with \without error bars
        if err_style == "bar":
            ax0.errorbar(dts, rs, yerr=errs, lw=2, color=color, label=label)
        elif err_style == "fill":
            ax0.plot(dts, rs, lw=2, color=color, label=label)
            ax0.fill_between(dts, rs-errs, rs+errs, color=color, alpha=0.4)
        elif err_style == None:
            ax0.plot(dts, rs, lw=2, color=color, label=label)
        else:
            ax0.plot(dts, dts*0, color=color)  
            
    # ===== layout ===== #
    if ylim != [None, None]:
        new_ylim = ylim
    else:
        new_ylim = [None, None]

    ax0.set_ylim(new_ylim)

        
    if ylabel == None:
        ax0.set_ylabel(r"$\langle {:s}, {:s}\rangle$".format(x,y))
    else:
        ax0.set_ylabel(ylabel)
    
    if log:
        ax0.set_yscale('log')
        
        yticks = np.ravel([np.array([1, 5])*(10.**i) for i in np.arange(-10,10)])
        yticks = yticks[(yticks>=new_ylim[0]) * (yticks<=new_ylim[1])]
        
        ax0.set_yticks(yticks)
        ax0.set_yticklabels(yticks)
        
    ax0.axvline(0.25, color="tab:grey", lw=2, ls=':')
    ax0.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(base=0.25))

    ax0.set_xlim(xlim)

    if legend:
        ax0.legend(*legend_without_duplicate_labels(ax0), loc='lower center', bbox_to_anchor=(0.5, 1), ncol=2)

    ax0.set_xlabel(r'time lag $(t^\prime-t)/$(mdt)')
    
    # plt.grid(axis="y")
    
    fig.tight_layout()
    if plot_file != None:
        print("Saved in", plot_file)
        fig.savefig(plot_file, dpi=300, facecolor="white", bbox_inches="tight")
        plt.show()
        
    else:
        plt.show()
    plt.close()
    
    
    
plot_dir_main = os.path.join("../../figures_for_method_paper/plots_main/correlations" ,"{:s}_{:s}_{:s}_{:s}.pdf")
    
### ACETATE
condition=conditions_to_plot[0]

plot_xy_correlation_list(get_filtered_values(corr_dict, condition), 
                        "l(t+dt)", "l(t)", 
                        plot_file=plot_dir_main.format(condition, "corr", "l", "l"), 
                        label=get_filtered_keys(corr_dict, condition, keep="both"), 
                        gamma=get_filtered_values(gamma_lambdas, condition), 
                        mean_lambda = get_filtered_values(mean_lambdas, condition),
                        ylabel=r"$C[\lambda(t^\prime),\lambda(t)]$",
                        scale_t=True, 
                        ylim=[0.04,1], 
                        xlim=[0, 2], 
                        err_style="fill", 
                        log=True, 
                        min_joint_number=1, 
                        legend=False, 
                        color_by="condition",
                        title=condition)


### GLYCEROL
condition=conditions_to_plot[1]

plot_xy_correlation_list(get_filtered_values(corr_dict, condition), 
                        "l(t+dt)", "l(t)", 
                        plot_file=plot_dir_main.format(condition, "corr", "l", "l"), 
                        label=get_filtered_keys(corr_dict, condition, keep="both"), 
                        gamma=get_filtered_values(gamma_lambdas, condition), 
                        mean_lambda = get_filtered_values(mean_lambdas, condition),
                        ylabel=r"$C[\lambda(t^\prime),\lambda(t)]$",
                        scale_t=True, 
                        ylim=[0.04,1], 
                        xlim=[0, 2], 
                        err_style="fill", 
                        log=True, 
                        min_joint_number=1, 
                        legend=False, 
                        color_by="condition",
                        title=condition)

### GLUCOSE
condition=conditions_to_plot[2]

plot_xy_correlation_list(get_filtered_values(corr_dict, condition), 
                        "l(t+dt)", "l(t)", 
                        plot_file=plot_dir_main.format(condition, "corr", "l", "l"), 
                        label=get_filtered_keys(corr_dict, condition, keep="both"), 
                        gamma=get_filtered_values(gamma_lambdas, condition), 
                        mean_lambda = get_filtered_values(mean_lambdas, condition),
                        ylabel=r"$C[\lambda(t^\prime),\lambda(t)]$",
                        scale_t=True, 
                        ylim=[0.04,1], 
                        xlim=[0, 2], 
                        err_style="fill", 
                        log=True, 
                        min_joint_number=1, 
                        legend=False, 
                        color_by="condition",
                        title=condition)


for condition in conditions_to_plot[:3]:

    plot_xy_correlation_list(get_filtered_values(corr_dict, condition), 
                            "l(t+dt)", "l(t)", 
                            plot_file=plot_dir_main.format(condition, "corr", "l", "l").replace("corr_l_l", "corr_l_l_long4"), 
                            label=get_filtered_keys(corr_dict, condition, keep="both"), 
                            gamma=get_filtered_values(gamma_lambdas, condition), 
                            mean_lambda = get_filtered_values(mean_lambdas, condition),
                            ylabel=r"$C[\lambda(t^\prime),\lambda(t)]$",
                            scale_t=True, 
                            ylim=[-0.2,1], 
                            xlim=[0, 4], 
                            err_style="fill", 
                            log=False, 
                            min_joint_number=1, 
                            legend=False, 
                            color_by="condition",
                            title=condition)


    plot_xy_correlation_list(get_filtered_values(corr_dict, condition), 
                            "l(t+dt)", "l(t)", 
                            plot_file=plot_dir_main.format(condition, "corr", "l", "l").replace("corr_l_l", "corr_l_l_long6"), 
                            label=get_filtered_keys(corr_dict, condition, keep="both"), 
                            gamma=get_filtered_values(gamma_lambdas, condition), 
                            mean_lambda = get_filtered_values(mean_lambdas, condition),
                            ylabel=r"$C[\lambda(t^\prime),\lambda(t)]$",
                            scale_t=True, 
                            ylim=[-0.2,1], 
                            xlim=[0, 6], 
                            err_style="fill", 
                            log=False, 
                            min_joint_number=1, 
                            legend=False, 
                            color_by="condition",
                            title=condition)



    plot_xy_correlation_list(get_filtered_values(corr_dict, condition), 
                            "l(t+dt)", "l(t)", 
                            plot_file=plot_dir_main.format(condition, "corr", "l", "l").replace("corr_l_l", "corr_l_l_long7"), 
                            label=get_filtered_keys(corr_dict, condition, keep="both"), 
                            gamma=get_filtered_values(gamma_lambdas, condition), 
                            mean_lambda = get_filtered_values(mean_lambdas, condition),
                            ylabel=r"$C[\lambda(t^\prime),\lambda(t)]$",
                            scale_t=True, 
                            ylim=[-0.2,1], 
                            xlim=[0, 7], 
                            err_style="fill", 
                            log=False, 
                            min_joint_number=1, 
                            legend=False, 
                            color_by="condition",
                            title=condition)

    plot_xy_correlation_list(get_filtered_values(corr_dict, condition), 
                            "l(t+dt)", "l(t)", 
                            plot_file=plot_dir_main.format(condition, "corr", "l", "l").replace("corr_l_l", "corr_l_l_long8"), 
                            label=get_filtered_keys(corr_dict, condition, keep="both"), 
                            gamma=get_filtered_values(gamma_lambdas, condition), 
                            mean_lambda = get_filtered_values(mean_lambdas, condition),
                            ylabel=r"$C[\lambda(t^\prime),\lambda(t)]$",
                            scale_t=True, 
                            ylim=[-0.2,1], 
                            xlim=[0, 8], 
                            err_style="fill", 
                            log=False, 
                            min_joint_number=1, 
                            legend=False, 
                            color_by="condition",
                            title=condition)


    plot_xy_correlation_list(get_filtered_values(corr_dict, condition), 
                            "l(t+dt)", "l(t)", 
                            plot_file=plot_dir_main.format(condition, "corr", "l", "l").replace("corr_l_l", "corr_l_l_long9"), 
                            label=get_filtered_keys(corr_dict, condition, keep="both"), 
                            gamma=get_filtered_values(gamma_lambdas, condition), 
                            mean_lambda = get_filtered_values(mean_lambdas, condition),
                            ylabel=r"$C[\lambda(t^\prime),\lambda(t)]$",
                            scale_t=True, 
                            ylim=[-0.2,1], 
                            xlim=[0, 9], 
                            err_style="fill", 
                            log=False, 
                            min_joint_number=1, 
                            legend=False, 
                            color_by="condition",
                            title=condition)

