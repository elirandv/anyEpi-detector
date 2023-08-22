import matplotlib.pyplot as plt
import pandas as pd


def scatterplot_current_time_by_base(location, sequence):
    barcodes = [
        {"num": "09", "color": 'b'},
        {"num": "10", "color": 'orange'},
        {"num": "11", "color": 'g'},
        {"num": "12", "color": 'r'}
    ]
    N = location.split('_')
    N = int(N[0])
    bc11 = pd.read_csv(
        r"./220299_5hmc_training/barcode11_epi_basecaller_outputs/101217-3303900/aligned_reads_barcode11_forward_" + location + ".eventalign.tsv",
        header=None, sep="\t")
    bc10 = pd.read_csv(
        r"./220299_5hmc_training/barcode10_epi_basecaller_outputs/101217-3303900/aligned_reads_barcode10_forward_" + location + ".eventalign.tsv",
        header=None, sep="\t")
    bc09 = pd.read_csv(
        r"./220299_5hmc_training/barcode09_epi_basecaller_outputs/101217-3303900/aligned_reads_barcode09_forward_" + location + ".eventalign.tsv",
        header=None, sep="\t")
    bc12 = pd.read_csv(
        r"./220299_5hmc_training/barcode12_epi_basecaller_outputs/101217-3303900/aligned_reads_barcode12_forward_" + location + ".eventalign.tsv",
        header=None, sep="\t")
    bc09["ones"] = 1
    bc10["ones"] = 1
    bc11["ones"] = 1
    bc12["ones"] = 1
    grouped09 = bc09.groupby([1, 3]).agg({8: 'sum', "ones": 'sum', 6: 'mean'})
    grouped10 = bc10.groupby([1, 3]).agg({8: 'sum', "ones": 'sum', 6: 'mean'})
    grouped11 = bc11.groupby([1, 3]).agg({8: 'sum', "ones": 'sum', 6: 'mean'})
    grouped12 = bc12.groupby([1, 3]).agg({8: 'sum', "ones": 'sum', 6: 'mean'})
    grouped09 = grouped09.reset_index()
    grouped10 = grouped10.reset_index()
    grouped11 = grouped11.reset_index()
    grouped12 = grouped12.reset_index()
    mapped09 = grouped09.groupby([1]).agg({8: 'mean', "ones": 'mean', 6: 'mean'})
    for index, row in grouped09.iterrows():
        grouped09.at[index, "current diff"] = (grouped09.loc[index, 6]) - (mapped09.loc[(grouped09.iloc[index, 0]), 6])
        grouped09.at[index, "dwell diff"] = (grouped09.loc[index, 8]) - (mapped09.loc[(grouped09.iloc[index, 0]), 8])
    for index, row in grouped10.iterrows():
        grouped10.at[index, "current diff"] = (grouped10.loc[index, 6]) - (mapped09.loc[(grouped10.iloc[index, 0]), 6])
        grouped10.at[index, "dwell diff"] = (grouped10.loc[index, 8]) - (mapped09.loc[(grouped10.iloc[index, 0]), 8])
    for index, row in grouped11.iterrows():
        grouped11.at[index, "current diff"] = (grouped11.loc[index, 6]) - (mapped09.loc[(grouped11.iloc[index, 0]), 6])
        grouped11.at[index, "dwell diff"] = (grouped11.loc[index, 8]) - (mapped09.loc[(grouped11.iloc[index, 0]), 8])
    for index, row in grouped12.iterrows():
        grouped12.at[index, "current diff"] = (grouped12.loc[index, 6]) - (mapped09.loc[(grouped12.iloc[index, 0]), 6])
        grouped12.at[index, "dwell diff"] = (grouped12.loc[index, 8]) - (mapped09.loc[(grouped12.iloc[index, 0]), 8])
    for i, g in grouped09.groupby(1):
        globals()['grouped09_' + str(i)] = g
    for i, g in grouped10.groupby(1):
        globals()['grouped10_' + str(i)] = g
    for i, g in grouped11.groupby(1):
        globals()['grouped11_' + str(i)] = g
    for i, g in grouped12.groupby(1):
        globals()['grouped12_' + str(i)] = g
        print(str(i))
    # Put here the string as in excle, in each of the following line that comes after "#" line, you need to change the base according to this sequence order:

    fig, axs = plt.subplots(2, 5, figsize=(35, 10))

    for i in range(10):
        if i < 5:
            loc = (0, i)
        else:
            loc = (1, i - 5)
        cur_axs = axs[loc[0], loc[1]]
        for bc in barcodes:
            cur_axs.scatter(eval("grouped" + bc["num"] + "_" + str(N + i + 1))['current diff'],
                            eval("grouped" + bc["num"] + "_" + str(N + i + 1))['dwell diff'],
                            color=bc["color"], alpha=0.2, s=5, linewidths=None)
        cur_axs.set_ylim(-0.005, 0.02)
        cur_axs.set_xlim(-15, 10)
        cur_axs.set_ylabel("delta_time")
        cur_axs.set_xlabel("delta_current")
        #
        cur_axs.set_title((sequence[:i] + r'$\bf{{{' + sequence[i] + '}}}$' + sequence[i + 1:]))
        cur_axs.axhline(y=0, color='k', linewidth=.5)
        cur_axs.axvline(x=0, color='k', linewidth=.5)

    plt.show()


if __name__ == "__main__":
    #Put here name of chunk:
    location = "2201053_2201072"
    sequence = "CTACATGTCGCATG"
    scatterplot_current_time_by_base(location, sequence)
    #
    # def scatterplot_current_time_by_base(location, sequence):
    #     barcodes = [
    #         {"num": "09", "color": 'b'},
    #         {"num": "10", "color": 'orange'},
    #         {"num": "11", "color": 'g'},
    #         {"num": "12", "color": 'r'}
    #     ]
    #     N = location.split('_')
    #     N = int(N[0])
    #     bc11 = pd.read_csv(r"./220299_5hmc_training/barcode11_epi_basecaller_outputs/101217-3303900/aligned_reads_barcode11_forward_" + location + ".eventalign.tsv", header = None, sep ="\t")
    #     bc10 = pd.read_csv(r"./220299_5hmc_training/barcode10_epi_basecaller_outputs/101217-3303900/aligned_reads_barcode10_forward_" + location + ".eventalign.tsv", header = None, sep ="\t")
    #     bc09 = pd.read_csv(r"./220299_5hmc_training/barcode09_epi_basecaller_outputs/101217-3303900/aligned_reads_barcode09_forward_" + location + ".eventalign.tsv", header = None, sep ="\t")
    #     bc12 = pd.read_csv(r"./220299_5hmc_training/barcode12_epi_basecaller_outputs/101217-3303900/aligned_reads_barcode12_forward_" + location + ".eventalign.tsv", header = None, sep ="\t")
    #     bc09["ones"] = 1
    #     bc10["ones"] = 1
    #     bc11["ones"] = 1
    #     bc12["ones"] = 1
    #     grouped09 = bc09.groupby([1,3]).agg({8:'sum', "ones":'sum',6:'mean'})
    #     grouped10 = bc10.groupby([1,3]).agg({8:'sum', "ones":'sum',6:'mean'})
    #     grouped11 = bc11.groupby([1,3]).agg({8:'sum', "ones":'sum',6:'mean'})
    #     grouped12 = bc12.groupby([1,3]).agg({8:'sum', "ones":'sum',6:'mean'})
    #     grouped09 = grouped09.reset_index()
    #     grouped10 = grouped10.reset_index()
    #     grouped11 = grouped11.reset_index()
    #     grouped12 = grouped12.reset_index()
    #     mapped09 = grouped09.groupby([1]).agg({8:'mean', "ones":'mean',6:'mean'})
    #     for index, row in grouped09.iterrows():
    #         grouped09.at[index,"current diff"] = (grouped09.loc[index,6])-(mapped09.loc[(grouped09.iloc[index,0]),6])
    #         grouped09.at[index,"dwell diff"] = (grouped09.loc[index,8])-(mapped09.loc[(grouped09.iloc[index,0]),8])
    #     for index, row in grouped10.iterrows():
    #         grouped10.at[index,"current diff"] = (grouped10.loc[index,6])-(mapped09.loc[(grouped10.iloc[index,0]),6])
    #         grouped10.at[index,"dwell diff"] = (grouped10.loc[index,8])-(mapped09.loc[(grouped10.iloc[index,0]),8])
    #     for index, row in grouped11.iterrows():
    #         grouped11.at[index,"current diff"] = (grouped11.loc[index,6])-(mapped09.loc[(grouped11.iloc[index,0]),6])
    #         grouped11.at[index,"dwell diff"] = (grouped11.loc[index,8])-(mapped09.loc[(grouped11.iloc[index,0]),8])
    #     for index, row in grouped12.iterrows():
    #         grouped12.at[index,"current diff"] = (grouped12.loc[index,6])-(mapped09.loc[(grouped12.iloc[index,0]),6])
    #         grouped12.at[index,"dwell diff"] = (grouped12.loc[index,8])-(mapped09.loc[(grouped12.iloc[index,0]),8])
    #     for i, g in grouped09.groupby(1):
    #         globals()['grouped09_' + str(i)] = g
    #     for i, g in grouped10.groupby(1):
    #         globals()['grouped10_' + str(i)] = g
    #     for i, g in grouped11.groupby(1):
    #         globals()['grouped11_' + str(i)] = g
    #     for i, g in grouped12.groupby(1):
    #         globals()['grouped12_' + str(i)] = g
    #         print(str(i))
    #     #Put here the string as in excle, in each of the following line that comes after "#" line, you need to change the base according to this sequence order:
    #
    #     fig, axs = plt.subplots(2, 5, figsize=(35, 10))
    #
    #     for i in range(10):
    #         if i < 5:
    #             loc = (0, i)
    #         else:
    #             loc = (1, i-5)
    #         cur_axs = axs[loc[0], loc[1]]
    #         for bc in barcodes:
    #             cur_axs.scatter(eval("grouped"+bc["num"]+"_"+str(N+i+1))['current diff'],
    #                             eval("grouped"+bc["num"]+"_"+str(N+i+1))['dwell diff'],
    #                             color=bc["color"], alpha=0.2, s=5, linewidths=None)
    #         cur_axs.set_ylim(-0.005, 0.02)
    #         cur_axs.set_xlim(-15, 10)
    #         cur_axs.set_ylabel("delta_time")
    #         cur_axs.set_xlabel("delta_current")
    #         #
    #         cur_axs.set_title((sequence[:i] + r'$\bf{{{' + sequence[i] + '}}}$' + sequence[i + 1:]))
    #         cur_axs.axhline(y=0, color='k', linewidth=.5)
    #         cur_axs.axvline(x=0, color='k', linewidth=.5)
    #
    #     plt.show()
    #     exit()
    #
    #     bold_index  =  1
    #     axs[0,1].scatter(eval("grouped09_"+str(N+2))['current diff'],eval("grouped09_"+str(N+2))['dwell diff'], color = 'b', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,1].scatter(eval("grouped10_"+str(N+2))['current diff'],eval("grouped10_"+str(N+2))['dwell diff'], color = 'orange', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,1].scatter(eval("grouped11_"+str(N+2))['current diff'],eval("grouped11_"+str(N+2))['dwell diff'], color = 'g', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,1].scatter(eval("grouped12_"+str(N+2))['current diff'],eval("grouped12_"+str(N+2))['dwell diff'], color = 'r', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,1].set_ylim(-0.005,0.02)
    #     axs[0,1].set_xlim(-15,10)
    #     axs[0,1].set_ylabel("delta_time")
    #     axs[0,1].set_xlabel("delta_current")
    #     #
    #     axs[0,1].set_title((sequence[:bold_index] + r'$\bf{{{' + sequence[bold_index] + '}}}$' + sequence[bold_index + 1:]))
    #     axs[0,1].axhline(y = 0, color = 'k', linewidth = .5)
    #     axs[0,1].axvline(x = 0, color = 'k', linewidth = .5)
    #     bold_index  =  2
    #     axs[0,2].scatter(eval("grouped09_"+str(N+3))['current diff'],eval("grouped09_"+str(N+3))['dwell diff'], color = 'b', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,2].scatter(eval("grouped10_"+str(N+3))['current diff'],eval("grouped10_"+str(N+3))['dwell diff'], color = 'orange', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,2].scatter(eval("grouped11_"+str(N+3))['current diff'],eval("grouped11_"+str(N+3))['dwell diff'], color = 'g', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,2].scatter(eval("grouped12_"+str(N+3))['current diff'],eval("grouped12_"+str(N+3))['dwell diff'], color = 'r', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,2].set_ylim(-0.005,0.02)
    #     axs[0,2].set_xlim(-15,10)
    #     axs[0,2].set_ylabel("delta_time")
    #     axs[0,2].set_xlabel("delta_current")
    #     #
    #     axs[0,1].set_title((sequence[:bold_index] + r'$\bf{{{' + sequence[bold_index] + '}}}$' + sequence[bold_index + 1:]))
    #     axs[0,2].axhline(y = 0, color = 'k', linewidth = .5)
    #     axs[0,2].axvline(x = 0, color = 'k', linewidth = .5)
    #     bold_index  =  3
    #     axs[0,3].scatter(eval("grouped09_"+str(N+4))['current diff'],eval("grouped09_"+str(N+4))['dwell diff'], color = 'b', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,3].scatter(eval("grouped10_"+str(N+4))['current diff'],eval("grouped10_"+str(N+4))['dwell diff'], color = 'orange', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,3].scatter(eval("grouped11_"+str(N+4))['current diff'],eval("grouped11_"+str(N+4))['dwell diff'], color = 'g', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,3].scatter(eval("grouped12_"+str(N+4))['current diff'],eval("grouped12_"+str(N+4))['dwell diff'], color = 'r', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,3].set_ylim(-0.005,0.02)
    #     axs[0,3].set_xlim(-15,10)
    #     axs[0,3].set_ylabel("delta_time")
    #     axs[0,3].set_xlabel("delta_current")
    #     #
    #     axs[0,3].set_title((sequence[:bold_index] + r'$\bf{{{C}}}$' + sequence[bold_index + 1:]))
    #     axs[0,3].axhline(y = 0, color = 'k', linewidth = .5)
    #     axs[0,3].axvline(x = 0, color = 'k', linewidth = .5)
    #     bold_index  =  4
    #     axs[0,4].scatter(eval("grouped09_"+str(N+5))['current diff'],eval("grouped09_"+str(N+5))['dwell diff'], color = 'b', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,4].scatter(eval("grouped10_"+str(N+5))['current diff'],eval("grouped10_"+str(N+5))['dwell diff'], color = 'orange', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,4].scatter(eval("grouped11_"+str(N+5))['current diff'],eval("grouped11_"+str(N+5))['dwell diff'], color = 'g', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,4].scatter(eval("grouped12_"+str(N+5))['current diff'],eval("grouped12_"+str(N+5))['dwell diff'], color = 'r', alpha = 0.2, s = 5, linewidths = None)
    #     axs[0,4].set_ylim(-0.005,0.02)
    #     axs[0,4].set_xlim(-15,10)
    #     axs[0,4].set_ylabel("delta_time")
    #     axs[0,4].set_xlabel("delta_current")
    #     #
    #     axs[0,4].set_title((sequence[:bold_index] + r'$\bf{{{A}}}$' + sequence[bold_index + 1:]))
    #     axs[0,4].axhline(y = 0, color = 'k', linewidth = .5)
    #     axs[0,4].axvline(x = 0, color = 'k', linewidth = .5)
    #     bold_index  =  5
    #     axs[1,0].scatter(eval("grouped09_"+str(N+6))['current diff'],eval("grouped09_"+str(N+6))['dwell diff'], color = 'b', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,0].scatter(eval("grouped10_"+str(N+6))['current diff'],eval("grouped10_"+str(N+6))['dwell diff'], color = 'orange', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,0].scatter(eval("grouped11_"+str(N+6))['current diff'],eval("grouped11_"+str(N+6))['dwell diff'], color = 'g', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,0].scatter(eval("grouped12_"+str(N+6))['current diff'],eval("grouped12_"+str(N+6))['dwell diff'], color = 'r', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,0].set_ylim(-0.005,0.02)
    #     axs[1,0].set_xlim(-15,10)
    #     axs[1,0].set_ylabel("delta_time")
    #     axs[1,0].set_xlabel("delta_current")
    #     #
    #     axs[1,0].set_title((sequence[:bold_index] + r'$\bf{{{T}}}$' + sequence[bold_index + 1:]))
    #     axs[1,0].axhline(y = 0, color = 'k', linewidth = .5)
    #     axs[1,0].axvline(x = 0, color = 'k', linewidth = .5)
    #     bold_index  =  6
    #     axs[1,1].scatter(eval("grouped09_"+str(N+7))['current diff'],eval("grouped09_"+str(N+7))['dwell diff'], color = 'b', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,1].scatter(eval("grouped10_"+str(N+7))['current diff'],eval("grouped10_"+str(N+7))['dwell diff'], color = 'orange', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,1].scatter(eval("grouped11_"+str(N+7))['current diff'],eval("grouped11_"+str(N+7))['dwell diff'], color = 'g', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,1].scatter(eval("grouped12_"+str(N+7))['current diff'],eval("grouped12_"+str(N+7))['dwell diff'], color = 'r', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,1].set_ylim(-0.005,0.02)
    #     axs[1,1].set_xlim(-15,10)
    #     axs[1,1].set_ylabel("delta_time")
    #     axs[1,1].set_xlabel("delta_current")
    #     #
    #     axs[1,1].set_title((sequence[:bold_index] + r'$\bf{{{G}}}$' + sequence[bold_index + 1:]))
    #     axs[1,1].axhline(y = 0, color = 'k', linewidth = .5)
    #     axs[1,1].axvline(x = 0, color = 'k', linewidth = .5)
    #     bold_index  =  7
    #     axs[1,2].scatter(eval("grouped09_"+str(N+8))['current diff'],eval("grouped09_"+str(N+8))['dwell diff'], color = 'b', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,2].scatter(eval("grouped10_"+str(N+8))['current diff'],eval("grouped10_"+str(N+8))['dwell diff'], color = 'orange', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,2].scatter(eval("grouped11_"+str(N+8))['current diff'],eval("grouped11_"+str(N+8))['dwell diff'], color = 'g', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,2].scatter(eval("grouped12_"+str(N+8))['current diff'],eval("grouped12_"+str(N+8))['dwell diff'], color = 'r', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,2].set_ylim(-0.005,0.02)
    #     axs[1,2].set_xlim(-15,10)
    #     axs[1,2].set_ylabel("delta_time")
    #     axs[1,2].set_xlabel("delta_current")
    #     #
    #     axs[1,2].set_title((sequence[:bold_index] + r'$\bf{{{T}}}$' + sequence[bold_index + 1:]))
    #     axs[1,2].axhline(y = 0, color = 'k', linewidth = .5)
    #     axs[1,2].axvline(x = 0, color = 'k', linewidth = .5)
    #     bold_index  =  8
    #     axs[1,3].scatter(eval("grouped09_"+str(N+9))['current diff'],eval("grouped09_"+str(N+9))['dwell diff'], color = 'b', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,3].scatter(eval("grouped10_"+str(N+9))['current diff'],eval("grouped10_"+str(N+9))['dwell diff'], color = 'orange', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,3].scatter(eval("grouped11_"+str(N+9))['current diff'],eval("grouped11_"+str(N+9))['dwell diff'], color = 'g', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,3].scatter(eval("grouped12_"+str(N+9))['current diff'],eval("grouped12_"+str(N+9))['dwell diff'], color = 'r', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,3].set_ylim(-0.005,0.02)
    #     axs[1,3].set_xlim(-15,10)
    #     axs[1,3].set_ylabel("delta_time")
    #     axs[1,3].set_xlabel("delta_current")
    #     #
    #     axs[1,3].set_title((sequence[:bold_index] + r'$\bf{{{C}}}$' + sequence[bold_index + 1:]))
    #     axs[1,3].axhline(y = 0, color = 'k', linewidth = .5)
    #     axs[1,3].axvline(x = 0, color = 'k', linewidth = .5)
    #
    #     bold_index  =  9
    #     axs[1,4].scatter(eval("grouped09_"+str(N+10))['current diff'],eval("grouped09_"+str(N+10))['dwell diff'], color = 'b', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,4].scatter(eval("grouped10_"+str(N+10))['current diff'],eval("grouped10_"+str(N+10))['dwell diff'], color = 'orange', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,4].scatter(eval("grouped11_"+str(N+10))['current diff'],eval("grouped11_"+str(N+10))['dwell diff'], color = 'g', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,4].scatter(eval("grouped12_"+str(N+10))['current diff'],eval("grouped12_"+str(N+10))['dwell diff'], color = 'r', alpha = 0.2, s = 5, linewidths = None)
    #     axs[1,4].set_ylim(-0.005,0.02)
    #     axs[1,4].set_xlim(-15,10)
    #     axs[1,4].set_ylabel("delta_time")
    #     axs[1,4].set_xlabel("delta_current")
    #     axs[1,4].set_title((sequence[:bold_index] + r'$\bf{{{G}}}$' + sequence[bold_index + 1:]))
    #     axs[1,4].axhline(y = 0, color = 'k', linewidth = .5)
    #     axs[1,4].axvline(x = 0, color = 'k', linewidth = .5)
    #
    #     plt.show()
