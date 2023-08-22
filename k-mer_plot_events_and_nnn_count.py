# importing the required module
import os
import re
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from seabornfig2grid import SeabornFig2Grid as sfg

from functools import wraps
import time
import math


def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f'Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds')
        return result

    return timeit_wrapper


@timeit
def calculate_something(num):
    """
    Simple function that returns sum of all numbers up to the square of num.
    """
    total = sum((x for x in range(0, num ** 2)))
    return total


def scatterplot_current_time_by_base(args: dict, back: int, front: int):
    sns.set()
    assert not args.back_window < back, 'back window limit exceeded'
    assert not 9 < front, 'front window limit exceeded'
    seq_location = args.location
    sequence = args.sequence
    output_file = args.output
    print(f"scatterplot_current_time_by_base:(seq_location={seq_location}, sequence={sequence})")
    full_dict = {"09": {"i": 0, "color": 'b', "name": "Cytosine"},
                 "10": {"i": 1, "color": 'orange', "name": "5hmC"},
                 "11": {"i": 2, "color": 'g', "name": "5gmc"},
                 "12": {"i": 3, "color": 'r', "name": "5gmc-N3"}}
    barcode_dict = {key: full_dict[key] for key in full_dict if key in args.barcodes}
    start_loc = int(seq_location.split('_')[0])
    mod_loc = start_loc + args.back_window 1
    bc_df = {}
    grouped = {}
    for bc in args.barcodes:
        bc_output_dir = f'{args.output}/aligned_reads_barcode{bc}_{args.output.split("/")[-1]}_{args.strand}/'
        file = f'{bc_output_dir}{args.location}.eventalign.tsv'
        bc_df[bc] = pd.read_csv(file, header=None, sep="\t", usecols=[1, 3, 6, 7, 8, 13])
        bc_df[bc].columns = ['pos', 'read', 'curr', 'std', 'd_time', 'samples']

        bc_df[bc] = bc_df[bc][(bc_df[bc]['pos'] >= mod_loc - back) & (bc_df[bc]['pos'] < mod_loc + front)]

        bc_df[bc]['ones'] = 1

        ### Jonathans window agg suggestion!!! -8 to +3 known pos to known pos
        hole_window_agg[bc] = bc_df[bc].groupby(['read']).agg({'curr': 'mean', 'd_time': 'sum', 'samples': lambda x: ','.join(x)})
        hole_window_agg[bc]['samples'] = hole_window_agg[bc]['samples'].str.split(",")
        hole_window_agg[bc]['samples_count'] = hole_window_agg[bc]['samples'].str.len()
        hole_window_agg[bc]['samples_count'].value_counts().sort_values().plot(kind='bar')\
            .savefig(f'{bc_output_dir}/{seq_location}_samples_count_hist')

        hole_window_agg09 = hole_window_agg[min(args.barcodes)].groupby(lambda _: True)\
            .agg({'curr': 'mean', 'd_time': 'sum'})
        hole_window_agg[bc]['curr diff'] = hole_window_agg[bc]['curr']-hole_window_agg09['curr'][1]
        hole_window_agg[bc]['d_time diff'] = hole_window_agg[bc]['d_time']-hole_window_agg09['d_time'][1]
        # cur_axs.scatter(hole_window_agg[bc_num]['curr diff'], hole_window_agg[bc_num]['d_time diff'],
        #     #                             color=bc_dict["color"], alpha=0.2, s=5, linewidths=None)
        ### make samples_count & mean_curr distributions!!!

        grouped[bc] = bc_df[bc].groupby(['pos', 'read']).agg({'curr': 'mean', 'd_time': 'sum', 'ones': 'sum'})
        grouped[bc] = grouped[bc].reset_index()

    mapped09 = grouped[min(args.barcodes)].groupby(['pos']).agg({'curr': 'mean', 'd_time': 'mean', 'ones': 'mean'})

    def get_diff(row, col_name):
        return row[col_name] - mapped09.loc[(row['pos']), [col_name]]

    df = pd.concat(grouped)

    @timeit
    def diffs():
        df['curr_diff'] = df.apply(lambda row: get_diff(row, 'curr'), axis=1)
        df['d_time_diff'] = df.apply(lambda row: get_diff(row, 'd_time'), axis=1)

    @timeit
    def diffs2():
        for bc in args.barcodes:
            for index, row in df.iterrows():
                grouped[bc].at[index, "current diff"] = (grouped[bc].loc[index, 'curr']) - (
                mapped09.loc[(grouped[bc].iloc[index, 0]), 'curr'])
                grouped[bc].at[index, "dwell diff"] = (grouped[bc].loc[index, 'd_time']) - (
                mapped09.loc[(grouped[bc].iloc[index, 0]), 'd_time'])
        # for index, row in grouped[bc].iterrows():
        #     grouped[bc].at[index, "current diff"] = (grouped[bc].loc[index, 'curr']) - (mapped09.loc[(grouped[bc].iloc[index, 0]), 'curr'])
        #     grouped[bc].at[index, "dwell diff"] = (grouped[bc].loc[index, 'd_time']) - (mapped09.loc[(grouped[bc].iloc[index, 0]), 'd_time'])
        # for i, g in grouped[bc].groupby(1):
        #     globals()[f'grouped{bc}_' + str(i)] = g

    diffs()
    # diffs2()

    positions = [i for i in range(mod_loc - back, mod_loc + front)]

    fig = plt.figure(figsize=(24, 16))
    gs = gridspec.GridSpec(math.ceil(front + back / 4), 4)
    i = 0
    for pos in positions:
        plot = sns.jointplot(data=df.loc[df['pos'] == pos], x='curr_diff', y='d_time_diff', alpha=0.3, s=5)
        #, linewidths=None, legend=False, marginal_kws={'common_norm': False}, hue="Modification")
        plot.ax_joint.axvline(x=0, color='k', linewidth=.5)
        plot.ax_joint.axhline(y=0, color='k', linewidth=.5)
        plot.ax_joint.set_xlim(-15, 10)
        plot.ax_joint.set_ylim(-0.005, 0.02)
        # plot.set_title(pos)
        # plot.tick_params(axis='both', which='major', labelsize=8)
        mgkey = sfg(plot, fig, gs[i])
        i = i + 1
        print(pos)
    # gs = gs.annotate(stats.pearsonr, fontsize=10)
    gs.tight_layout(fig)
    # for j in range(2):
    #     fig, axs = plt.subplots(1, 5, figsize=(30, 5))
    #
    #     for i in range(5):
    #         if j < 1:
    #             loc = i
    #         else:
    #             loc = i + 5
    #         cur_axs = axs[i]
    #         for bc_num, bc_dict in barcode_dict.items():
    #
    #             # plot scatter
    #             cur_axs.scatter(eval("grouped" + bc_num + "_" + str(start_loc + loc + 1))['current diff'],
    #                             eval("grouped" + bc_num + "_" + str(start_loc + loc + 1))['dwell diff'],
    #                             color=bc_dict["color"], alpha=0.2, s=5, linewidths=None)
    #
    #         cur_axs.set_ylim(-0.005, 0.02)
    #         cur_axs.set_xlim(-15, 10)
    #         cur_axs.set_ylabel("delta_time")
    #         cur_axs.set_xlabel("delta_current")
    #         cur_axs.set_title((sequence[:loc] + r'$\bf{{{' + sequence[loc] + '}}}$' + sequence[loc + 1:]))
    #         cur_axs.axhline(y=0, color='k', linewidth=.5)
    #         cur_axs.axvline(x=0, color='k', linewidth=.5)
    #     plt.savefig(f'{output_file}/{seq_location}_{j}')
    #     print(f'{output_file}/pngs/{seq_location}_{j} png was saved')

    fig2, axes = plt.subplots(1, 8, figsize=(50, 1))
    i = 0
    for pos in positions:
        ax = sns.kdeplot(data=df.loc[df[1] == pos], x='current diff', hue="Modification", ax=axes[i], legend=False)
        ax.set(xlabel=None)
        ax.set_title(pos)
        ax.set(ylabel=None)
        ax.set_xlim(-15, 10)
        i += 1


def kmer_plot_events_and_nnn_count(args: dict):
    def input_tsv(bc):
        import os
        file = f'{args.output}/aligned_reads_barcode{bc}_{args.output.split("/")[-2]}_{args.strand}/{args.location}_events_count_norm.tsv'
        if not os.path.exists(file):
            print(f"tsv file does not exists: {file}")
            return None
        return file

    if None not in [args.output]:
        # plotting the points

        barcode_dict = {"09": {"i": 0, "color": 'tab:blue', "name": "Cytosine"},
                        "10": {"i": 1, "color": 'tab:blue', "name": "5hmC"},
                        "11": {"i": 2, "color": 'tab:blue', "name": "5gmc"},
                        "12": {"i": 3, "color": 'tab:blue', "name": "5gmc-N3"}}

        arr = [-i for i in range(args.back_window - 1, -10, -1)]

        fig = plt.figure(figsize=(15, 10))
        ax1 = fig.add_subplot(111)
        sequence = args.sequence[:5] + 'm' + args.sequence[5 + 1:]
        ax1.set_xlabel(('                ' + '         '.join(sequence) + '\nModification offset'))
        ax1.set_ylabel('events')
        lns = [[None for b in barcode_dict] for i in [1, 2]]
        for barcode in args.barcodes:
            color = barcode_dict[barcode]["color"]
            pd_real = pd.read_csv(input_tsv(barcode), sep=" ", header=0)
            # axis values
            x = pd_real["mod_dist"]
            z = pd_real["event_norm"]
            lns[0][barcode_dict[barcode]["i"]] = ax1.plot(x, z, label=f'{barcode_dict[barcode]["name"]}_event_norm')
        ax1.tick_params(axis='y')

        # Adding Twin Axes to plot using dataset_2
        ax2 = ax1.twinx()

        ax2.set_ylabel(
            'NNN/events                                                                                                  ')
        for barcode in args.barcodes:
            color = barcode_dict[barcode]["color"]
            pd_real = pd.read_csv(input_tsv(barcode), sep=" ", header=0)
            # axis values
            x = pd_real["mod_dist"]
            y = pd_real["NNN_ratio"]
            lns[1][barcode_dict[barcode]["i"]] = ax2.plot(x, y,
                                                          label=f'{barcode_dict[barcode]["name"]}_NNN_ratio',
                                                          linestyle="dotted")
        ax2.tick_params(axis='y')

        plt.xticks(arr)
        ax2.set_ylim(0, 2)
        ax1.set_ylim(-0.05, 0.15)

        # Naming the x-axis, y-axis and the whole graph
        plt.title(f"Sequence\n{args.location}\n{args.sequence}")

        # Adding legend, which helps us recognize the curve according to it's color
        # added these three lines
        lnss = None
        for i in [0, 1]:
            for b in args.barcodes:
                if lnss == None:
                    lnss = lns[i][barcode_dict[b]["i"]]
                else:
                    lnss += lns[i][barcode_dict[b]["i"]]

        labs = [ln.get_label() for ln in lnss]
        plt.legend(lnss, labs, loc=0)

        # function to show the plot
        plt.savefig(f'{args.output}/events & NNN plots/{args.back_window}_back_{args.location}')
        print(f'{args.output}/events & NNN plots/{args.back_window}_back_{args.location}')

    else:
        raise argparse.ArgumentTypeError('The required arguments are: input and output')

    print("DONE!")


def arguments_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output')
    parser.add_argument('-b', '--barcodes')
    parser.add_argument('-i', '--input_cpg_cites')
    parser.add_argument('-s', '--strand')
    parser.add_argument('-seq', '--sequence')
    parser.add_argument('-w', '--back_window')
    return parser.parse_args()


if __name__ == "__main__":
    args = arguments_parser()
    args.barcodes = args.barcodes.split(' ')
    args.back_window = int(args.back_window)
    if not os.path.exists(f'{args.output}/events & NNN plots'):
        os.mkdir(f'{args.output}/events & NNN plots')

    with open(args.input_cpg_cites) as cites:
        for line in cites:
            line = line.split()
            args.location = "_".join(str(n) for n in [int(line[1])-args.back_window, int(line[2])+8])
            print(args.location)
            args.sequence = line[3]
            kmer_plot_events_and_nnn_count(args=args)
     # scatterplot_current_time_by_base(args=args, back=1, front=0)

