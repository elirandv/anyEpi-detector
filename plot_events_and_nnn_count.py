# importing the required module
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def arguments_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-b', '--barcode')
    parser.add_argument('-l', '--location')
    parser.add_argument('-s', '--strand')
    parser.add_argument('-seq', '--sequence')
    parser.add_argument('-w', '--back_window')
    return parser.parse_args()


if __name__ == "__main__":
    args = arguments_parser()
    if None not in [args.input, args.output]:
        back_window = int(args.back_window)
        arr = [-i for i in range(back_window - 1, -10, -1)]
        arr2 = [i / 10 for i in range(11)]

        print("Load the input file")
        pd_real = pd.read_csv(args.input, sep=" ", header=0)
        # x axis values
        x = pd_real["mod_dist"]
        y = pd_real["NNN_ratio"]
        z = pd_real["event_norm"]
        # corresponding y axis values
        # plotting the points

        fig, ax1 = plt.subplots()

        color = 'tab:green'
        ax1.set_xlabel('Modification offset')
        ax1.set_ylabel('events', color=color)
        lns1 = ax1.plot(x, z, color=color, label='event_count')
        ax1.tick_params(axis='y', labelcolor=color)

        # Adding Twin Axes to plot using dataset_2
        ax2 = ax1.twinx()

        color = 'tab:red'
        ax2.set_ylabel('NNN/events', color=color)
        lns2 = ax2.plot(x, y, color=color, label='NNN_ratio')
        ax2.tick_params(axis='y', labelcolor=color)

        plt.xticks(arr)
        ax2.set_ylim(0, 2)
        ax1.set_ylim(0, pd_real["event_count"].max() + 2000)

        # Naming the x-axis, y-axis and the whole graph
        plt.title(f"Barcode{args.barcode}\n{args.strand}\n{args.location}\n{args.sequence}")

        # Adding legend, which helps us recognize the curve according to it's color
        # added these three lines
        lns = lns1 + lns2
        labs = [l.get_label() for l in lns]
        plt.legend(lns, labs, loc=0)

        # function to show the plot
        plt.savefig(args.output)

    else:
        raise argparse.ArgumentTypeError('The required arguments are: input and output')

    print("DONE!")
    exit()
