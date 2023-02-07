from Bio import SeqIO
import pandas as pd
from os import path
import csv
import random
from statistics import mean
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
from sklearn.model_selection import cross_val_score
from collections import Counter
from imblearn.combine import SMOTEENN, SMOTETomek
import pickle
import argparse
import tqdm


CONSOLIDATED_EVENTS_FILE_HEADER = [
    "position",
    "model_kmer",
    "read_name",
    "events_number0",
    "events_number-1",
    "events_number-2",
    "events_number-3",
    "events_number-4",
    "events_number-5",
    "six_n_number0",
    "six_n_number-1",
    "six_n_number-2",
    "six_n_number-3",
    "six_n_number-4",
    "six_n_number-5",
    "last_kmer_six_n0",
    "last_kmer_six_n-1",
    "last_kmer_six_n-2",
    "last_kmer_six_n-3",
    "last_kmer_six_n-4",
    "last_kmer_six_n-5",
    "shift_from_model0",
    "shift_from_model-1",
    "shift_from_model-2",
    "shift_from_model-3",
    "shift_from_model-4",
    "shift_from_model-5",
    "total_event_length0",
    "total_event_length-1",
    "total_event_length-2",
    "total_event_length-3",
    "total_event_length-4",
    "total_event_length-5",
    "kmer_skipped0",
    "kmer_skipped-1",
    "kmer_skipped-2",
    "kmer_skipped-3",
    "kmer_skipped-4",
    "kmer_skipped-5",
    "sum_six_n_number_exists",
    "sum_last_kmer_six_n",
    "sum_abs_shift_from_model",
    "sum_skipped",
    "type_score"
]

FEATURES = [
    "events_number0",
    "events_number-1",
    "events_number-2",
    "events_number-3",
    "events_number-4",
    "events_number-5",
    "six_n_number0",
    "six_n_number-1",
    "six_n_number-2",
    "six_n_number-3",
    "six_n_number-4",
    "six_n_number-5",
    "last_kmer_six_n0",
    "last_kmer_six_n-1",
    "last_kmer_six_n-2",
    "last_kmer_six_n-3",
    "last_kmer_six_n-4",
    "last_kmer_six_n-5",
    "shift_from_model0",
    "shift_from_model-1",
    "shift_from_model-2",
    "shift_from_model-3",
    "shift_from_model-4",
    "shift_from_model-5",
    "total_event_length0",
    "total_event_length-1",
    "total_event_length-2",
    "total_event_length-3",
    "total_event_length-4",
    "total_event_length-5",
    "kmer_skipped0",
    "kmer_skipped-1",
    "kmer_skipped-2",
    "kmer_skipped-3",
    "kmer_skipped-4",
    "kmer_skipped-5",
    "sum_six_n_number_exists",
    "sum_last_kmer_six_n",
    "sum_abs_shift_from_model",
    "sum_skipped",
]

ACTION_FEATURES = 'extract-features'
ACTION_DETECT = 'detect'
ACTION_FEATURES_DETECT = 'extract-features-and-detect'
ACTION_TRAIN = 'train'

# consts
KMER_LENGTH = 6
CG_SEQ = "CG"
SIX_N_KMER = "NNNNNN"

# Files suffixes
OUTPUT_FILE_SUFFIX = ".tsv"
BEDGRAPH_FILE_SUFFIX = ".bedgraph"

# The values of the positions in the "eventalign" file
POSITION = 1
READ_NAME = 3
EVENT_LEVEL_MEAN = 6
EVENT_LENGTH = 8
MODEL_KMER = 9
MODEL_MEAN = 10
EVENT_STANDARDIZED_LEVEL = 12

# Data frame column names
KMER_POSITION = "position"
KMER_BASE = "base"
KMER_SEQ = "model_kmer"
KMER_READ_NAME = "read_name"
KMER_EVENTS_NUMBER = "events_number"
KMER_SIX_N_NUMBER = "six_n_number"
KMER_TOTAL_SHIFT_FROM_MODEL = "shift_from_model"
KMER_TOTAL_EVENT_LENGTH = "total_event_length"
KMER_SIX_N_TOTAL_EVENT_LENGTH = "total_event_length_six_n"
KMER_WITH_LAST_NNNNN = "last_six_n"
KMER_SKKIPED = "skipped_kmer"
SKIPPED_EVENTS_RATIO = "skipped_events_ratio"
SIX_N_RATIO = "six_n_events_ratio"
SIX_N_IN_LAST_EVENT_RATIO = "six_n_last_event_ratio"
SUM_ABS_TOTAL_SHIFT_FROM_MODEL = "abs_shift_from_model"
DROP_ME = "drop_me"
SUM_FEATURES = "sum_features"
MOD_RATIO = "modification_ratio"
MOD_COUNT = "modification_ratio"
COVERAGE = "coverage"
CHR = "chr"
POSITION_START = "position_start"
POSITION_END = "position_end"
ANALYSIS_RESULT = "analysis_result"
TYPE_SCORE = "type_score"

DF_BASIC_COLUMNS = [
    KMER_READ_NAME,
    KMER_EVENTS_NUMBER,
    KMER_SIX_N_NUMBER,
    KMER_TOTAL_SHIFT_FROM_MODEL,
    KMER_TOTAL_EVENT_LENGTH,
    KMER_SIX_N_TOTAL_EVENT_LENGTH,
    KMER_WITH_LAST_NNNNN,
    KMER_SKKIPED,
    TYPE_SCORE,
    DROP_ME
]

READS_NUMBER = "reads_number"
ABS_TOTAL_SHIFT_FROM_MODEL = "abs_shift_from_model"

RESULTS_COLUMNS = [
    KMER_POSITION,
    SKIPPED_EVENTS_RATIO,
    SIX_N_RATIO,
    SIX_N_IN_LAST_EVENT_RATIO,
    READS_NUMBER,
    ABS_TOTAL_SHIFT_FROM_MODEL
]

KMER_ANALYSIS_RESULT = "analysis_result"

DF_COLUMNS = [
    KMER_READ_NAME,
    KMER_ANALYSIS_RESULT,
    DROP_ME
]

CONSOLIDATED_EVENTS_F_OUTPUT_FILE = "_consolidated_events_forward_strand.tsv"
CONSOLIDATED_EVENTS_R_OUTPUT_FILE = "_consolidated_events_reverse_strand.tsv"
CONSOLIDATED_EVENTS_F_OUTPUT_FILE_RELEVANT = "_consolidated_events_forward_strand_relevant.tsv"
CONSOLIDATED_EVENTS_R_OUTPUT_FILE_RELEVANT = "_consolidated_events_reverse_strand_relevant.tsv"


def arguments_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=[ACTION_FEATURES, ACTION_DETECT, ACTION_FEATURES_DETECT, ACTION_TRAIN],
                        help="Choose the desired action to perform")
    # EXTRACT FEATURES & DETECT
    parser.add_argument('-i', '--input', nargs="+",
                        help="The path to the 'eventalign' result file")
    parser.add_argument('--ref', dest='ref_fasta',
                        help="The path to reference genome in a FASTA file format")
    parser.add_argument('-f', '--full-seq', dest='is_full_seq', action='store_true',
                        help="Should be added if the FASTA file includes the full sequence")
    parser.add_argument('-sp', '--start-position', type=int)
    parser.add_argument('-ep', '--end-position', type=int)
    parser.add_argument('-s', '--strand', choices=['+', '-'])
    parser.add_argument('-rn', '--reads-number',
                        help="Should be used in order to perform the analysis with a specific reads number, "
                             "the reads will be chosen randomly")
    parser.add_argument('-ml', '--min-length', type=int,
                        help="Should be used in order to perform the analysis with reads that are longer then this value")
    parser.add_argument('-o', '--output')
    # DETECT
    parser.add_argument('--input-pkl',
                        help="The path to the 'pickle' file that will be used for the detection")
    # TRAIN
    parser.add_argument("--output-pkl",
                        help="The path to the 'pickle' file that will be generate after training a new model")

    return parser.parse_args()


def extract_features(input_file, ref_fasta, is_full_seq, start_position, end_position, strand, reads_number, output):
    # Save the parsed arguments to relevant variables
    ref_seq = SeqIO.read(ref_fasta, "fasta")
    chr_name = (path.basename(ref_fasta)).split(".")[0]
    seq_length = len(ref_seq)

    # If the fasta file contains the full sequence of the chromosome/genome
    # Cut out just the relevant region
    if is_full_seq:
        ref_seq = ref_seq.seq[start_position :end_position + 1]
        seq_length = len(ref_seq)
    else:
        ref_seq = ref_seq.seq[:]
        # If the fasta file contains just the relevant sequence
        # Check if the start and the end positions are match the length of this sequence
        if seq_length != end_position - start_position + 1:
            print(seq_length)
            raise Exception(
                "The length of the sequence does not match the length obtained between start and end position (try '-f’")

    current_reads_list = []
    with open(input_file) as eventalign_file:
        eventalign_tsv_file = pd.read_csv(eventalign_file, delimiter="\t", usecols=[READ_NAME])
        current_reads_list = eventalign_tsv_file.iloc[:, 0].unique()
        eventalign_tsv_file_len = len(eventalign_tsv_file)

    if reads_number is not None:
        print("[*] Choose " + str(reads_number) + " random reads to work with")
        while len(current_reads_list) > args.reads_number:
            current_reads_list.remove(random.choice(current_reads_list))
    else:
        reads_number = len(current_reads_list)
        print("[*] current_reads_list length: " + str(reads_number))

    print("[1] Open the 'eventalign' input file: \n" + str(input_file))
    with open(input_file) as eventalign_file:
        eventalign_tsv_file = csv.reader(eventalign_file, delimiter="\t")
        consolidated_events_file_writer = None
        consolidated_events_file_path = None

        if strand == "+":
            consolidated_events_file_path = output + CONSOLIDATED_EVENTS_F_OUTPUT_FILE
            consolidated_events_f_tsv_file = open(consolidated_events_file_path, 'w', newline='')
            consolidated_events_file_writer = csv.writer(consolidated_events_f_tsv_file, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            print("[*] Start writing the consolidated events from the forward strand to file: \n" + consolidated_events_file_path)

        if strand == "-":
            consolidated_events_file_path = output + CONSOLIDATED_EVENTS_R_OUTPUT_FILE
            consolidated_events_r_tsv_file = open(consolidated_events_file_path, 'w', newline='')
            consolidated_events_file_writer = csv.writer(consolidated_events_r_tsv_file, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            print("[*] Start writing the consolidated events from the reverse strand to file: \n" + consolidated_events_file_path)

        consolidated_events_file_writer.writerow(CONSOLIDATED_EVENTS_FILE_HEADER)

        # Create data frame to store the "per-kmer" information
        current_df = pd.DataFrame(columns=DF_BASIC_COLUMNS)

        header = next(eventalign_tsv_file)
        first_line = next(eventalign_tsv_file)
        while start_position >= int(first_line[POSITION]) or \
                int(first_line[POSITION]) >= end_position or \
                first_line[READ_NAME] not in current_reads_list:
            first_line = next(eventalign_tsv_file)

        drop_me = 0
        kmer_index = 0
        score = 0

        current_read = first_line[READ_NAME]
        current_position = int(first_line[POSITION])
        current_kmer = first_line[MODEL_KMER]
        total_event_length = 0
        total_event_length_six_n = 0
        total_shift_from_model = 0
        total_event_level_mean = 0
        sum_six_n_number = 0
        events_number = 1
        last_kmer_six_n = 0

        if current_kmer == SIX_N_KMER:
            sum_six_n_number = 1
            last_kmer_six_n = 1
            total_event_length_six_n = float(first_line[EVENT_LENGTH])
        else:
            total_event_length = float(first_line[EVENT_LENGTH])
            total_shift_from_model = float(first_line[EVENT_LENGTH]) * float(first_line[EVENT_STANDARDIZED_LEVEL])
            total_event_level_mean = float(first_line[EVENT_LENGTH]) * float(first_line[EVENT_LEVEL_MEAN])

        print("Extracting features")
        print("<< s t a r t >>")
        with tqdm.tqdm(smoothing=0.1, total=eventalign_tsv_file_len) as pbar:
            for current_event in eventalign_tsv_file:
                pbar.update()
                if current_event[READ_NAME] not in current_reads_list:
                    continue
                if start_position <= int(current_event[POSITION]) <= end_position:
                    # If this event is from the same position and read as the previous event (so it is the same k-mer),
                    # continue the current consolidation calculation
                    if current_event[READ_NAME] == current_read and int(current_event[POSITION]) == current_position:
                        if current_event[MODEL_KMER] == SIX_N_KMER:
                            sum_six_n_number += 1
                            total_event_length_six_n = total_event_length_six_n + float(current_event[EVENT_LENGTH])
                            last_kmer_six_n = 1
                        else:
                            total_event_length = total_event_length + float(current_event[EVENT_LENGTH])
                            total_shift_from_model = total_shift_from_model + (
                                        float(current_event[EVENT_LENGTH]) * float(current_event[EVENT_STANDARDIZED_LEVEL]))
                            total_event_level_mean = total_event_level_mean + (
                                    float(current_event[EVENT_LENGTH]) * float(current_event[EVENT_LEVEL_MEAN]))
                            last_kmer_six_n = 0

                            if current_kmer == SIX_N_KMER:
                                current_kmer = current_event[MODEL_KMER]

                        events_number += 1

                    # The 'current_event' is not from the same read as the previous one,
                    # or there were skipped events so the current position is not the same as the previous one
                    # Calc everything that is needed about the "current_kmer" and write it to the consolidation output file
                    # The 'current_event' will be part of the "current_kmer" in the next iteration
                    else:
                        shift_from_model = 0
                        if not total_event_length == 0:
                            shift_from_model = total_shift_from_model / total_event_length

                        # Add the information about this k-mer to the data frame
                        if (kmer_index > 5) and (current_df.at[kmer_index - 5, KMER_READ_NAME] == current_read):
                            drop_me = 0
                        else:
                            # If there isn't enough information about the current read,
                            # use this k-mer for future calculations, but don't write it to the output file
                            drop_me = 1

                        if abs(shift_from_model) > 1.5:
                            score = 1.2
                        elif sum_six_n_number > 0:
                            score = 1.3
                        else:
                            score = -1

                        # Add it to current data frame
                        current_df.loc[kmer_index] = [
                            current_read,
                            events_number,
                            sum_six_n_number,
                            shift_from_model,
                            total_event_length,
                            total_event_length_six_n,
                            last_kmer_six_n,
                            0,  # The kmer wasn't skipped
                            score,
                            drop_me
                        ]

                        if drop_me == 0:
                            # Add it to the output file
                            new_line = [
                                current_position,
                                current_kmer,
                                current_read,
                                events_number,
                                current_df.loc[kmer_index - 1][KMER_EVENTS_NUMBER],
                                current_df.loc[kmer_index - 2][KMER_EVENTS_NUMBER],
                                current_df.loc[kmer_index - 3][KMER_EVENTS_NUMBER],
                                current_df.loc[kmer_index - 4][KMER_EVENTS_NUMBER],
                                current_df.loc[kmer_index - 5][KMER_EVENTS_NUMBER],
                                sum_six_n_number,
                                current_df.loc[kmer_index - 1][KMER_SIX_N_NUMBER],
                                current_df.loc[kmer_index - 2][KMER_SIX_N_NUMBER],
                                current_df.loc[kmer_index - 3][KMER_SIX_N_NUMBER],
                                current_df.loc[kmer_index - 4][KMER_SIX_N_NUMBER],
                                current_df.loc[kmer_index - 5][KMER_SIX_N_NUMBER],
                                last_kmer_six_n,
                                current_df.loc[kmer_index - 1][KMER_WITH_LAST_NNNNN],
                                current_df.loc[kmer_index - 2][KMER_WITH_LAST_NNNNN],
                                current_df.loc[kmer_index - 3][KMER_WITH_LAST_NNNNN],
                                current_df.loc[kmer_index - 4][KMER_WITH_LAST_NNNNN],
                                current_df.loc[kmer_index - 5][KMER_WITH_LAST_NNNNN],
                                abs(shift_from_model),
                                abs(current_df.loc[kmer_index - 1][KMER_TOTAL_SHIFT_FROM_MODEL]),
                                abs(current_df.loc[kmer_index - 2][KMER_TOTAL_SHIFT_FROM_MODEL]),
                                abs(current_df.loc[kmer_index - 3][KMER_TOTAL_SHIFT_FROM_MODEL]),
                                abs(current_df.loc[kmer_index - 4][KMER_TOTAL_SHIFT_FROM_MODEL]),
                                abs(current_df.loc[kmer_index - 5][KMER_TOTAL_SHIFT_FROM_MODEL]),
                                total_event_length + total_event_length_six_n,
                                current_df.loc[kmer_index - 1][KMER_TOTAL_EVENT_LENGTH] +
                                current_df.loc[kmer_index - 1][KMER_SIX_N_TOTAL_EVENT_LENGTH],
                                current_df.loc[kmer_index - 2][KMER_TOTAL_EVENT_LENGTH] +
                                current_df.loc[kmer_index - 2][KMER_SIX_N_TOTAL_EVENT_LENGTH],
                                current_df.loc[kmer_index - 3][KMER_TOTAL_EVENT_LENGTH] +
                                current_df.loc[kmer_index - 3][KMER_SIX_N_TOTAL_EVENT_LENGTH],
                                current_df.loc[kmer_index - 4][KMER_TOTAL_EVENT_LENGTH] +
                                current_df.loc[kmer_index - 4][KMER_SIX_N_TOTAL_EVENT_LENGTH],
                                current_df.loc[kmer_index - 5][KMER_TOTAL_EVENT_LENGTH] +
                                current_df.loc[kmer_index - 5][KMER_SIX_N_TOTAL_EVENT_LENGTH],
                                0,  # The kmer wasn't skipped
                                current_df.loc[kmer_index - 1][KMER_SKKIPED],
                                current_df.loc[kmer_index - 2][KMER_SKKIPED],
                                current_df.loc[kmer_index - 3][KMER_SKKIPED],
                                current_df.loc[kmer_index - 4][KMER_SKKIPED],
                                current_df.loc[kmer_index - 5][KMER_SKKIPED],
                                sum(
                                    [int(bool(six_n_number)) for six_n_number in
                                     current_df.loc[kmer_index - 5:kmer_index][KMER_SIX_N_NUMBER]],
                                    int(bool(sum_six_n_number))
                                ),
                                sum(
                                    current_df.loc[kmer_index - 5:kmer_index][KMER_WITH_LAST_NNNNN],
                                    last_kmer_six_n
                                ),
                                sum(
                                    [abs(shift) for shift in
                                     current_df.loc[kmer_index - 5:kmer_index][KMER_TOTAL_SHIFT_FROM_MODEL]],
                                    abs(shift_from_model)
                                ),
                                sum(current_df.loc[kmer_index - 5:kmer_index][KMER_SKKIPED]),
                                score
                            ]

                            # print(new_line)
                            consolidated_events_file_writer.writerow(new_line)

                        if kmer_index > 5:
                            current_df = current_df.drop(kmer_index - 6)

                        kmer_index += 1

                        # These are all the skipped k-mers
                        if current_event[READ_NAME] == current_read and int(current_event[POSITION]) != current_position:

                            # Count all the skipped events between these positions
                            for skipped_position in range(current_position + 1, int(current_event[POSITION])):

                                # Add the information about this k-mer to the data frame
                                if (kmer_index > 5) and (current_df.at[kmer_index - 5, KMER_READ_NAME] == current_read):
                                    drop_me = 0
                                else:
                                    # If there isn't enough information about the current read,
                                    # use this k-mer for future calculations, but don't write it to the output file
                                    drop_me = 1

                                score = 1.1

                                # Add it to current data frame
                                current_df.loc[kmer_index] = [
                                    current_read,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    1,  # The kmer wasn't skipped
                                    score,
                                    drop_me
                                ]

                                if drop_me == 0:
                                    # Add it to the output file
                                    new_line = [
                                        current_position,
                                        current_kmer,
                                        current_read,
                                        0,
                                        current_df.loc[kmer_index - 1][KMER_EVENTS_NUMBER],
                                        current_df.loc[kmer_index - 2][KMER_EVENTS_NUMBER],
                                        current_df.loc[kmer_index - 3][KMER_EVENTS_NUMBER],
                                        current_df.loc[kmer_index - 4][KMER_EVENTS_NUMBER],
                                        current_df.loc[kmer_index - 5][KMER_EVENTS_NUMBER],
                                        0,
                                        current_df.loc[kmer_index - 1][KMER_SIX_N_NUMBER],
                                        current_df.loc[kmer_index - 2][KMER_SIX_N_NUMBER],
                                        current_df.loc[kmer_index - 3][KMER_SIX_N_NUMBER],
                                        current_df.loc[kmer_index - 4][KMER_SIX_N_NUMBER],
                                        current_df.loc[kmer_index - 5][KMER_SIX_N_NUMBER],
                                        0,
                                        current_df.loc[kmer_index - 1][KMER_WITH_LAST_NNNNN],
                                        current_df.loc[kmer_index - 2][KMER_WITH_LAST_NNNNN],
                                        current_df.loc[kmer_index - 3][KMER_WITH_LAST_NNNNN],
                                        current_df.loc[kmer_index - 4][KMER_WITH_LAST_NNNNN],
                                        current_df.loc[kmer_index - 5][KMER_WITH_LAST_NNNNN],
                                        0,
                                        abs(current_df.loc[kmer_index - 1][KMER_TOTAL_SHIFT_FROM_MODEL]),
                                        abs(current_df.loc[kmer_index - 2][KMER_TOTAL_SHIFT_FROM_MODEL]),
                                        abs(current_df.loc[kmer_index - 3][KMER_TOTAL_SHIFT_FROM_MODEL]),
                                        abs(current_df.loc[kmer_index - 4][KMER_TOTAL_SHIFT_FROM_MODEL]),
                                        abs(current_df.loc[kmer_index - 5][KMER_TOTAL_SHIFT_FROM_MODEL]),
                                        0,
                                        current_df.loc[kmer_index - 1][KMER_TOTAL_EVENT_LENGTH] +
                                        current_df.loc[kmer_index - 1][KMER_SIX_N_TOTAL_EVENT_LENGTH],
                                        current_df.loc[kmer_index - 2][KMER_TOTAL_EVENT_LENGTH] +
                                        current_df.loc[kmer_index - 2][KMER_SIX_N_TOTAL_EVENT_LENGTH],
                                        current_df.loc[kmer_index - 3][KMER_TOTAL_EVENT_LENGTH] +
                                        current_df.loc[kmer_index - 3][KMER_SIX_N_TOTAL_EVENT_LENGTH],
                                        current_df.loc[kmer_index - 4][KMER_TOTAL_EVENT_LENGTH] +
                                        current_df.loc[kmer_index - 4][KMER_SIX_N_TOTAL_EVENT_LENGTH],
                                        current_df.loc[kmer_index - 5][KMER_TOTAL_EVENT_LENGTH] +
                                        current_df.loc[kmer_index - 5][KMER_SIX_N_TOTAL_EVENT_LENGTH],
                                        1,  # The kmer wasn't skipped
                                        current_df.loc[kmer_index - 1][KMER_SKKIPED],
                                        current_df.loc[kmer_index - 2][KMER_SKKIPED],
                                        current_df.loc[kmer_index - 3][KMER_SKKIPED],
                                        current_df.loc[kmer_index - 4][KMER_SKKIPED],
                                        current_df.loc[kmer_index - 5][KMER_SKKIPED],
                                        sum(
                                            [int(bool(six_n_number)) for six_n_number in
                                             current_df.loc[kmer_index - 5:kmer_index][KMER_SIX_N_NUMBER]],
                                            int(bool(sum_six_n_number))
                                        ),
                                        sum(
                                            current_df.loc[kmer_index - 5:kmer_index][KMER_WITH_LAST_NNNNN],
                                            last_kmer_six_n
                                        ),
                                        sum(
                                            [abs(shift) for shift in
                                             current_df.loc[kmer_index - 5:kmer_index][KMER_TOTAL_SHIFT_FROM_MODEL]],
                                            abs(shift_from_model)
                                        ),
                                        sum(current_df.loc[kmer_index - 5:kmer_index][KMER_SKKIPED]),
                                        score
                                    ]

                                    consolidated_events_file_writer.writerow(new_line)

                                if kmer_index > 5:
                                    current_df = current_df.drop(kmer_index - 6)

                                kmer_index += 1

                        current_position = int(current_event[POSITION])
                        current_read = current_event[READ_NAME]
                        current_kmer = current_event[MODEL_KMER]
                        events_number = 1

                        if current_kmer == SIX_N_KMER:
                            total_event_length = 0
                            total_event_level_mean = 0
                            sum_six_n_number = 1
                            total_event_length_six_n = float(current_event[EVENT_LENGTH])
                            last_kmer_six_n = 1
                        else:
                            total_event_length = float(current_event[EVENT_LENGTH])
                            total_shift_from_model = (
                                        float(current_event[EVENT_LENGTH]) * float(current_event[EVENT_STANDARDIZED_LEVEL]))
                            total_event_level_mean = (
                                        float(current_event[EVENT_LENGTH]) * float(current_event[EVENT_LEVEL_MEAN]))
                            sum_six_n_number = 0
                            total_event_length_six_n = 0
                            last_kmer_six_n = 0

    print("Reads number: " + str(reads_number))
    return consolidated_events_file_path


def detect(input_file, ref_fasta, is_full_seq, start_position, end_position, strand, reads_number, min_length, output, input_pkl):
    ref_seq = SeqIO.read(ref_fasta, "fasta")
    with open(ref_fasta) as fasta:
        chr_name = fasta.readline().split(" ")[0]

    # If the fasta file contains the full sequence of the chromosome/genome
    # Cut out just the relevant sequence
    if is_full_seq:
        ref_seq = ref_seq.seq[start_position:end_position + 1]
    else:
        ref_seq = ref_seq.seq[:]
    print("ref len = " + str(len(ref_seq)))

    print("Load the input file")
    pd_real = pd.read_csv(input_file, sep="\t")

    reads_names = list(set(pd_real["read_name"]))
    current_reads_number = len(reads_names)
    print("initial reads number: " + str(current_reads_number))
    relevant_reads_names = []
    if min_length is not None:
        for read_name in reads_names:
            if len(pd_real[pd_real["read_name"] == read_name]) >= int(min_length):
                relevant_reads_names.append(read_name)
    else:
        relevant_reads_names = reads_names

    if reads_number is not None:
        print("relevant reads number: " + str(len(relevant_reads_names)))
        current_reads_number = int(reads_number)
        if current_reads_number > len(relevant_reads_names):
            raise ("Pick a higher reads number to work with")
        print("[*] Choose " + str(current_reads_number) + " random reads to work with")
        while len(relevant_reads_names) > current_reads_number:
            relevant_reads_names.remove(random.choice(relevant_reads_names))
        pd_real = pd_real[pd_real["read_name"].isin(relevant_reads_names) == True]

    # Loading the saved random forest model (pickle)
    random_forest_model_pkl = open(input_pkl, 'rb')
    random_forest_model = pickle.load(random_forest_model_pkl)
    print("Loaded Random Forest model :: " + str(random_forest_model))

    X_real = pd_real[FEATURES].values  # Features
    y_pred = random_forest_model.predict(X_real)
    pd_real["label"] = y_pred

    pd_results = pd_real[["read_name", "position", "label"]]
    pd_results.to_csv(output + "_" + str(current_reads_number) + "_reads.tsv", index=False, sep="\t")

    print("Create bedgraph file")
    results_df_columns = [CHR, POSITION_START, POSITION_END, MOD_RATIO, MOD_COUNT, COVERAGE]
    results_df_ratio = pd.DataFrame(columns=results_df_columns)

    for current_position in range(start_position, end_position + 1):
        current_df = pd_results[pd_results["position"] == current_position][pd_results["label"] == 1]
        reads_count = pd_results["position"].shape[0]
        if current_df.empty:
            modification_count = -1
        else:
            modification_count = current_df.shape[0]  # reads count for current_position with ["label"] == 1

        results_df_ratio.loc[current_position] = [
            chr_name,
            current_position,
            current_position + 1,
            modification_count/reads_count,
            modification_count,
            reads_count
        ]

    results_df_ratio.to_csv(
        output + "_" + str(current_reads_number) + "_reads.bedgraph",
        index=False,
        header=False,
        sep='\t'
    )

    print("Create only CG output files")

    cg_1_results_df_ratio = results_df_ratio.copy()
    cg_3_results_df_ratio = results_df_ratio.copy()
    cg_5_results_df_ratio = results_df_ratio.copy()

    index = 0
    cg_index_list = []
    ref_seq_upper = ref_seq.upper()

    if args.strand == "+":
        while index != -1:
            index = ref_seq_upper.find("CG", index + 1)
            cg_index_list.append(start_position + index)

    if args.strand == "-":
        while index != -1:
            index = ref_seq_upper.find("CG", index + 1)
            cg_index_list.append(start_position + index + 1)

    for current_position in range(start_position, end_position + 1):
        if current_position not in range(start_position + 5, end_position - 1):
            cg_1_results_df_ratio[MOD_RATIO][current_position] = 0
            cg_3_results_df_ratio[MOD_RATIO][current_position] = 0
            cg_5_results_df_ratio[MOD_RATIO][current_position] = 0
        elif current_position not in cg_index_list:
            cg_1_results_df_ratio[MOD_RATIO][current_position] = 0
            if not ((current_position + 1 in cg_index_list) or
                    (current_position - 1 in cg_index_list)):
                cg_3_results_df_ratio[MOD_RATIO][current_position] = 0
                if not ((current_position + 2 in cg_index_list) or
                        (current_position - 2 in cg_index_list)):
                    cg_5_results_df_ratio[MOD_RATIO][current_position] = 0

    cg_1_results_df_ratio.to_csv(
        output + "_cg_1_" + str(current_reads_number) + "_reads.bedgraph",
        index=False,
        header=False,
        sep='\t'
    )

    cg_3_results_df_ratio.to_csv(
        output + "_cg_3_" + str(current_reads_number) + "_reads.bedgraph",
        index=False,
        header=False,
        sep='\t'
    )

    cg_5_results_df_ratio.to_csv(
        output + "_cg_5_" + str(current_reads_number) + "_reads.bedgraph",
        index=False,
        header=False,
        sep='\t'
    )


def train(inputs, output_pkl):
    print("Load dataset")
    input_files = []
    for input_file in inputs:
        input_files.append(pd.read_csv(input_file, sep="\t"))
    pd_train = pd.concat(input_files)
    feature_cols = FEATURES

    print("Set training set")
    X_train = pd_train[feature_cols]  # Features
    y_train = pd_train.label  # Target variable
    print(sorted(Counter(y_train).items()))

    smote_tomek = SMOTETomek(random_state=0)
    X_resampled, y_resampled = smote_tomek.fit_resample(X_train, y_train)
    print(sorted(Counter(y_resampled).items()))

    X_train = X_resampled
    y_train = y_resampled

    print("Create a RandomForestClassifier")
    clf = RandomForestClassifier(n_estimators=500)

    print("Train the model using the training sets")
    clf.fit(X_train, y_train)

    print("Predict the response for test set")
    y_pred = clf.predict(X_test)

    # Print the model Accuracy, how often is the classifier correct?
    print("Accuracy:", accuracy_score(y_test, y_pred))
    print("------------------")
    print(confusion_matrix(y_test, y_pred))
    print("------------------")
    print(classification_report(y_test, y_pred))

    print("Save the trained model")
    # Dump the random forest classifier with Pickle
    random_forest_model_pkl = open(output_pkl, 'wb')
    pickle.dump(clf, random_forest_model_pkl)
    # Close the pickle instances
    random_forest_model_pkl.close()


if __name__ == "__main__":
    args = arguments_parser()

    if args.action == ACTION_FEATURES:
        if None not in [args.input, args.ref_fasta, args.start_position, args.end_position, args.strand, args.output]:
            if len(args.input) == 1:
                extract_features(args.input[0], args.ref_fasta, args.is_full_seq, args.start_position, args.end_position,
                                 args.strand, args.reads_number, args.output)
        else:
            raise argparse.ArgumentTypeError('The required arguments are: '
                                             'input, ref-fasta, start-position, end-position, strand and output')

    elif args.action == ACTION_DETECT:
        if None not in [args.input, args.ref_fasta, args.start_position, args.end_position, args.strand, args.output, args.input_pkl]:
            if len(args.input) == 1:
                detect(args.input[0], args.ref_fasta, args.is_full_seq, args.start_position, args.end_position, args.strand,
                       args.reads_number, args.min_length, args.output, args.input_pkl)
        else:
            raise argparse.ArgumentTypeError('The required arguments are: '
                                             'input, ref-fasta, start-position, end-position, strand, input-pkl and output')

    elif args.action == ACTION_FEATURES_DETECT:
        if None not in [args.input, args.ref_fasta, args.start_position, args.end_position, args.strand, args.output, args.input_pkl]:
            if len(args.input) == 1:
                results = extract_features(args.input[0], args.ref_fasta, args.is_full_seq,
                                           args.start_position, args.end_position, args.strand, args.reads_number, args.output)
                detect(results, args.ref_fasta, args.is_full_seq, args.start_position, args.end_position, args.strand,
                       args.reads_number, args.min_length, args.output, args.input_pkl)
        else:
            raise argparse.ArgumentTypeError('The required arguments are: '
                                             'input, ref-fasta, start-position, end-position, strand, input-pkl and output')

    elif args.action == ACTION_TRAIN:
        if None not in [args.input, args.output_pkl]:
            train(args.input, args.output_pkl)
        else:
            raise argparse.ArgumentTypeError('The required arguments are: input and output-pkl')

    print("DONE!")
