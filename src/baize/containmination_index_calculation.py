import logging
import os
import pickle
import re
import shutil
import sys
import warnings

import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from pandas import ExcelWriter

matplotlib.use('Agg')


def setup_logging(log_file):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    fonttools_logger = logging.getLogger("fontTools")
    fonttools_logger.propagate = False

    fonttools_logger.setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)  # Python Imaging Library

    fonttools_logger.addHandler(logging.NullHandler())

    return logger


def read_file(file_path):
    file_extension = file_path.split('.')[-1].lower()
    if file_extension == 'xlsx':
        return pd.read_excel(file_path)
    elif file_extension in ['csv', 'tsv', 'txt']:
        return pd.read_csv(file_path, sep='\t' if file_extension in ['tsv', 'txt'] else ',')
    else:
        raise ValueError(f"Unsupported file type: {file_extension}")


def contamination_index(marker, pro, reverse=False, model=1, n=1):
    pattern = '|'.join(re.escape(m) for m in marker)
    marker_pro = pro[pro.index.str.contains(pattern)]

    marker_exp = marker_pro.sum(axis=0)
    total_exp = pro.sum(axis=0)
    if (reverse == True):
        res = pd.DataFrame(total_exp / marker_exp)
    else:
        res = pd.DataFrame(marker_exp / total_exp)

    res.columns = ['contamination_index']
    res['Sample_ID'] = res.index
    res = res[['Sample_ID', 'contamination_index']].sort_values(by='Sample_ID')

    if (model != 1):
        res['contaminated'] = res['contamination_index'] > n

        res['log10_predicted_cell_count_per_ul_plasma_sample'] = res['contamination_index'].apply(model)

        res['predicted_cell_count_per_ul_plasma_sample'] = round(
            10 ** (res['log10_predicted_cell_count_per_ul_plasma_sample']) / 15)
        res['log10_predicted_cell_count_per_ul_plasma_sample'] = np.log10(
            res['predicted_cell_count_per_ul_plasma_sample'])
        mask = res['contaminated'] == False

        res.loc[mask, res.columns[-2:]] = None
    return res


def plot_dataframe(df, outfile_prefix, title='', width=1, color="blue"):
    fig, ax = plt.subplots(figsize=(width * 6, 5))

    ax = sns.barplot(x=df.index, y='contamination_index', data=df.reset_index(), color=color)

    mean_value = df['contamination_index'].mean()
    mean_line = ax.axhline(y=mean_value, color='red', linestyle='--',
                           label=f'Mean contamination index: {mean_value:.5f}')

    has_text = False
    handles, labels = [], []

    if 'log10_predicted_cell_count_per_ul_plasma_sample' in df.columns:
        for i, v in enumerate(df['log10_predicted_cell_count_per_ul_plasma_sample']):
            if pd.notna(v):
                has_text = True
                ax.text(i, df['contamination_index'][i], f'{v:.2f}', color='green', ha='center', va='bottom')
    if has_text:
        from matplotlib.lines import Line2D
        text_legend = Line2D([0], [0],
                             marker='$\mathsf{T}$',
                             markersize=12,
                             markerfacecolor='none',
                             markeredgecolor='green',
                             markeredgewidth=2,
                             linestyle='None',
                             label='log(predicted cell count per μl plasma sample)')
        handles = [mean_line, text_legend]
        labels = [line.get_label() for line in handles]
    ax.legend(handles=handles, labels=labels, loc='best', frameon=True)

    ax.set_title(title)
    ax.set_xlabel('')
    ax.set_ylabel('Contamination Index')

    plt.xticks(rotation=90)
    plt.tight_layout()

    return fig


def is_numeric_column(series):
    if pd.api.types.is_numeric_dtype(series):
        return True

    try:
        converted = pd.to_numeric(series, errors='coerce')
        return converted.notna().all()
    except:
        return False


def main(base_save_dir, pg_matrix):
    font = FontProperties(family='Arial')
    plt.rcParams['font.family'] = font.get_name()
    mpl.rcParams['pdf.fonttype'] = 42

    warnings.filterwarnings("ignore", category=UserWarning)

    # file_name = os.path.splitext(os.path.basename(pg_matrix))[0]
    pg_matrix = read_file(pg_matrix)

    outfile_prefix = os.path.join(base_save_dir, 'result')
    if os.path.exists(outfile_prefix):
        shutil.rmtree(outfile_prefix)

    os.makedirs(outfile_prefix, exist_ok=True)

    # os.makedirs(f"{file_name}_OmniProt_contamination_Calculation",exist_ok=True)
    log_path = os.path.join(base_save_dir, 'info.log')
    logger = setup_logging(log_path)

    platelet_marker = ["P48059", "P04899", "Q15286", "P06576", "Q93084", "P05556", "P48735", "P40197", "P24557",
                       "P05141", "P30048", "P14625", "P23229", "O60610", "O43182", "P10809", "P10720", "Q8TC12",
                       "Q9NR12", "P61106", "Q86YW5", "Q96AP7", "Q14165", "Q9BX67", "Q9Y277", "P40926", "P45880",
                       "P16615", "P62873", "P21796"]
    rbc_marker = ["O75955", "P00558", "P00915", "P02042", "P02549", "P02730", "P04406", "P11142", "P11171", "P11277",
                  "P13987", "P16157", "P16452", "P23528", "P27105", "P30043", "P35611", "P35612", "P37840", "P48426",
                  "P53396", "P60891", "P61224", "P62826", "P62937", "P68871", "P69905", "Q00013", "Q9NTK5", "Q9UQ80"]
    coagulation_marker = ['A6NMZ7', 'O95399', 'P02462', 'P02671', 'P02675', 'P02679', 'P02786', 'P12259', 'Q08397',
                          'Q08830', 'Q6UXI9', 'Q6ZMP0', 'Q8N302', 'Q8NI99', 'Q9BWP8', 'Q9UPU3', 'Q9Y546', 'Q9Y6L7',
                          'Q9Y6Y1', 'Q9Y6Z7']

    merged_markers = list(set(platelet_marker + rbc_marker + coagulation_marker))
    lower_merged_markers = {marker.lower() for marker in merged_markers}
    matching_columns = []

    for col in pg_matrix.columns:
        col_values = pg_matrix[col].dropna().astype(str).str.lower().unique()
        if any(marker in col_values for marker in lower_merged_markers):
            matching_columns.append(col)

    if len(matching_columns) == 0:
        platelet_marker = ["LIMS1", "GNAI2", "RAB35", "ATP5F1B", "ATP2A3", "ITGB1", "IDH2", "GP5", "TBXAS1", "SLC25A5",
                           "PRDX3", "HSP90B1", "ITGA6", "DIAPH1", "ARHGAP6", "HSPD1", "PF4V1", "RDH11", "PDLIM7",
                           "RAB14", "TREML1", "ESAM", "MLEC", "JAM3", "VDAC3", "MDH2", "VDAC2", "ATP2A2", "GNB1",
                           "VDAC1"]
        rbc_marker = ["FLOT1", "PGK1", "CA1", "HBD", "SPTA1", "SLC4A1", "GAPDH", "HSPA8", "EPB41", "SPTB", "CD59",
                      "ANK1", "EPB42", "CFL1", "STOM", "BLVRB", "ADD1", "ADD2", "SNCA", "PIP4K2A", "ACLY", "PRPS1",
                      "RAP1B", "RAN", "PPIA", "HBB", "HBA1", "MPP1", "OLA1", "PA2G4"]
        coagulation_marker = ['COL6A6', 'UTS2', 'COL4A1', 'FGA', 'FGB', 'FGG', 'TFRC', 'F5', 'LOXL1', 'FGL1', 'NPNT',
                              'THSD4', 'AGGF1', 'ANGPTL6', 'COLEC11', 'SORCS3', 'LRRC42', 'TLL2', 'CAMTA1', 'COLEC10']

        merged_markers = list(set(platelet_marker + rbc_marker + coagulation_marker))
        lower_merged_markers = {marker.lower() for marker in merged_markers}
        matching_columns = []

        for col in pg_matrix.columns:
            col_values = pg_matrix[col].dropna().astype(str).str.lower().unique()
            if any(marker in col_values for marker in lower_merged_markers):
                matching_columns.append(col)

        if len(matching_columns) == 0:
            logger.info(
                "No Uniprot ID or Gene ID column detected. Using the first column as the Uniprot ID column; please check whether the first column of the input file contains Uniprot IDs!")
            platelet_marker = ["P48059", "P04899", "Q15286", "P06576", "Q93084", "P05556", "P48735", "P40197", "P24557",
                               "P05141", "P30048", "P14625", "P23229", "O60610", "O43182", "P10809", "P10720", "Q8TC12",
                               "Q9NR12", "P61106", "Q86YW5", "Q96AP7", "Q14165", "Q9BX67", "Q9Y277", "P40926", "P45880",
                               "P16615", "P62873", "P21796"]
            rbc_marker = ["O75955", "P00558", "P00915", "P02042", "P02549", "P02730", "P04406", "P11142", "P11171",
                          "P11277", "P13987", "P16157", "P16452", "P23528", "P27105", "P30043", "P35611", "P35612",
                          "P37840", "P48426", "P53396", "P60891", "P61224", "P62826", "P62937", "P68871", "P69905",
                          "Q00013", "Q9NTK5", "Q9UQ80"]
            coagulation_marker = ['A6NMZ7', 'O95399', 'P02462', 'P02671', 'P02675', 'P02679', 'P02786', 'P12259',
                                  'Q08397', 'Q08830', 'Q6UXI9', 'Q6ZMP0', 'Q8N302', 'Q8NI99', 'Q9BWP8', 'Q9UPU3',
                                  'Q9Y546', 'Q9Y6L7', 'Q9Y6Y1', 'Q9Y6Z7']
            col = pg_matrix.columns[0]
        else:
            logger.info(f"Detected column {matching_columns[0]} as Gene Name.")
            col = matching_columns[0]
    else:
        logger.info(f"Detected column {matching_columns[0]} as Uniprot ID.")
        col = matching_columns[0]

    with open('resource/baize/platelet_model.pkl', 'rb') as f:
        platelet_spline = pickle.load(f)
    with open('resource/baize/rbc_model.pkl', 'rb') as f:
        rbc_spline = pickle.load(f)

    pg_matrix[col] = pg_matrix[col].str.upper()
    pg_matrix.index = pg_matrix[col]
    pg_matrix.index.name = None
    pg_matrix = pg_matrix[~pg_matrix[col].str.contains(';')]

    numeric_cols = [col for col in pg_matrix.columns if is_numeric_column(pg_matrix[col])]
    pg_matrix = pg_matrix.loc[:, numeric_cols]
    logger.info(f"Detected sample columns: {numeric_cols}")

    missing_rate = pg_matrix.iloc[:, :].isnull().mean(axis=1)
    pg_matrix = pg_matrix[missing_rate < 0.5]

    png_name_list = ['Platelet_contamination_index_profile.png', 'Erythrocyte_contamination_index_profile.png',
                     'Coagulation_contamination_index_profile.png']
    with PdfPages(os.path.join(outfile_prefix, 'OmniProt-based Plasma Contamination Index Profile.pdf')) as pdf:

        analyses = [
            ("platelet", platelet_marker, platelet_spline, 0.009, "#2d56a6", "Platelet contamination index"),
            ("rbc", rbc_marker, rbc_spline, 0.1, "#8e2122", "Red blood cell contamination index"),
            ("coagulation", coagulation_marker, 1, None, "#e19273", "Coagulation contamination index", True)
        ]

        results = {}

        for dd_index, analysis in enumerate(analyses):
            name, markers, model, n, color, title = analysis[:6]
            reverse = len(analysis) > 6 and analysis[6]

            df = contamination_index(
                markers,
                pg_matrix,
                model=model,
                n=n,
                reverse=reverse
            )

            results[name] = df

            fig = plot_dataframe(
                df,
                None,
                title=title,
                width=max(1, len(numeric_cols) / 20),
                color=color
            )

            pdf.savefig(fig)
            fig.savefig(os.path.join(outfile_prefix, png_name_list[dd_index]))
            plt.close(fig)

        excel_path = os.path.join(outfile_prefix, 'OmniProt-based Plasma Contamination Index Profile.xlsx')
        with ExcelWriter(excel_path) as writer:
            for sheet, df in results.items():
                df.to_excel(writer, sheet_name=sheet.capitalize(), index=False)
