import os
import re
import shutil

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages  # 导入 PdfPages 类
from matplotlib.backends.backend_agg import FigureCanvasAgg    # 导入 PdfPages 类
from matplotlib.font_manager import FontProperties
import matplotlib as mpl

font = FontProperties(family='Arial')
plt.rcParams['font.family'] = font.get_name()
mpl.rcParams['pdf.fonttype'] = 42


def read_file(file_path):
    file_extension = file_path.split('.')[-1].lower()
    if file_extension == 'xlsx':
        return pd.read_excel(file_path)
    elif file_extension in ['csv', 'tsv', 'txt']:
        return pd.read_csv(file_path, sep='\t' if file_extension in ['tsv', 'txt'] else ',')
    else:
        raise ValueError(f"Unsupported file type: {file_extension}")


def contamination_index(marker,pro,reverse = False):
    pattern = '|'.join(re.escape(m) for m in marker)
    marker_pro = pro[pro.index.str.contains(pattern)]

    marker_exp = marker_pro.sum(axis=0)
    total_exp = pro.sum(axis=0)
    if(reverse ==True):
        res = pd.DataFrame(total_exp/marker_exp)
    else:
        res = pd.DataFrame(marker_exp/total_exp)

    res.columns = ['contamination_index']
    res['Sample_ID'] = res.index

    return res[['Sample_ID','contamination_index']].sort_values(by = 'Sample_ID')


def plot_dataframe(df, outfile_prefix, title='',width = 1,color = "blue"):
    ## 20个样本，需要宽度为6画出来好看，根绝样本量调节宽度
    fig, ax = plt.subplots(figsize=(width * 6, 5))

    # 绘制柱状图
    ax = sns.barplot(x=df.index, y='contamination_index', data=df.reset_index(),color = color)

    # 添加平均值和标准差的水平线
    mean_value = df['contamination_index'].mean()
    #std_dev = df['contamination_index'].std()
    ax.axhline(y=mean_value, color='red', linestyle='--', label=f'Mean contamination index: {mean_value:.5f}')

    # 设置标题和标签
    ax.set_title(title)
    ax.set_xlabel('')
    ax.set_ylabel('Contamination Index')
    ax.legend()

    # 显示图形
    plt.xticks(rotation=90)
    plt.tight_layout()

    return fig


def main(base_save_dir, pg_matrix, sampleinfo, contaminationType, uniprot, gene):
    # 创建 ArgumentParser 对象
    # outfile_prefix = os.path.join(base_save_dir, 'result')
    # if not os.path.exists(outfile_prefix):
    #     os.makedirs(outfile_prefix)

    if uniprot:
        col = uniprot
        index_type = "uniprot"
    elif gene:
        col = gene
        index_type = "gene"

    ## 用uniprot或者genen name读入血小板/红细胞/凝血marker
    if index_type == "uniprot":
        platelet_marker = ["P48059", "P04899", "Q15286", "P06576", "Q93084", "P05556", "P48735", "P40197", "P24557", "P05141", "P30048", "P14625", "P23229", "O60610", "O43182", "P10809", "P10720", "Q8TC12", "Q9NR12", "P61106", "Q86YW5", "Q96AP7", "Q14165", "Q9BX67", "Q9Y277", "P40926", "P45880", "P16615", "P62873", "P21796"]
        rbc_marker = ["O75955", "P00558", "P00915", "P02042", "P02549", "P02730", "P04406", "P11142", "P11171", "P11277", "P13987", "P16157", "P16452", "P23528", "P27105", "P30043", "P35611", "P35612", "P37840", "P48426", "P53396", "P60891", "P61224", "P62826", "P62937", "P68871", "P69905", "Q00013", "Q9NTK5", "Q9UQ80"]
        coagulation_marker = ["P02679", "P02675", "P02671"]
    elif index_type == "gene":
        platelet_marker = ["LIMS1", "GNAI2", "RAB35", "ATP5F1B", "ATP2A3", "ITGB1", "IDH2", "GP5", "TBXAS1", "SLC25A5", "PRDX3", "HSP90B1", "ITGA6", "DIAPH1", "ARHGAP6", "HSPD1", "PF4V1", "RDH11", "PDLIM7", "RAB14", "TREML1", "ESAM", "MLEC", "JAM3", "VDAC3", "MDH2", "VDAC2", "ATP2A2", "GNB1", "VDAC1"]
        rbc_marker = ["FLOT1", "PGK1", "CA1", "HBD", "SPTA1", "SLC4A1", "GAPDH", "HSPA8", "EPB41", "SPTB", "CD59", "ANK1", "EPB42", "CFL1", "STOM", "BLVRB", "ADD1", "ADD2", "SNCA", "PIP4K2A", "ACLY", "PRPS1", "RAP1B", "RAN", "PPIA", "HBB", "HBA1", "MPP1", "OLA1", "PA2G4"]
        coagulation_marker = ["FGA", "FGB", "FGG"]

    ## 处理蛋白矩阵和样本信息表 ##
    pg_matrix = read_file(pg_matrix)
    sampleinfo = read_file(sampleinfo)
    ## 矩阵预处理
    pg_matrix[col] = pg_matrix[col].str.upper()
    pg_matrix.index = pg_matrix[col]
    pg_matrix.index.name = None
    ## 去除单个肽段对应到多个蛋白的情况
    pg_matrix = pg_matrix[~pg_matrix[col].str.contains(';')]
    ## 这一部分可以删掉
    pg_matrix = pg_matrix[sampleinfo['SampleID']]

    ## 移除缺失率>50%的蛋白
    missing_rate = pg_matrix.iloc[:, :].isnull().mean(axis=1)
    pg_matrix = pg_matrix[missing_rate < 0.5]
    ## 这一部分可以删掉

    outfile_prefix = os.path.join(base_save_dir, 'result')
    if os.path.exists(outfile_prefix):
        shutil.rmtree(outfile_prefix)

    os.makedirs(outfile_prefix, exist_ok=True)

    contaminationType = contaminationType.split(",")

    with PdfPages(os.path.join(outfile_prefix, 'OmniProt-based Plasma Contamination Index Profile.pdf')) as pdf:
        if "platelet" in contaminationType:
            ## 血小板/红细胞的marker gene是随着血小板/红细胞数量增加而丰度增加的蛋白
            platelet_contamination_index = contamination_index(platelet_marker, pg_matrix)
            platelet_contamination_index.to_excel(os.path.join(outfile_prefix, 'Platelet_contamination_index_profile.xlsx')  ,index = False)
            platelet_fig = plot_dataframe(platelet_contamination_index, None, title='Platelet contamination index',width = len(sampleinfo) / 20,color = "#2d56a6")
            pdf.savefig(platelet_fig)  # 将图形写入 PDF 文件
            platelet_fig.savefig(os.path.join(outfile_prefix, 'Platelet_contamination_index_profile.png'))
            plt.close(platelet_fig)

        if "rbc" in contaminationType:
            rbc_contamination_index = contamination_index(rbc_marker, pg_matrix)
            rbc_contamination_index.to_excel(os.path.join(outfile_prefix, 'Erythrocyte_contamination_index_profile.xlsx'),index = False)
            rbc_fig = plot_dataframe(rbc_contamination_index, None, title='Erythrocyte contamination index',width = len(sampleinfo) / 20,color = "#8e2122")
            pdf.savefig(rbc_fig)  # 将图形写入 PDF 文件
            rbc_fig.savefig(os.path.join(outfile_prefix, 'Erythrocyte_contamination_index_profile.png'))
            plt.close(rbc_fig)

        if "coagulation" in contaminationType:
            ## 凝血的三个基marker在凝血发生后显著丰度下降；因此这三个基因在样本中表达越低，凝血越严重
            coagulation_contamination_index = contamination_index(coagulation_marker, pg_matrix,reverse = True)
            coagulation_contamination_index.to_excel(os.path.join(outfile_prefix, 'Coagulation_contamination_index_profile.xlsx'),index = False)
            coag_fig = plot_dataframe(coagulation_contamination_index, None, title='Coagulation contamination index',width = len(sampleinfo) / 20,color = "#e19273" )
            pdf.savefig(coag_fig)  # 将图形写入 PDF 文件
            coag_fig.savefig(os.path.join(outfile_prefix, 'Coagulation_contamination_index_profile.png'))
            plt.close(coag_fig)

if __name__ == "__main__":
    base_save_dir = r'E:\data\omniProt\upload\5b58eaec72b4609d1bbaef9fda42df3c'
    pg_matrix = r'E:\data\omniProt\upload\5b58eaec72b4609d1bbaef9fda42df3c\wosp24266_pg_matrix.csv'
    sampleinfo = r'E:\data\omniProt\upload\5b58eaec72b4609d1bbaef9fda42df3c\wosp24266_sampleinfo.xlsx'
    contaminationType = 'platelet,rbc,coagulation'
    uniprot = 'Protein.Group'
    main(base_save_dir, pg_matrix, sampleinfo, contaminationType, uniprot, '')
