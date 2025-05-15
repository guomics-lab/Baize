
<p align="center" style="margin-bottom: 0px !important;">
  <img src="https://guomics.com/wp-content/uploads/2025/02/Baize-01-1-800x448.jpg" width="200" height="112">
</p>
Baize, a web-based software tool for rapid evaluation of sample contamination in three critical dimensions: platelet contamination, erythrocyte lysis, and residual coagulation protein carryover. The algorithm calculates cell-type-specific contamination indices by normalizing the summed intensity of marker proteins against the total intensity of plasma proteome (Contamination Index = Marker Protein Intensities / Σ All Plasma Protein Intensities). Users can submit a protein matrix with a sample annotation table, and then Baize will output an evaluation of contamination for every plasma sample. Baize is freely accessible at https://www.guomics.com/Baize.

### Output
The output provides each sample contamination profile, enabling rapid quality control decisions without requiring manual intervention.

## Installation
If you want to use Baize by source code, you can install python and install requirements package.

### Prerequisites
Please make sure you have a valid installation of conda or miniconda. We recommend setting up miniconda as described on their website.

```shell
git clone https://github.com/guomics-lab/Baize.git
cd Baize
```

```shell
conda create -n Baize python=3.9
conda activate Baize
```

```shell
pip install -r requirements.txt
```

Run
```shell
python web_main.py
```

Open your browser and go to http://127.0.0.1:5000/Baize