# PCmaster
Plant Cell master    
    
Available now：    
PCmaster_anno: Plant Cell master for automatic annotation of cell types      

#### ####
#### ####
------------------------------------------------------------------------------------------------------------------------------------------
#### ####
#### ####
### PCmaster_anno ###
### Version 1.0.0 ###
---
#### Installation ####
```
# The versions may not be exactly the same.
# Please install conda first
conda create -n PCmaster_anno_0 --offline
conda activate PCmaster_anno_0
conda install python=3.8

pip install jupyter d2l torch torchvision
# jupyter                   1.0.0                    pypi_0    pypi
# d2l                       0.17.6                   pypi_0    pypi
# torch                     1.13.1                   pypi_0    pypi
# torchvision               0.14.1                   pypi_0    pypi

conda install -c conda-forge scanpy
# scanpy                    1.9.2              pyhd8ed1ab_0    conda-forge

python -m ipykernel install --user --name=PCmaster_anno_0 --display-name='Environment (PCmaster_anno_0)'
conda install -c conda-forge pytables
# pytables                  3.7.0            py38hf19a122_1

conda install -c conda-forge leidenalg
# leidenalg                 0.9.1            py38h8dc9893_0    conda-forge

conda install -c bioconda harmonypy
# harmonypy                 0.0.6              pyhdfd78af_0    bioconda

pip install optuna -i https://pypi.tuna.tsinghua.edu.cn/simple
# optuna                   3.1.0

pip install optuna-dashboard -i https://pypi.tuna.tsinghua.edu.cn/simple
# optuna-dashboard         0.8.1

pip install plotly -i https://pypi.tuna.tsinghua.edu.cn/simple
# plotly                   5.13.1

conda install -c bioconda -c conda-forge scrublet
# scrublet                  0.2.3              pyh5e36f6f_1    bioconda

pip install doubletdetection
# doubletdetection         4.2

conda install -c r -c conda-forge r-irkernel
# r-irkernel                1.3               r40hc72bb7e_0    conda-forge

conda install -c conda-forge conda-pack
# conda-forge/noarch::conda-pack-0.7.0-pyh6c4a22f_0

# maybe needed
# pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

```
#### Installation with docker  ####
```
# Please install docker first
# Download the docker image file from https://drive.google.com/file/d/1236fiXdtY4WwtrU3zPVqE-eJcJf5m5vI/view?usp=drive_link
docker load --input bioinplant_pcmaster_anno_0_23_9_28.tar
docker images
docker run -it --name pcmaster_anno_0_23_9_28 --gpus all -p 8996:8997 bioinplant/pcmaster_anno_0_23_9_28 /bin/bash

ls
cd home
conda activate pcmaster_anno_0
jupyter-notebook --ip=xxx.yyy.zzz.aaa --no-browser
# Open the web browser and go to http://xxx.yyy.zzz.aaa:8996 or https://xxx.yyy.zzz.aaa:8996
```
---
#### Usage ####
#### See ipynbs in 'Tmp_tutorial'. ####
#### The refrence datasets can be downloaded from CNGBdb (https://db.cngb.org/), EBI, NCBI, scPlantDB (https://biobigdata.nju.edu.cn/scplantdb/dataset) and STOmics DB (https://db.cngb.org/stomics/). ####
#### The 'resnet18.pth' can be downloaded from https://drive.google.com/file/d/1dZful0MsOm73hodk2Nt1zBV7li7Ypmhb/view?usp=drive_link. ####
---
SCAPP is the old version of PCmaster_anno.    
https://github.com/shlin0415/SCAPP    

For PCmaster_anno (v.1.0.0), Miniconda3 (v.23.7.3), Jupyter (v.1.0.0) and Python (v.3.8.16) are utilized to build the analysis platform on CentOS Linux (release 7.4.1708) with one NVIDIA Tesla V100-SXM2-32GB GPU (also workable on Windows with one NVIDIA GeForce RTX 3080 Laptop GPU). D2l (v.0.17.6), Numpy (v.1.23.5), Pandas (v.2.0.3), SCANPY (v.1.9.2), Scikit-learn (v.1.2.1), Torch (v.1.13.1) and their dependent packages are mainly used for analysis and annotation. Matplotlib (v.3.5.3), Plotly (v.5.9.0), Plottable (v.0.1.5), Seaborn (v.0.13.2) and their dependent packages are mainly utilized for visualization.

Sincerely thanks to the contributors of these packages: cell blast, ciform, conda, d2l, docker, doubletdetection, harmonypy, jupyter, matplotlib, nrtpredictor, numpy, openmmlab, pandas, plotly, plottable, python, scanpy, scgpt, scikit-learn, scplant, scrublet, scvi, seaborn, seurat, singler, torch, tosica, etc.

If you have any questions, please send the email to 12216017@zju.edu.cn.
