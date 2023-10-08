# PCmaster
Plant Cell master    
    
Available now：    
PCmaster_anno: Plant Cell master for automatic annotation of cell types    
    
In preparation：    
PCmaster_seg: Plant Cell master for cell segment of embyro/leaf/...    
#### ####
#### ####
------------------------------------------------------------------------------------------------------------------------------------------
#### ####
#### ####
### PCmaster_anno ###
### Version 1.0.0 ###
---
#### [2023-10-08] ####
#### Installation with TOSICA  ####
```
# Please install conda first
conda create -n pcmaster_anno_0 --offline
conda activate pcmaster_anno_0
conda install -c conda-forge python=3.8 scanpy
conda install pytorch=1.7.1 torchvision=0.8.2 torchaudio=0.7.2 cudatoolkit=10.1 -c pytorch
# Download TOSICA-main from https://github.com/JackieHanLab/TOSICA
cd TOSICA-main
pip install .
pip install d2l jupyter==1.0.0
python -m ipykernel install --user --name=pcmaster_anno_0 --display-name='Environment (pcmaster_anno_0)'
jupyter-notebook --ip=xxx.yyy.zzz.aaa --no-browser
# Open the web browser and go to http://xxx.yyy.zzz.aaa:cccc or https://xxx.yyy.zzz.aaa:cccc
```
#### Installation with docker  ####
```
# Please install docker first
# Download the docker image file from https://drive.google.com/file/d/1236fiXdtY4WwtrU3zPVqE-eJcJf5m5vI/view?usp=drive_link
docker load --input bioinplant_pcmaster_anno_0_23_9_28.tar
docker images
docker run -it --name pcmaster_anno_0_23_9_28 --gpus all -p 8996:8997 bioinplant/pcmaster_anno_0_23_9_28 /bin/bash
jupyter-notebook --ip=xxx.yyy.zzz.aaa --no-browser
# Open the web browser and go to http://xxx.yyy.zzz.aaa:8996 or https://xxx.yyy.zzz.aaa:8996
```
#### Simple usage (auto annotation with ref datasets and deep learning models)  ####
```
# Other files such as PCmaster_anno_0_23_10_3.py, plant_marker_gene_list.txt and the pth file of resnet in this git project are also needed
# More details in PCmaster_anno_0_guide_1_mainly_in_Chinese.docx
with open('PCmaster_anno_0_23_10_3.py','r') as f:
    exec(f.read())
pcma = PCmaster_anno_0()
pcma.auto_annotation_with_deep_learning_0(original_obj=the_original_obj,gpu_code = gpu_code_n,
                                          learning_rate=the_learning_rate,
                                          epochs=the_epochs,
                                          batch_size = the_batch_size,dropout = the_dropout,
                                          num_workers = the_num_workers
                                          )
```
```
On the test dataset
Accuracy: 0.9025
Macro Precision: 0.901893589299713
Macro Recall: 0.9025
Macro F1 Score: 0.9016473881006721

On a dataset which is not involved in model training
Accuracy: 0.6825028968713789
Macro Precision: 0.7359756337721365
Macro Recall: 0.6869098209694935
Macro F1 Score: 0.6703407567251473

The score is higher than SingleR
[1] "obj"
An object of class Seurat 
5859 features across 7462 samples within 1 assay 
Active assay: RNA (5859 features, 0 variable features)
[1] 100
[1] "count/length(big_df$true)"
[1] 0.5844504

```
#### Cell type true ####
![image](https://github.com/bioinplant/PCmaster/blob/main/celltype-true.png)
#### Cell type pred ####
![image](https://github.com/bioinplant/PCmaster/blob/main/celltype-pred.png)
---
#### [2023-09-11] ####
#### PCmaster_anno_0_23_9_11.py has been uploaded. The default 'cluster_n_neighbors' has been changed from 20 to 10, which is the same as scanpy. And the function of auto annotation with reference datasets has been improved. You can experience improved functions by replacing old contents in original files with new contents.  ####
---
#### [2023-08-15] ####
#### PCmaster_anno_0_23_8_15.py and plant_marker_gene_list_23_8_15.txt have been uploaded. You can experience improved functions by replacing old contents in original files with new contents.  ####
---
#### [2023-07-05] ####
#### Currently, it is recommended to follow the guidance in PCmaster_anno_0_guide_1_mainly_in_Chinese.docx and install the environment from scratch. ####
#### Installing via requirements.txt seems to have compatibility issues at the moment. (It seems to be affected by conda version and network environment, too.) ####
#### Installing through the compressed file PCmaster_anno_test_1.tar.gz (There are also guides in PCmaster_anno_0_guide_1_mainly_in_Chinese.docx.) via [this link](https://pan.baidu.com/s/1p1im8oCfebzGjzptk7PSmQ?pwd=nnvr) is the most convenient way, but there may also be compatibility issues. ####
---
#### [2023-06-05] ####
#### It is recommended to put your ipynb, "PCmaster/PCmaster_anno/PCmaster_anno_0_copy1.py" and "PCmaster/PCmaster_anno/plant_marker_gene_list.txt" into the same folder. ####
---      
#### If you don't want to reinstall the conda environment from scratch, you can download the compressed package through [this link](https://pan.baidu.com/s/1p1im8oCfebzGjzptk7PSmQ?pwd=nnvr), which also contains some files and codes that can be used for testing. ####
#### You can temporarily get some instructions through "PCmaster/PCmaster_anno/PCmaster_anno_0_guide_1_mainly_in_Chinese.docx". We are trying to make a guide website. #### 
You can see some analysis examples in "PCmaster/PCmaster_anno/PCmaster_anno_example_1.ipynb".    
If you want to use the latest version, please download the latest 'PCmaster_anno_0_XX_XX_XX.py' file and copy the code in it to replace the code in 'PCmaster_anno_0_copy1.py'.    
Please note that the latest version of the code may cause some undetected bugs.     

SCAPP is the old version of PCmaster_anno.    
https://github.com/shlin0415/SCAPP    

Sincerely thanks to the contributors of python packages such as d2l, doubletdetection, harmonypy, numpy, openmmlab, pandas, scanpy, seaborn, scrublet, scikit-learn, torch, tosica, etc.    

If you have some questions, please send email to 12216017@zju.edu.cn.    
