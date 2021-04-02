# Processing single-cell RNA-seq data for dimension reduction-based analyses using open-source tools

### Install Dependencies
The environment for Step 1 of this protocol require Singularity. The details of this installation can be found in the [Section 1 notebook](https://github.com/Ken-Lau-Lab/STAR_Protocol/blob/master/Single-cell_read_alignment_and_DropEst_library_quantification.ipynb).


The environment for Sections 2 ([Variant 1](https://github.com/Ken-Lau-Lab/STAR_Protocol/blob/master/Variant_1_Heuristic_droplet_filtering.ipynb) and [Variant 2](https://github.com/Ken-Lau-Lab/STAR_Protocol/STAR_Protocol/blob/master/Variant_2_Automated_droplet_filtering_with_dropkick.ipynb)) and [3](https://github.com/Ken-Lau-Lab/STAR_Protocol/blob/master/Post-processing_and_dimension_reduction_structure_preservation_analysis.ipynb) of this protocol can be installed through the included .yml with the following, which will also initialize a jupyter kernel: 
```
conda env create -f qc_pipe_env.yml
python -m ipykernel install --user --name=qc_pipe
```

Otherwise, the requirements.txt file can be manually installed In a working Python virtual environment using `pip`:
```
pip install -r requirements.txt
```
