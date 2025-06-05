# launching:
projects=/storage/replicated/cirad/projects/GE2POP
images_path=${projects}/APPTAINER_IMAGES
singularity run ${images_path}/python3.11_libs/v0.1/env_python3.11_libs.sif ./check_gene_overlaps.py test_files/test.list

