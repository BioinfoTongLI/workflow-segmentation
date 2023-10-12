docker build -t bioinfotongli/workflow-segmentation:latest .
singularity build /lustre/scratch126/cellgen/team283/imaging_sifs/bioinfotongli-workflow-segmentation-latest.sif docker-daemon://bioinfotongli/workflow-segmentation:latest
