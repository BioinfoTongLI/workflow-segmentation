docker build -t workflow-segmentation:latest .
singularity build /lustre/scratch126/cellgen/team283/imaging_sifs/workflow-segmentation.sif docker-daemon://workflow-segmentation:latest
