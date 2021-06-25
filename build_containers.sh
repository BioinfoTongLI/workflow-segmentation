set -e pipefail

docker build -t eu.gcr.io/imaging-gpu-eval/extract_channels:latest -f containers/Dockerfile.extract_channels .

docker build -t eu.gcr.io/imaging-gpu-eval/expand_labels:latest -f containers/Dockerfile.expand_labels .

docker build -t eu.gcr.io/imaging-gpu-eval/filter_labels:latest -f containers/Dockerfile.filter_labels .

docker build -t eu.gcr.io/imaging-gpu-eval/cellpose:latest -f containers/Dockerfile.cellpose .

docker build -t eu.gcr.io/imaging-gpu-eval/segmentation_helper:latest -f containers/Dockerfile.helper .
