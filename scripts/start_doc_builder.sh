#!/bin/bash
if ! docker image inspect cfd-doc-image &> /dev/null; then
    echo "Building cfd-doc-image"
    docker build doc -t cfd-doc-image
fi
docker run -it --rm --name cfd-doc -v "$(pwd):/workspace" cfd-doc-image /bin/bash
