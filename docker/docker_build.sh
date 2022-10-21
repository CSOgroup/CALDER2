#!/bin/bash

DOCKERHUB_USERNAME="lucananni93"
PACKAGE_NAME="calder2"
GITHUB_REPO="CSOgroup/CALDER2"
DOCKER_FILE="docker/Dockerfile"
LATEST_TAG=$(curl --silent "https://api.github.com/repos/${GITHUB_REPO}/releases/latest" | grep '"tag_name"' | cut -d':' -f2 | tr -d '," ')

echo "Building Docker for ${GITHUB_REPO} at ${DOCKERHUB_USERNAME}/${PACKAGE_NAME}:${LATEST_TAG}"

docker build . \
            -t ${DOCKERHUB_USERNAME}/${PACKAGE_NAME}:${LATEST_TAG} \
            -t ${DOCKERHUB_USERNAME}/${PACKAGE_NAME}:latest \
            -f ${DOCKER_FILE} \
            --build-arg tag_name=${LATEST_TAG} \
            --build-arg repo_name=${GITHUB_REPO}

docker push ${DOCKERHUB_USERNAME}/${PACKAGE_NAME}:${LATEST_TAG}
docker push ${DOCKERHUB_USERNAME}/${PACKAGE_NAME}:latest