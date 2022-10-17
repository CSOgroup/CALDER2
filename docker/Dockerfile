# We need Bioconductor, so we inherit from their repositorycd .
FROM bioconductor/bioconductor_docker

ARG tag_name
ARG repo_name

ENV github_tag=$tag_name
ENV github_repo=$repo_name

RUN echo "Building repo ${github_repo} with tag: ${github_tag}"

################
# CALDER SETUP #
################

# Get latest release tag
RUN wget https://github.com/${github_repo}/archive/refs/tags/${github_tag}.tar.gz
# Untar the release
RUN tar xvfs ${github_tag}.tar.gz
# Adding calder script to bin
ENV calder_cmd=CALDER2-${github_tag}/scripts/calder
RUN ln -s $(realpath ${calder_cmd}) /bin/
# Removing tar file
RUN rm ${github_tag}.tar.gz

# # Install the package
RUN R -e "remotes::install_local('CALDER2-${github_tag}', dependencies=T); if (!library(CALDER, logical.return=T)) quit(status=10)"

# # Install other packages
RUN pip install cooler 

CMD ["bash"]
