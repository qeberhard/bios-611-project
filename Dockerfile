FROM rocker/verse
 RUN apt update && apt-get install -y openssh-server python3-pip
 RUN pip3 install sklearn pandas numpy
 RUN R -e "install.packages(\"shiny\")"
 RUN R -e "install.packages(\"deSolve\")"
 RUN R -e "install.packages(\"signal\")"
 RUN R -e "install.packages(\"reticulate\")";