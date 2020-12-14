FROM scilus/docker-scilpy:latest

ADD JHU_template_GIN_dil.tar.bz2 /JHU_template_GIN_dil
ADD filtering_lists.tar.bz2 /filtering_lists

WORKDIR /
RUN wget http://trackvis.org/bin/TrackVis_v0.6.1_x86_64.tar.gz
run mkdir TrackVis
RUN tar -xzf TrackVis_v0.6.1_x86_64.tar.gz -C TrackVis
ENV PATH="/TrackVis:${PATH}"
