FROM continuumio/miniconda3:4.5.11

RUN apt-get update && \
    apt-get -y dist-upgrade && \
    apt-get -y --no-install-recommends install build-essential && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN useradd -m charge_assign && \
    echo '. /opt/conda/etc/profile.d/conda.sh' >>/home/charge_assign/.bashrc

WORKDIR /home/charge_assign
USER charge_assign

RUN bash -c ' \
    . /opt/conda/etc/profile.d/conda.sh && \
    conda create -n charge_assign && \
    conda install -c rdkit -c conda-forge -n charge_assign rdkit nauty \
    '

USER root
COPY charge /home/charge_assign/charge
COPY charge_server /home/charge_assign/charge_server
COPY scripts /home/charge_assign/scripts
COPY setup.py /home/charge_assign/setup.py
COPY setup.cfg /home/charge_assign/setup.cfg
RUN find /home/charge_assign/charge /home/charge_assign/charge_server -name '__pycache__' -depth -exec rm -r \{\} \;
RUN chown -R charge_assign:charge_assign charge charge_server setup.py setup.cfg

USER charge_assign
RUN bash -c ' \
    source activate charge_assign && \
    cd /home/charge_assign && \
    pip install . \
    '

EXPOSE 8080
WORKDIR /home/charge_assign
CMD bash -c 'source activate charge_assign && python3 -m charge_server -r /home/charge_assign/repo.zip'
