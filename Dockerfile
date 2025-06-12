FROM python:3-slim

WORKDIR /app

RUN pip install --no-cache fitparse
RUN pip install --no-cache argparse
RUN pip install --no-cache pyyaml
RUN pip install --no-cache pandas
RUN pip install --no-cache matplotlib
RUN pip install --no-cache seaborn
RUN pip install --no-cache scipy

COPY fitcompare.py /app
COPY fitcompare_advanced.py /app
COPY config.ini /app

ENTRYPOINT ["python", "-u", "fitcompare.py"]
