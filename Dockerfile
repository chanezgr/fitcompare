FROM python:3-slim

WORKDIR /app

COPY fitcompare.py /app
COPY config.ini /app

RUN pip install --no-cache fitparse
RUN pip install --no-cache argparse
RUN pip install --no-cache pyyaml
RUN pip install --no-cache pandas
RUN pip install --no-cache matplotlib
RUN pip install --no-cache seaborn
RUN pip install --no-cache scipy

ENTRYPOINT ["python", "-u", "fitcompare.py"]
