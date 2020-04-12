FROM osgeo/gdal:ubuntu-small-3.0.3 AS compile-image

RUN apt-get install -y python3-setuptools

RUN useradd --create-home appuser
WORKDIR /home/appuser
USER appuser

RUN mkdir dempy
ADD dempy ./dempy
COPY setup.py ./dempy

RUN python3 ./dempy/setup.py install --user

ADD example.py ./dempy
CMD ["python3"]
