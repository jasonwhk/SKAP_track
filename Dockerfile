FROM python:3.6-alpine
RUN apk add -U --no-cache gcc build-base \
        python3 python3-dev
RUN pip install --no-cache-dir astropy
ADD track.py /
ADD finals2000A.all /
