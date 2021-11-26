# Dockerfile
# syntax=docker/dockerfile:1

FROM python:3.8-slim-buster

RUN mkdir /app
WORKDIR /app
RUN chmod -R 777 /app

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY . .

ENV FLASK_APP main.py

ENTRYPOINT ["python", "-m", "flask", "run", "--host=0.0.0.0"]
