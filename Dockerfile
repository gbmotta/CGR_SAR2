FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements-streamlit.txt requirements-api.txt ./
RUN pip install --no-cache-dir -r requirements-streamlit.txt -r requirements-api.txt

COPY setup.sh ./
COPY cgr_sar2 ./cgr_sar2
COPY streamlit_app.py ./
RUN chmod +x setup.sh && ./setup.sh

EXPOSE 8000 8501
CMD ["python", "-m", "cgr_sar2", "api", "--host", "0.0.0.0", "--port", "8000"]
