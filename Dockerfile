FROM feiler98/pyomics_fedora

WORKDIR /scratch/tmp/feiler/dbenchInferCNVpy
COPY . .

RUN pip install --no-cache-dir -r requirements.txt
CMD ["python3", "/home/f/feiler/dbenchInferCNVpy/run_infercnv.py"]