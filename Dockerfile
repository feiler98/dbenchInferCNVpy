FROM feiler98/pyomics_fedora

RUN mkdir -p /scratch/tmp/feiler/dbenchInferCNVpy
WORKDIR /scratch/tmp/feiler/dbenchInferCNVpy
COPY . .

RUN pip install --no-cache-dir -r requirements.txt
CMD ["python3", "/scratch/tmp/feiler/dbenchInferCNVpy/run_infercnv.py"]