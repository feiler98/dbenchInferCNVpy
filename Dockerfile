FROM feiler98/pyomics_fedora

WORKDIR /scratch/mpt/feiler/dbenchInfercnvpy
COPY . .

RUN pip install --no-cache-dir -r requirements.txt
CMD ["python3", "/home/f/feiler/dbenchCopyKat/run_infercnv.py"]