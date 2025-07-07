import os

# External Service URLs
PHYTOZOME_SERVICE_URL = os.getenv(
    "PHYTOZOME_SERVICE_URL",
    "https://phytozome-next.jgi.doe.gov/phytomine/service"
)

# Magic Numbers & Defaults
DEFAULT_SLEEP_SECONDS = float(os.getenv("PHYTO_SLEEP_SECONDS", "1.5"))
DEFAULT_CHUNK_SIZE = int(os.getenv("PHYTO_CHUNK_SIZE", "16"))
DEFAULT_MAX_WORKERS = int(os.getenv("PHYTO_MAX_WORKERS", "8"))
HOMOLOGS_OUTPUT_FILE = 'step1output.csv'
TSV_DIR = r'C:\Users\toffe\OneDrive - Monash University\1 Pupu\Plant Energy\homologs\genes NdhA-NhdO'
JOIN2_OUTPUT_FILE = 'step2output.csv' 