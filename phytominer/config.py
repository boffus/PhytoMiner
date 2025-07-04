import os

# External Service URLs
PHYTOZOME_SERVICE_URL = os.getenv(
    "PHYTOZOME_SERVICE_URL",
    "https://phytozome-next.jgi.doe.gov/phytomine/service"
)

# Magic Numbers and Defaults
DEFAULT_SLEEP_SECONDS = float(os.getenv("PHYTO_SLEEP_SECONDS", "1.5"))
DEFAULT_CHUNK_SIZE = int(os.getenv("PHYTO_CHUNK_SIZE", "16"))
DEFAULT_MAX_WORKERS = int(os.getenv("PHYTO_MAX_WORKERS", "8"))
