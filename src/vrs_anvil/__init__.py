import json
import logging
import os
import subprocess
from typing import Optional, Generator, Any
import zipfile

import psutil
from diskcache import Cache
from dotenv import load_dotenv
from ga4gh.vrs import models as VRS
from ga4gh.vrs.dataproxy import _DataProxy, create_dataproxy
from ga4gh.vrs.extras.translator import AlleleTranslator
from pathlib import Path
from glom import glom
from pydantic import BaseModel, model_validator
import requests
import yaml

load_dotenv()

_logger = logging.getLogger(__name__)
LOGGED_ALREADY = set()
METAKB_API = "https://pediatric.metakb.org/api"


manifest: "Manifest" = None

# TODO - read from manifest
gigabytes = 20
bytes_in_a_gigabyte = 1024**3  # 1 gigabyte = 1024^3 bytes
cache_size_limit = gigabytes * bytes_in_a_gigabyte


def get_cache_directory(cache_dir: str, cache_name: str) -> str:
    """Return the cache directory."""
    return str(Path(cache_dir) / cache_name)


class CachingAlleleTranslator(AlleleTranslator):
    """A subclass of AlleleTranslator that uses cache results and adds a method to run in a threaded fashion."""

    _cache: Cache = None

    def __init__(self, data_proxy: _DataProxy, normalize: bool = False):
        super().__init__(data_proxy)
        self.normalize = normalize
        self._cache = None
        if manifest and manifest.cache_enabled:
            self._cache = Cache(
                directory=get_cache_directory(
                    manifest.cache_directory, "allele_translator"
                ),
                size_limit=cache_size_limit,
            )
        else:
            _logger.info("Cache is not enabled")

    def translate_from(self, var, fmt=None, **kwargs):
        """Check and update cache"""

        if self._cache is not None:
            key = f"{var}-{fmt}"
            if key in self._cache:
                return self._cache[key]

        allele = super().translate_from(var, fmt=fmt, **kwargs)

        assert isinstance(
            allele, VRS.Allele
        ), f"Allele is not the expected Pydantic Model {type(allele)}: {allele}"

        if self._cache is not None:
            self._cache[key] = allele.id

        return allele.id


def caching_allele_translator_factory(
    normalize: bool = True, seqrepo_uri: str | None = None
):
    """Return a CachingAlleleTranslator instance with a SeqRepo dataproxy."""
    if manifest is not None:
        dp = create_dataproxy(seqrepo_uri or manifest.seqrepo_uri)
    else:
        dp = create_dataproxy(seqrepo_uri)
    translator = CachingAlleleTranslator(dp)
    translator.normalize = normalize
    return translator


def generate_gnomad_ids(vcf_line, compute_for_ref: bool = True) -> list[str]:
    """Assuming a standard VCF format with tab-separated fields, generate a gnomAD-like ID from a VCF line.
    see https://github.com/ga4gh/vrs-python/blob/main/src/ga4gh/vrs/extras/vcf_annotation.py#L386-L411
    """
    fields = vcf_line.strip().split("\t")
    gnomad_ids = []
    # Extract relevant information (you may need to adjust these indices based on your VCF format)
    chromosome = fields[0]
    position = fields[1]
    reference_allele = fields[3]
    alternate_allele = fields[4]

    gnomad_loc = f"{chromosome}-{position}"
    if compute_for_ref:
        gnomad_ids.append(f"{gnomad_loc}-{reference_allele}-{reference_allele}")
    for alt in alternate_allele.split(","):
        alt = alt.strip()
        # TODO - Should we be raising a ValueError hear and let the caller do the logging?
        # TODO - Should this be a config in the manifest?
        # ['<INS>', '<DEL>', '<DUP>', '<INV>', '<CNV>', '<DUP:TANDEM>', '<DUP:INT>', '<DUP:EXT>', '*']
        invalid_alts = ["INS", "DEL", "DUP", "INV", "CNV", "TANDEM", "INT", "EXT", "*"]
        is_valid = True
        for invalid_alt in invalid_alts:
            if invalid_alt in alt:
                is_valid = False
                _ = f"Invalid alt found: {alt}"
                if _ not in LOGGED_ALREADY:
                    LOGGED_ALREADY.add(_)
                    _logger.error(_)
                break
        if is_valid:
            gnomad_ids.append(f"{gnomad_loc}-{reference_allele}-{alt}")

    return gnomad_ids


def params_from_vcf(path, limit=None) -> Generator[dict, None, None]:
    """Open the vcf file, skip headers, yield the first lines as gnomad-like IDs"""
    from vrs_anvil.translator import VCFItem

    c = 0
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            gnomad_ids = generate_gnomad_ids(line)
            for gnomad_id in gnomad_ids:
                yield VCFItem(
                    fmt="gnomad",
                    var=gnomad_id,
                    file_name=path,
                    line_number=c,
                    identifier=None,
                )  # TODO - add identifier
            c += 1
            if limit and c > limit:
                break


def find_items_with_key(dictionary, key_to_find):
    """Find all items in a dictionary that have a specific key."""
    result = glom(dictionary, f"**.{key_to_find}")
    return result


class MetaKBProxy(BaseModel):
    """A proxy for the MetaKB, maintains a cache of VRS ids."""

    metakb_path: Path
    cache_path: Path
    _cache: Optional[Cache] = None

    def __init__(self, metakb_path: Path, cache_path: Path, cache: Cache = None):
        super().__init__(metakb_path=metakb_path, cache_path=cache_path, _cache=cache)
        if cache is None:
            reload_cache = False
            if not (metakb_path / "cache").is_dir():
                reload_cache = True
            cache = Cache(directory=get_cache_directory(cache_path, "metakb"))
            # cache.stats(enable=True) # drives up disk usage
            if reload_cache:
                for _ in metakb_ids(metakb_path):
                    cache.set(_, True)
        self._cache = cache

    def get(self, vrs_id: str) -> bool:
        """Get the vrs_id from the cache."""
        return self._cache.get(vrs_id, False)


def metakb_ids(metakb_path) -> Generator[str, None, None]:
    """Find all the applicable vrs ids in the metakb files."""
    if len(list(Path(metakb_path).glob("*.json"))) == 0:
        _get_metakb_models(metakb_path)

    for file_name in Path(metakb_path).glob("*.json"):
        if file_name.is_file():
            with open(file_name, "r") as file:
                data = json.loads(file.read())
                yield from (
                    [
                        _
                        for _ in find_items_with_key(data, "id")
                        if _.startswith("ga4gh:VA")
                    ]
                )


def _get_metakb_models(metakb_path):
    def _download_s3(url: str, outfile_path: Path) -> None:
        """Download objects from public s3 bucket

        :param url: URL for metakb file in s3 bucket
        :param outfile_path: Path where file should be saved
        """
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(outfile_path, "wb") as h:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        h.write(chunk)

    Path(metakb_path).mkdir(exist_ok=True)

    date = "20240305"
    json_files = [f"civic_cdm_{date}.json", f"moa_cdm_{date}.json"]

    for json_file in json_files:
        json_path = f"{metakb_path}/{json_file}"

        url = (
            "https://vicc-metakb.s3.us-east-2.amazonaws.com"
            + f"/cdm/{date}/{json_file}.zip"
        )
        zip_path = f"{json_path}.zip"

        _download_s3(url, zip_path)
        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(metakb_path)
        os.remove(zip_path)


class Manifest(BaseModel):
    """
    A class to represent the manifest file.
    Defaults to the values below if not provided in the manifest.yaml
    """

    cache_directory: str = "cache/"
    """Path to the cache directory, defaults to cache/ (relative to the root of the repository)"""

    num_threads: int = 2
    """Number of threads to use for processing, defaults to 2"""

    # TODO: not implemented
    annotate_vcfs: bool = False
    """Should we create new VCFs with annotations. FOR FUTURE USE"""

    state_directory: str = "state/"
    """where to store the state of the application, log files, etc."""

    vcf_files: list[str]
    """The local file paths or URLs to vcf files to be processed"""
    # TODO - 2x check why local files need to be absolute paths

    work_directory: str = "work/"
    """The directory to store intermediate files"""

    seqrepo_uri: str | None = None
    """Description of the available SeqRepo resource.

    Currently accepted URI schemes:

    * seqrepo+file:///path/to/seqrepo/root
    * seqrepo+:../relative/path/to/seqrepo/root
    * seqrepo+http://localhost:5000/seqrepo
    * seqrepo+https://somewhere:5000/seqrepo
    """

    normalize: bool = True
    """Normalize the VRS ids"""

    limit: Optional[int] = None
    """Stop processing after this many lines"""

    cache_enabled: Optional[bool] = True
    """Cache results"""

    compute_for_ref: Optional[bool] = False
    """Compute reference allele"""

    estimated_vcf_lines: Optional[int] = 4000000
    """How many lines per vcf file?  Used for progress bar"""

    metakb_directory: str = "metakb/"
    """Where the CDM files are located.  This is a directory containing json files"""

    disable_progress_bars: Optional[bool] = False

    @model_validator(mode="after")
    def check_paths(self) -> "Manifest":
        """Post init method to set the cache directory."""
        self.work_directory = str(Path(self.work_directory).expanduser())
        self.cache_directory = str(Path(self.cache_directory).expanduser())
        self.state_directory = str(Path(self.state_directory).expanduser())
        self.metakb_directory = str(Path(self.metakb_directory).expanduser())

        if not Path(self.metakb_directory).exists():
            raise ValueError("MetaKB directory does not exist")

        for _ in ["work_directory", "cache_directory", "state_directory"]:
            if not Path(getattr(self, _)).exists():
                Path(getattr(self, _)).mkdir(parents=True, exist_ok=True)
                _logger.debug(f"Created directory {getattr(self, _)}")

        return self


def query_metakb(
    variation_str: str | None = None,
    disease: str | None = None,
    therapy: str | None = None,
    gene: str | None = None,
    statement_id: str | None = None,
    log: bool = False,
) -> dict:
    """
    Query MetaKB API for variant statements.

    This is a more Pythonic wrapper around the MetaKB /search/statements endpoint.
    While the metakb package provides a direct Python API (search_statements), it requires
    a Neo4j database connection. Since we're querying the hosted API, we use HTTP requests
    but validate the response using the Pydantic models.

    Args:
        variation_str: Variation identifier (VRS ID, HGVS, gene symbol, etc.)
        disease: Optional disease filter
        therapy: Optional therapy filter
        gene: Optional gene filter
        statement_id: Optional statement ID filter (e.g., "civic.eid:102", "civic.aid:7")
        log: Whether to log the request (for debugging)

    Returns:
        Dictionary containing the search results (compatible with SearchStatementsResponse)

    Example:
        >>> # Query by freetext variation
        >>> results = query_metakb(variation_str="BRAF V600E")
        >>>
        >>> # Query by VRS ID
        >>> results = query_metakb(variation_str="ga4gh:VA.j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L")
        >>>
        >>> # Query by statement ID
        >>> results = query_metakb(statement_id="civic.eid:102")
        >>>
        >>> # Combine filters
        >>> results = query_metakb(variation_str="BRAF V600E", disease="melanoma")

    Raises:
        ValueError: If no query parameters are provided or if API returns an error
    """
    # Build query parameters
    params = {}
    if variation_str:
        params["variation"] = variation_str
    if disease:
        params["disease"] = disease
    if therapy:
        params["therapy"] = therapy
    if gene:
        params["gene"] = gene
    if statement_id:
        params["statement_id"] = statement_id

    # Validate that at least one parameter is provided
    if not params:
        raise ValueError(
            "At least one query parameter must be provided (variation, disease, therapy, gene, or statement_id)"
        )

    # Make the API request
    response = requests.get(
        f"{METAKB_API}/search/statements",
        params=params,
        headers={"Accept": "application/json"},
    )

    # Handle errors
    if response.status_code >= 400:
        error_msg = f"API error: {response.text} ({response.status_code})"
        _logger.error(error_msg)
        raise ValueError(error_msg)

    # Parse and return JSON response
    # The response can be validated with SearchStatementsResponse if needed:
    # validated = SearchStatementsResponse(limit=None, **response.json())
    return response.json()


def run_command_in_background(command) -> Any:
    """Execute the command in the background, return pid."""
    # Detach the process from the parent process (this process)
    if not isinstance(command, list):
        command = command.split()
    return subprocess.Popen(
        command, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )


def get_process_info(pid):
    """Return the process information for the pid."""
    try:
        return psutil.Process(pid)
    except psutil.NoSuchProcess:
        return None


def save_manifest(manifest: Manifest, manifest_path: str):
    """pass in a Manifest and yaml path"""
    with open(manifest_path, "w") as stream:
        yaml.dump(manifest.model_dump(), stream)
