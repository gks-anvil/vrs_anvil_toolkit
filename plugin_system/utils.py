import io
import json
import os

from firecloud import api as fapi
import pandas as pd

####################################
#  utilities for transforming data #
####################################

WORKSPACE_ENV_KEYS = ["WORKSPACE_NAMESPACE", "WORKSPACE_NAME"]


def save_dict(d, path):
    # save to disk if specified
    if path:
        if os.path.exists(path):
            raise FileExistsError(
                f"index already exists at path: {path}. Please delete it before continuing"
            )
        else:
            with open(path, "w") as file:
                json.dump(d, file)


def load_dict(path):
    if os.path.exists(path):
        with open(path, "r") as file:
            d = json.load(file)
            return d
    else:
        raise FileNotFoundError(f"path to phenotype index does not exist: {path}")


def terra_data_table_to_dataframe(table_name) -> pd.DataFrame:
    # if unspecified, ensure valid Terra environment
    for env_key in WORKSPACE_ENV_KEYS:
        if env_key not in os.environ:
            raise OSError(
                f"ERROR: No {env_key} key found in environmental variables in the Terra workspace. If you are working in a Terra workspace, please ensure both a WORKSPACE_NAMESPACE and a WORKSPACE_NAME are specified."
            )

    # create dataframe from Terra data table
    # https://github.com/broadinstitute/fiss/blob/master/firecloud/api.py
    try:
        response = fapi.get_entities_tsv(
            os.environ["WORKSPACE_NAMESPACE"],
            os.environ["WORKSPACE_NAME"],
            table_name,
            model="flexible",
        )
        response.raise_for_status()
    except Exception as e:
        if response.json() and e in response.json():
            error_message = response.json()["message"]
        else:
            error_message = e
        print(
            f"Error while loading phenotype data table from workspace: \n{error_message}"
        )

    phenotype_tsv = io.StringIO(response.text)
    return pd.read_csv(phenotype_tsv, sep="\t")


def csv_to_dataframe(path):
    # table path specified, parse using that table
    with open(path, "r") as file:
        if path.endswith(".csv"):
            return pd.read_csv(file)
        elif path.endswith(".tsv"):
            return pd.read_csv(file, sep="\t")
        else:
            raise Exception(
                "Only csv and tsv file types implemented for phenotype table"
            )
