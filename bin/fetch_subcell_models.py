#!/usr/bin/env python 

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

from skimage import measure
from skimage import util
import argparse
import datetime
import logging
import os
import sys
import pandas as pd
import requests
import torch
import yaml
from skimage.io import imread
import boto3
from botocore import UNSIGNED
from botocore.config import Config
from botocore.exceptions import ClientError
from urllib.parse import urlparse

import inference
from vit_model import ViTPoolClassifier

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# This is the general configuration variable. We are going to use the special key "log" in the dictionary to use the log in our code
config = {"log": log}

# If you want to use constants with your script, add them here
config["model_channels"] = "rybg"
config["model_type"] = "mae_contrast_supcon_model"
config["update_model"] = False
config["create_csv"] = False
config["gpu"] = -1

with open(
    os.path.join(
        "models",
        config["model_channels"],
        config["model_type"],
        "model_config.yaml",
    ),
    "r",
) as config_buffer:
    model_config_file = yaml.safe_load(config_buffer)

classifier_paths = None

classifier_paths = None
if "classifier_paths" in model_config_file:
    classifier_paths = model_config_file["classifier_paths"]
encoder_path = model_config_file["encoder_path"]
needs_update = config["update_model"]
for curr_classifier in classifier_paths:
    needs_update = needs_update or not os.path.isfile(curr_classifier)
needs_update = needs_update or not os.path.isfile(encoder_path)


if needs_update:
    config["log"].info("- Downloading models...")
    with open("models_urls.yaml", "r") as urls_file:
        url_info = yaml.safe_load(urls_file)

        for index, curr_url_info in enumerate(url_info[config["model_channels"]][config["model_type"]]["classifiers"]):
            if curr_url_info.startswith("s3://"):
                try:
                    s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))
                    urlcomponents = urlparse(curr_url_info)
                    s3.download_file(urlcomponents.netloc, urlcomponents.path[1:], classifier_paths[index])
                    config["log"].info("  - " + classifier_paths[index] + " updated.")
                except ClientError as e:
                    config["log"].warning("  - " + classifier_paths[index] + " s3 url " + curr_url_info + " not working.")
            else:
                response = requests.get(curr_url_info)
                if response.status_code == 200:
                    with open(classifier_paths[index], "wb") as downloaded_file:
                        downloaded_file.write(response.content)
                    config["log"].info("  - " + classifier_paths[index] + " updated.")
                else:
                    config["log"].warning("  - " + classifier_paths[index] + " url " + curr_url_info + " not found.")

        curr_url_info = url_info[config["model_channels"]][config["model_type"]]["encoder"]
        if curr_url_info.startswith("s3://"):
            try:
                s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))
                urlcomponents = urlparse(curr_url_info)
                s3.download_file(urlcomponents.netloc, urlcomponents.path[1:], encoder_path)
                config["log"].info("  - " + encoder_path + " updated.")
            except ClientError as e:
                config["log"].warning("  - " + encoder_path + " s3 url " + curr_url_info + " not working.")
        else:
            response = requests.get(curr_url_info)
            if response.status_code == 200:
                with open(encoder_path, "wb") as downloaded_file:
                    downloaded_file.write(response.content)
                config["log"].info("  - " + encoder_path + " updated.")
            else:
                config["log"].warning("  - " + encoder_path + " url " + curr_url_info + " not found.")


