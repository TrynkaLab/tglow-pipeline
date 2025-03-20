#!/usr/bin/env python 
import logging
import os
import requests
import yaml
import boto3
import argparse
from botocore import UNSIGNED
from botocore.config import Config
from botocore.exceptions import ClientError
from urllib.parse import urlparse

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

if __name__ == "__main__":
    
    # CLI 
    parser = argparse.ArgumentParser(description="Download subcell model")
    parser.add_argument('-i','--input', help='Path to .yaml config for model', required=True)
    parser.add_argument('--model_urls', help='Path to file with model URLs', required=True)
    parser.add_argument('--model_channels', help='name of the model', required=True)
    parser.add_argument('--model_type', help='name of the model', required=True)
    parser.add_argument('-o','--output', help='Output prefix', default="./")
    args = parser.parse_args()
  
    # args.input = "/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/results/models/rybg/mae_contrast_supcon_model/model_config.yaml"
    # args.model_urls = "/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/results/test_subcell/models_urls.yaml"
    # args.output = "./dave"
    # args.model_channels = "rybg"
    # args.model_type = "mae_contrast_supcon_model"
    # This is the general configuration variable. We are going to use the special key "log" in the dictionary to use the log in our code
    config = {"log": log}

    if not os.path.exists(args.output):
        os.makedirs(f"{args.output}/models/{args.model_channels}/{args.model_type}", exist_ok=True)

    # Read config file
    with open(args.input) as config_buffer:
        model_config_file = yaml.safe_load(config_buffer)

    classifier_paths = model_config_file['classifier_paths']
    encoder_path = model_config_file['encoder_path']

    config["log"].info("- Downloading models...")
    with open(args.model_urls) as urls_file:
        url_info = yaml.safe_load(urls_file)

        for index, curr_url_info in enumerate(url_info[args.model_channels][args.model_type]["classifiers"]):
            if curr_url_info.startswith("s3://"):
                try:
                    s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))
                    urlcomponents = urlparse(curr_url_info)
                    s3.download_file(urlcomponents.netloc, urlcomponents.path[1:], f"{args.output}/{classifier_paths[index]}")
                    config["log"].info("  - " + classifier_paths[index] + " updated.")
                except ClientError as e:
                    config["log"].warning("  - " + classifier_paths[index] + " s3 url " + curr_url_info + " not working.")
            else:
                response = requests.get(curr_url_info)
                if response.status_code == 200:
                    with open(f"{args.output}/{classifier_paths[index]}", "wb") as downloaded_file:
                        downloaded_file.write(response.content)
                    config["log"].info("  - " + classifier_paths[index] + " updated.")
                else:
                    config["log"].warning("  - " + classifier_paths[index] + " url " + curr_url_info + " not found.")

        curr_url_info = url_info[args.model_channels][args.model_type]["encoder"]
        if curr_url_info.startswith("s3://"):
            try:
                s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))
                urlcomponents = urlparse(curr_url_info)
                s3.download_file(urlcomponents.netloc, urlcomponents.path[1:], f"{args.output}/{encoder_path}")
                config["log"].info("  - " + encoder_path + " updated.")
            except ClientError as e:
                config["log"].warning("  - " + encoder_path + " s3 url " + curr_url_info + " not working.")
        else:
            response = requests.get(curr_url_info)
            if response.status_code == 200:
                with open(f"{args.output}/{encoder_path}", "wb") as downloaded_file:
                    downloaded_file.write(response.content)
                config["log"].info("  - " + encoder_path + " updated.")
            else:
                config["log"].warning("  - " + encoder_path + " url " + curr_url_info + " not found.")


