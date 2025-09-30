#!/usr/bin/env python 

import os
import logging
import argparse
import tifffile
import numpy as np
import pandas as pd
import torch.nn.functional as F
import torch
import yaml
import copy
from tqdm import tqdm

from tglow.io.image_query import ImageQuery
from tglow.io.tglow_io import AICSImageReader
from skimage.transform import rescale

from skimage import measure
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union

from torch import nn, Tensor
from skimage.io import imsave
from torchvision.utils import make_grid
from transformers.modeling_outputs import BaseModelOutput
from transformers.models.vit.configuration_vit import ViTConfig
from transformers.models.vit.modeling_vit import (
    BaseModelOutputWithPooling,
    ViTAttention,
    ViTEmbeddings,
    ViTIntermediate,
    ViTOutput,
    ViTPatchEmbeddings,
    ViTPooler,
    ViTPreTrainedModel,
    ViTSdpaAttention,
)

os.environ["DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"

def get_crops(img, msk, outdim=640, dont_mask=False, scale_factor=None):
    regions = measure.regionprops(msk)
    
    crops = {}

    padding=outdim/2

    #for i, region in enumerate(regions):
    for region in regions:
        z, y, x = region.centroid

        id = f"{int(y / scale_factor)}:{int(x / scale_factor)}"
    
        crop = np.zeros((img.shape[0], 640, 640))
    
        minr = int(max(0, y - padding))
        minc = int(max(0, x - padding))
        
        if minr == 0:
            padding_y = int(y)
        else:
            padding_y = padding
            
        if minc == 0:
            padding_x = int(x)
        else:
            padding_x=padding
        
        maxr = int(min(img.shape[1], y + padding_y))
        maxc = int(min(img.shape[2], x + padding_x))
        
        if maxr == img.shape[1]:
            minr == y - (maxr-y)
        
        if maxc == img.shape[1]:
            minc == x - (maxr-x)
        
        curd = img[:,minr:maxr, minc:maxc]

        # Mask the area outside the cell
        if not dont_mask:
            curd = curd * (msk[:,minr:maxr, minc:maxc] == region.label)
        
        large_center = np.array(crop.shape) // 2
        small_center = np.array(curd.shape) // 2
        
        start = large_center - small_center
        end = start + curd.shape
        
        crop[start[0]:end[0], start[1]:end[1], start[2]:end[2]] = curd
                
        # Crop the image
        #crops.append(crop)
        crops[id] = crop
        
    return(crops)


CLASS2NAME = {
    0: "Actin filaments",
    1: "Aggresome",
    2: "Cell Junctions",
    3: "Centriolar satellite",
    4: "Centrosome",
    5: "Cytokinetic bridge",
    6: "Cytoplasmic bodies",
    7: "Cytosol",
    8: "Endoplasmic reticulum",
    9: "Endosomes",
    10: "Focal adhesion sites",
    11: "Golgi apparatus",
    12: "Intermediate filaments",
    13: "Lipid droplets",
    14: "Lysosomes",
    15: "Microtubules",
    16: "Midbody",
    17: "Mitochondria",
    18: "Mitotic chromosome",
    19: "Mitotic spindle",
    20: "Nuclear bodies",
    21: "Nuclear membrane",
    22: "Nuclear speckles",
    23: "Nucleoli",
    24: "Nucleoli fibrillar center",
    25: "Nucleoli rim",
    26: "Nucleoplasm",
    27: "Peroxisomes",
    28: "Plasma membrane",
    29: "Vesicles",
    30: "Negative",
}

CLASS2COLOR = {
    0: "#ffeb3b",
    1: "#76ff03",
    2: "#ff6d00",
    3: "#eb30c1",
    4: "#faadd4",
    5: "#795548",
    6: "#64ffda",
    7: "#00e676",
    8: "#03a9f4",
    9: "#4caf50",
    10: "#ffc107",
    11: "#00bcd4",
    12: "#cddc39",
    13: "#212121",
    14: "#8bc34a",
    15: "#ff9800",
    16: "#ae8c08",
    17: "#ffff00",
    18: "#31b61f",
    19: "#9e9e9e",
    20: "#2196f3",
    21: "#e91e63",
    22: "#3f51b5",
    23: "#9c27b0",
    24: "#673ab7",
    25: "#d3a50b",
    26: "#f44336",
    27: "#009688",
    28: "#ff9e80",
    29: "#242e4b",
    30: "#000000",
}


def min_max_standardize(im):
    min_val = torch.amin(im, dim=(1, 2, 3), keepdims=True)
    max_val = torch.amax(im, dim=(1, 2, 3), keepdims=True)

    im = (im - min_val) / (max_val - min_val + 1e-6)
    return im


def save_attention_map(attn, input_shape, output_path):
    attn = F.interpolate(attn, size=input_shape, mode="bilinear", align_corners=False)
    attn = make_grid(
        attn.permute(1, 0, 2, 3),
        normalize=True,
        nrow=attn.shape[1],
        padding=0,
        scale_each=True,
    )
    attn = (attn.permute(1, 2, 0).cpu().numpy() * 255).astype(np.uint8)
    imsave(output_path + "_attention_map.png", attn)


@torch.no_grad()
def run_model(model, cell_crop, device, output_path=None):
    #cell_crop = np.stack(cell_crop, axis=1)
    cell_crop = torch.from_numpy(cell_crop).float().to(device)
    cell_crop = min_max_standardize(cell_crop)

    output = model(cell_crop)

    probabilities = output.probabilities[0].cpu().numpy()
    embedding = output.pool_op[0].cpu().numpy()

    if output_path is not None:
        np.save(output_path + "_embedding.npy", embedding)
        np.save(output_path + "_probabilities.npy", probabilities)
        save_attention_map(
            output.pool_attn, (cell_crop.shape[2], cell_crop.shape[3]), output_path
        )
    return np.array(embedding), np.array(probabilities)


@dataclass
class ViTPoolModelOutput:
    attentions: Tuple[torch.FloatTensor] = None
    last_hidden_state: torch.FloatTensor = None

    pool_op: torch.FloatTensor = None
    pool_attn: torch.FloatTensor = None

    probabilities: torch.FloatTensor = None


class GatedAttentionPooler(nn.Module):
    def __init__(
        self, dim: int, int_dim: int = 512, num_heads: int = 1, out_dim: int = None
    ):
        super().__init__()

        self.num_heads = num_heads

        self.attention_v = nn.Sequential(nn.Linear(dim, int_dim), nn.Tanh())
        self.attention_u = nn.Sequential(nn.Linear(dim, int_dim), nn.GELU())
        self.attention = nn.Linear(int_dim, num_heads)

        self.softmax = nn.Softmax(dim=-1)

        if out_dim is None:
            self.out_dim = dim * num_heads
            self.out_proj = nn.Identity()
        else:
            self.out_dim = out_dim
            self.out_proj = nn.Linear(dim * num_heads, out_dim)

    def forward(self, x: torch.Tensor) -> Tuple[Tensor, Tensor]:
        v = self.attention_v(x)
        u = self.attention_u(x)

        attn = self.attention(v * u).permute(0, 2, 1)
        attn = self.softmax(attn)

        x = torch.bmm(attn, x)
        x = x.view(x.shape[0], -1)

        x = self.out_proj(x)
        return x, attn


class ViTLayer(nn.Module):
    """This corresponds to the Block class in the timm implementation."""

    def __init__(self, config: ViTConfig, sdpa_attn=False) -> None:
        super().__init__()
        self.chunk_size_feed_forward = config.chunk_size_feed_forward
        self.seq_len_dim = 1
        self.attention = (
            ViTAttention(config) if not sdpa_attn else ViTSdpaAttention(config)
        )
        self.intermediate = ViTIntermediate(config)
        self.output = ViTOutput(config)
        self.layernorm_before = nn.LayerNorm(
            config.hidden_size, eps=config.layer_norm_eps
        )
        self.layernorm_after = nn.LayerNorm(
            config.hidden_size, eps=config.layer_norm_eps
        )

    def forward(
        self,
        hidden_states: torch.Tensor,
        head_mask: Optional[torch.Tensor] = None,
        output_attentions: bool = False,
    ) -> Union[Tuple[torch.Tensor, torch.Tensor], Tuple[torch.Tensor]]:
        self_attention_outputs = self.attention(
            self.layernorm_before(
                hidden_states
            ),  # in ViT, layernorm is applied before self-attention
            head_mask,
            output_attentions=output_attentions,
        )
        attention_output = self_attention_outputs[0]
        outputs = self_attention_outputs[
            1:
        ]  # add self attentions if we output attention weights

        # first residual connection
        hidden_states = attention_output + hidden_states

        # in ViT, layernorm is also applied after self-attention
        layer_output = self.layernorm_after(hidden_states)
        layer_output = self.intermediate(layer_output)

        # second residual connection is done here
        layer_output = self.output(layer_output, hidden_states)

        outputs = (layer_output,) + outputs

        return outputs


class ViTEncoder(nn.Module):
    def __init__(self, config: ViTConfig) -> None:
        super().__init__()
        self.config = config
        # self.layer = nn.ModuleList(
        #     [ViTLayer(config) for _ in range(config.num_hidden_layers)]
        # )
        layer = []
        for i in range(config.num_hidden_layers):
            if i == config.num_hidden_layers - 1:
                layer.append(ViTLayer(config, sdpa_attn=False))
            else:
                layer.append(ViTLayer(config, sdpa_attn=True))
        self.layer = nn.ModuleList(layer)
        self.gradient_checkpointing = False

    def forward(
        self,
        hidden_states: torch.Tensor,
        head_mask: Optional[torch.Tensor] = None,
        output_attentions: bool = False,
        output_hidden_states: bool = False,
        return_dict: bool = True,
    ) -> Union[tuple, BaseModelOutput]:
        all_hidden_states = () if output_hidden_states else None
        all_self_attentions = () if output_attentions else None

        for i, layer_module in enumerate(self.layer):
            if output_hidden_states:
                all_hidden_states = all_hidden_states + (hidden_states,)

            layer_head_mask = head_mask[i] if head_mask is not None else None

            if self.gradient_checkpointing and self.training:
                layer_outputs = self._gradient_checkpointing_func(
                    layer_module.__call__,
                    hidden_states,
                    layer_head_mask,
                    output_attentions,
                )
            else:
                layer_outputs = layer_module(
                    hidden_states, layer_head_mask, output_attentions
                )

            hidden_states = layer_outputs[0]

            if output_attentions:
                all_self_attentions = all_self_attentions + (layer_outputs[1],)

        if output_hidden_states:
            all_hidden_states = all_hidden_states + (hidden_states,)

        if not return_dict:
            return tuple(
                v
                for v in [hidden_states, all_hidden_states, all_self_attentions]
                if v is not None
            )
        return BaseModelOutput(
            last_hidden_state=hidden_states,
            hidden_states=all_hidden_states,
            attentions=all_self_attentions,
        )


class ViTInferenceModel(ViTPreTrainedModel):
    def __init__(
        self,
        config: ViTConfig,
        add_pooling_layer: bool = True,
        use_mask_token: bool = False,
    ):
        super().__init__(config)
        self.config = config

        self.embeddings = ViTEmbeddings(config, use_mask_token=use_mask_token)
        self.encoder = ViTEncoder(config)

        self.layernorm = nn.LayerNorm(config.hidden_size, eps=config.layer_norm_eps)
        self.pooler = ViTPooler(config) if add_pooling_layer else None

        # Initialize weights and apply final processing
        self.post_init()

    def get_input_embeddings(self) -> ViTPatchEmbeddings:
        return self.embeddings.patch_embeddings

    def _prune_heads(self, heads_to_prune: Dict[int, List[int]]) -> None:
        """
        Prunes heads of the model. heads_to_prune: dict of {layer_num: list of heads to prune in this layer} See base
        class PreTrainedModel
        """
        for layer, heads in heads_to_prune.items():
            self.encoder.layer[layer].attention.prune_heads(heads)

    def forward(
        self,
        pixel_values: Optional[torch.Tensor] = None,
        bool_masked_pos: Optional[torch.BoolTensor] = None,
        head_mask: Optional[torch.Tensor] = None,
        output_attentions: Optional[bool] = None,
        output_hidden_states: Optional[bool] = None,
        interpolate_pos_encoding: Optional[bool] = None,
        return_dict: Optional[bool] = None,
    ) -> Union[Tuple, BaseModelOutputWithPooling]:
        r"""
        bool_masked_pos (`torch.BoolTensor` of shape `(batch_size, num_patches)`, *optional*):
            Boolean masked positions. Indicates which patches are masked (1) and which aren't (0).
        """
        output_attentions = (
            output_attentions
            if output_attentions is not None
            else self.config.output_attentions
        )
        output_hidden_states = (
            output_hidden_states
            if output_hidden_states is not None
            else self.config.output_hidden_states
        )
        return_dict = (
            return_dict if return_dict is not None else self.config.use_return_dict
        )

        if pixel_values is None:
            raise ValueError("You have to specify pixel_values")

        # Prepare head mask if needed
        # 1.0 in head_mask indicate we keep the head
        # attention_probs has shape bsz x n_heads x N x N
        # input head_mask has shape [num_heads] or [num_hidden_layers x num_heads]
        # and head_mask is converted to shape [num_hidden_layers x batch x num_heads x seq_length x seq_length]
        head_mask = self.get_head_mask(head_mask, self.config.num_hidden_layers)

        # TODO: maybe have a cleaner way to cast the input (from `ImageProcessor` side?)
        expected_dtype = self.embeddings.patch_embeddings.projection.weight.dtype
        if pixel_values.dtype != expected_dtype:
            pixel_values = pixel_values.to(expected_dtype)

        embedding_output = self.embeddings(
            pixel_values,
            bool_masked_pos=bool_masked_pos,
            interpolate_pos_encoding=interpolate_pos_encoding,
        )

        encoder_outputs = self.encoder(
            embedding_output,
            head_mask=head_mask,
            output_attentions=output_attentions,
            output_hidden_states=output_hidden_states,
            return_dict=return_dict,
        )
        sequence_output = encoder_outputs[0]
        sequence_output = self.layernorm(sequence_output)
        pooled_output = (
            self.pooler(sequence_output) if self.pooler is not None else None
        )

        if not return_dict:
            head_outputs = (
                (sequence_output, pooled_output)
                if pooled_output is not None
                else (sequence_output,)
            )
            return head_outputs + encoder_outputs[1:]

        return BaseModelOutputWithPooling(
            last_hidden_state=sequence_output,
            pooler_output=pooled_output,
            hidden_states=encoder_outputs.hidden_states,
            attentions=encoder_outputs.attentions,
        )


class ViTPoolClassifier(nn.Module):
    def __init__(self, config: Dict):
        super(ViTPoolClassifier, self).__init__()

        self.vit_config = ViTConfig(**config["vit_model"])

        self.encoder = ViTInferenceModel(self.vit_config, add_pooling_layer=False)

        pool_config = config.get("pool_model")
        self.pool_model = GatedAttentionPooler(**pool_config) if pool_config else None

        self.out_dim = (
            self.pool_model.out_dim if self.pool_model else self.vit_config.hidden_size
        )

        self.num_classes = config["num_classes"]
        self.sigmoid = nn.Sigmoid()

    def make_classifier(self):
        return nn.Sequential(
            nn.Linear(self.out_dim, 512),
            nn.ReLU(inplace=True),
            nn.Linear(512, 256),
            nn.ReLU(inplace=True),
            nn.Linear(256, self.num_classes),
        )

    def load_model_dict(
        self,
        encoder_path: str,
        classifier_paths: Union[str, List[str]],
        device="cpu",
    ):
        checkpoint = torch.load(encoder_path, map_location=device, weights_only=True)
        encoder_ckpt = {
            k[len("encoder.") :]: v for k, v in checkpoint.items() if "encoder." in k
        }

        status = self.encoder.load_state_dict(encoder_ckpt)
        print(f"Encoder status: {status}")

        pool_ckpt = {
            k.replace("pool_model.", ""): v
            for k, v in checkpoint.items()
            if "pool_model." in k
        }
        pool_ckpt = {k.replace("1.", "0."): v for k, v in pool_ckpt.items()}
        if pool_ckpt and self.pool_model:
            status = self.pool_model.load_state_dict(pool_ckpt)
            print(f"Pool model status: {status}")
        else:
            print("No pool model found in checkpoint")

        if isinstance(classifier_paths, str):
            classifier_paths = [classifier_paths]

        self.classifiers = nn.ModuleList(
            [self.make_classifier() for _ in range(len(classifier_paths))]
        )
        for i, classifier_path in enumerate(classifier_paths):
            classifier_ckpt = torch.load(classifier_path, map_location=device, weights_only=True)
            classifier_ckpt = {
                k.replace("3.", "2."): v for k, v in classifier_ckpt.items()
            }
            classifier_ckpt = {
                k.replace("6.", "4."): v for k, v in classifier_ckpt.items()
            }
            status = self.classifiers[i].load_state_dict(classifier_ckpt)
            print(f"Classifier {i+1} status: {status}")

    def forward(self, x: torch.Tensor) -> ViTPoolModelOutput:
        b, c, h, w = x.shape
        outputs = self.encoder(x, output_attentions=True, interpolate_pos_encoding=True)

        if self.pool_model:
            pool_op, pool_attn = self.pool_model(outputs.last_hidden_state)
        else:
            pool_op = torch.mean(outputs.last_hidden_state, dim=1)
            pool_attn = None

        probs = torch.stack(
            [self.sigmoid(classifier(pool_op)) for classifier in self.classifiers],
            dim=1,
        )
        probs = torch.mean(probs, dim=1)

        h_feat = h // self.vit_config.patch_size
        w_feat = w // self.vit_config.patch_size

        attentions = outputs.attentions[-1][:, :, 0, 1:].reshape(
            b, self.vit_config.num_attention_heads, h_feat, w_feat
        )

        pool_attn = pool_attn[:, :, 1:].reshape(
            b, self.pool_model.num_heads, h_feat, w_feat
        )

        return ViTPoolModelOutput(
            last_hidden_state=outputs.last_hidden_state,
            attentions=attentions,
            pool_op=pool_op,
            pool_attn=pool_attn,
            probabilities=probs,
        )




 # Logging


logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Main loop
if __name__ == "__main__":
                  
    parser = argparse.ArgumentParser(description="Run subcell inference on a plate/row/col/field.ome.tiff and associated masks ")

    # # Common argulents
    parser.add_argument('-i','--input', help='Base dir to raw input', required=True)
    parser.add_argument('-o','--output', help='Output prefix', default="./")
    parser.add_argument('-p','--plate', help='Subfolder in raw dir to process', required=True)
    parser.add_argument('-w','--well', help='Subfolder in raw dir to process', required=True)
    parser.add_argument('--model', help='Model dir', default="./")
    parser.add_argument('--ch_nucl', help='Nucleus channel', required=True)
    parser.add_argument('--ch_er', help='ER channel', required=True)
    parser.add_argument('--ch_tub', help='Tublin channel', required=True)
    parser.add_argument('--channels', help='Channels to infer localization for. <name>=<channel> [<name>=<channel>]', required=True, nargs="+")
    parser.add_argument('--scale_factor', help="The scale factor to get to 80nm per pixel", default=149/80)
    parser.add_argument('--gpu', help='Use gpu, uses the first device', default=False, action="store_true")
    parser.add_argument('--mask_pattern', help="The pattern to discover masks, defaults to match run_cellpose.py cell. <pattern>.tiff", default="*_cell_mask_*_cp_masks.tiff")
    parser.add_argument('--save_embeddings', help='Save the embeddings', default=False, action="store_true")
    parser.add_argument('--dont_mask', help='Do not mask the range outside the cell in the crop', default=False, action="store_true")
    args = parser.parse_args()
    
    # # #cd /lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/results/test_subcell
    
    # parser = argparse.ArgumentParser(description="Run subcell inference on a plate/row/col/field.ome.tiff and associated masks ")
    # args = parser.parse_args()
    # args.output="/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/results/test_subcell2"
    # args.input="/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/results/processed_images"
    # args.well = "C10"
    # args.plate = "mo13_240518_TGlow_drugperturb_72h_plate1_cycle1"
    # args.scale_factor=149/80
    # args.model="/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/results/models/rybg/mae_contrast_supcon_model"
    # args.gpu = False
    # args.save_embeddings = False
    # args.mask_pattern = "*_cell_mask_*_cp_masks.tiff"
    # args.ch_nucl = 7
    # args.ch_er = 5
    # args.ch_tub = 8
    # args.channels = ['cd25_ki67=0', 'ki67=10', 'cd25=11', 'mito=1', 'actin=3']
    # args.dont_mask = False
    
    # Covert channels to int
    tmp_channels = {}
    for channel in args.channels:
        key, value = channel.split("=")
        tmp_channels[key] = int(value)
    args.channels = tmp_channels
    log.info(f"Running SubCell localization for channels: {args.channels}")
    
    args.ch_nucl = int(args.ch_nucl)
    args.ch_er = int(args.ch_er)
    args.ch_tub = int(args.ch_tub)
    args.scale_factor = float(args.scale_factor)
    
    #--------------------------------------------------------------------------------------
    # Setup the image readers
    img_reader = AICSImageReader(args.input, plates_filter=args.plate)
    msk_reader = AICSImageReader(args.input, plates_filter=args.plate, pattern=args.mask_pattern)
    
    #--------------------------------------------------------------------------------------
    # Read the model config file
    with open(f"{args.model}/model_config.yaml") as config_buffer:
        model_config_file = yaml.safe_load(config_buffer)

    classifier_paths = [f"{args.model}/{os.path.basename(x)}" for x in model_config_file["classifier_paths"]]
    encoder_path = f"{args.model}/{os.path.basename(model_config_file['encoder_path'])}"

    #--------------------------------------------------------------------------------------
    # Setup the device to use
    num_devices = torch.cuda.device_count()

    # Print information about each device
    for i in range(num_devices):
        device = torch.cuda.get_device_properties(i)
        log.debug(f"Device {i}: {device.name}")
        log.debug(f"  Compute Capability: {device.major}.{device.minor}")
        log.debug(f"  Total Memory: {device.total_memory / 1024**3:.2f} GB")

    if torch.cuda.is_available() and args.gpu :
        log.info("Using the first GPU")
        device = torch.device("cuda:" + str(0))
    else:
        log.warning("CUDA not available. Using CPU.")
        device = torch.device("cpu")

    #--------------------------------------------------------------------------------------
    # Setup the model object
    model_config = model_config_file.get("model_config")
    model_config['vit_model']['num_channels'] = 4# 3 + len(args.channels)
    model = ViTPoolClassifier(model_config)
    model.load_model_dict(encoder_path, classifier_paths)
    model.eval()
      
    model.to(device)
    #--------------------------------------------------------------------------------------
    # Setup the output DFs
    #out_header = ["ObjectNumber", "ImageNumber", "Parent_cell", "Location_x", "Location_y"]
    #out_header = ["ObjectNumber", "ImageNumber", "Parent_cell"]
    out_header = []
    out_header.extend([f"Localization_{CLASS2NAME[x].title().replace(' ', '')}" for x in range(0, len(CLASS2NAME))])
    
    if args.save_embeddings:
        out_header.extend([f"Embedding_{x}" for x in range(0, 1536)])
  
    #out_df = pd.DataFrame(columns=out_header)
  
    meta_header = ["ImageNumber"] + [f"Metadata_SubCell_{x}" for x in ['plate', 'well', 'row', 'col', 'field']]
    meta_df = pd.DataFrame(columns=meta_header)
    cell_df = pd.DataFrame(columns=["ObjectNumber", "ImageNumber", "Location_x", "Location_y"])
    locl_df = {}

    # Start the extraction loop
    row, col = ImageQuery.well_id_to_index(args.well)
    
    log.info("-----------------------------------------------------------")
    args.output = f"{args.output}/{args.plate}/{ImageQuery.ID_TO_ROW[str(row)]}/{col}/"
    
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        log.info(f"Folder created: {args.output}")
        log.info("-----------------------------------------------------------")

    
    pb = tqdm(range(len(img_reader.fields[args.plate]) * len(args.channels)))
    
    img_id = 0
    cell_id = 0
    for field in img_reader.fields[args.plate][str(row)][str(col)]:
        iq = ImageQuery(args.plate, row, col, field)
        
        pb.set_description(f"Reading field {field}")
        # Assume images are already MP, so drop the Z plate
        img = img_reader.read_stack(iq)[:,0,:,:]
        msk = msk_reader.read_stack(iq)[:,0,:,:]
        
        pb.set_description(f"Scaling field {field}")
        # Rescale to the right factor
        img = rescale(img, args.scale_factor, channel_axis=0, order=0)   
        msk = rescale(msk, args.scale_factor, channel_axis=0, order=0)       
        crops = get_crops(img, msk, dont_mask=args.dont_mask, scale_factor=args.scale_factor)    
        
        
        
        # Add to the image level df
        meta_df.loc[img_id] = [img_id+1, iq.plate, iq.get_well_id(), iq.get_row_letter(), iq.col, iq.field]

        # Get the list of crop keys, so order is maintained
        cell_keys = [x for x in crops.keys()]
    
        j = cell_id
        for key in cell_keys:
            x, y = [int(x) for x in key.split(":")]
            cell_df.loc[j] = [j+1, img_id+1, x, y]
            j += 1
        
        # Loop over the channels        
        for channel_name in args.channels.keys():
            
            if channel_name not in locl_df.keys():
                cur_header = copy.deepcopy(out_header)
                cur_header = [f"{x}_{channel_name}" for x in cur_header]  
                locl_df[channel_name] = pd.DataFrame(columns=cur_header)
                
            cur_out_df = locl_df[channel_name]
            channel = args.channels[channel_name]
            
            pb.set_description(f"Processing {channel_name} field {field}")
            # Loop over images here
            j = cell_id
            for key in cell_keys:
                
                pb.set_description(f"Processing {channel_name} field {field} cell {key}")
                #cell_data=crops[key][np.newaxis,([args.ch_nucl, args.ch_er, args.ch_tub] + args.channels), :,:]
                cell_data=crops[key][np.newaxis,(args.ch_tub, args.ch_er, args.ch_nucl,  channel), :,:]
                embedding, probabilities = run_model(
                    model,
                    cell_data,
                    device
                )                
                if args.save_embeddings:
                    cur_out_df.loc[j] = probabilities.tolist() + embedding.tolist()
                else:
                    cur_out_df.loc[j] = probabilities.tolist()
                j += 1   
            pb.update(1)
            
        cell_id = j
        img_id += 1
    pb.close()
    
    log.info("Done extracting, saving")
    meta_df.to_csv(f'{args.output}/SubCell_Image.txt', sep="\t", index=False)    
    
    # Bind the cell dfs together
    df_out  = pd.concat([cell_df] + [x for x in locl_df.values()], axis=1)
    df_out.to_csv(f'{args.output}/SubCell_cell.txt', sep="\t", index=False)    

    # Write all items in the Namespace to a file
    with open(f'{args.output}/SubCell_Experiment.txt', "w") as f:
        for key, value in vars(args).items():
            f.write(f"{key}: {value}\n")
            
        
        
    # for channel_name in args.channels.keys():
    #     channel = args.channels[channel_name]
    #     locl_df[channel_name].to_csv(f'{args.output}/SubCell_cellCh{channel}.txt', sep="\t", index=False)   
        
    # curr_probs_l = probabilities.tolist()
    # max_location_class = curr_probs_l.index(max(curr_probs_l))
    # max_location_name = CLASS2NAME[max_location_class]
    # max_3_location_classes = sorted(
    #     range(len(curr_probs_l)), key=lambda sub: curr_probs_l[sub]
    # )[-3:]
    # max_3_location_classes.reverse()
    # max_3_location_names = (
    #     CLASS2NAME[max_3_location_classes[0]]
    #     + ","
    #     + CLASS2NAME[max_3_location_classes[1]]
    #     + ","
    #     + CLASS2NAME[max_3_location_classes[2]]
    # )
    # cell_crop = np.stack(cell_data, axis=1)
    # cell_crop = torch.from_numpy(cell_crop).float().to(device)
    # cell_crop = min_max_standardize(cell_crop)

    # output = model(cell_crop)
    # save_attention_map(output.attentions, (cell_crop.shape[2], cell_crop.shape[3]), "testing_1_12_attn")
    
    


