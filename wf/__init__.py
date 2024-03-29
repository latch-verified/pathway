import csv
import functools
import json
import os
import re
import subprocess
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union
from xml.etree import ElementTree

import imagesize
from flytekit import LaunchPlan
from jinja2 import Template
from latch import medium_task, message, workflow
from latch.types import (LatchAuthor, LatchDir, LatchFile, LatchMetadata,
                         LatchParameter)

print = functools.partial(print, flush=True)


class Species(Enum):
    HUMAN = "Homo sapiens"
    MOUSE = "Mus musculus"
    YEAST = "Saccharomyces cerevisiae"
    RAT = "Rattus norvegicusz"


species_2_msig = {
    Species.HUMAN: "Homo sapiens",
    Species.MOUSE: "Mus musculus",
    Species.YEAST: "Saccharomyces cerevisiae",
    Species.RAT: "Rattus norvegicus",
}

species_2_go = {
    Species.HUMAN: "org.Hs.eg.db",
    Species.MOUSE: "org.Mm.eg.db",
    Species.YEAST: "org.Sc.sgd.db",
    Species.RAT: "org.Rn.eg.db",
}

species_2_kegg = {
    Species.HUMAN: "hsa",
    Species.MOUSE: "mmu",
    Species.YEAST: "sce",
    Species.RAT: "rno",
}


def run_and_capture_output(command: List[str]) -> Tuple[int, str]:
    captured_stdout = []

    with subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=10,
        universal_newlines=True,
    ) as process:
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="")
            captured_stdout.append(line)
        process.wait()
        returncode = process.returncode

    return returncode, "".join(captured_stdout)


def slugify(value: str, keep_spaces: bool = False) -> str:
    value = value.replace("/", "_")
    if not keep_spaces:
        value = value.replace(" ", "_")
    return value


GENESETS_PATH = "/root/tempres/genesets.txt"


def pathway_to_genesets_mapping(
    relevant_pathway_ids: Set[str],
) -> Dict[str, Tuple[str, str]]:
    path = Path(GENESETS_PATH).resolve()
    pathway_ids = []
    entrez_ids = []
    gene_names = []
    with path.open("r") as f:
        current = None
        split = None
        for line in f.readlines():
            if line.startswith("PATHWAYIDS"):
                current, split = pathway_ids, False
            elif line.startswith("ENTREZIDS"):
                current, split = entrez_ids, True
            elif line.startswith("NAMES"):
                current, split = gene_names, True
            elif current is not None:
                data = line.strip()
                if split:
                    data = data.split(" ")
                current.append(data)
    mapping = dict(zip(pathway_ids, zip(entrez_ids, gene_names)))
    return {k: v for k, v in mapping.items() if k in relevant_pathway_ids}


def js_injectable(data):
    # Replace quotes with escaped quotes to prevent premature string termination
    # when injecting this data into JavaScript code – see usage of the method
    # injection() in /report/src/index.js
    return json.dumps(data).replace('"', '\\"')


def make_message_pattern(flag: str) -> re.Pattern:
    return re.compile(f"__LATCH_{flag}_START__(.*?)__LATCH_{flag}_END__", re.DOTALL)


ERROR_PATTERN = make_message_pattern("ERROR")
WARNING_PATTERN = make_message_pattern("WARNING")
PATHVIEW_IMAGE_WIDTH = 800


def emit_messages(pattern: re.Pattern, message_type: str, stdout_logs: str) -> None:
    for item in re.findall(pattern, stdout_logs):
        title = f"Pathway enrichment analysis {message_type}"
        message(message_type, {"title": title, "body": item})


@dataclass
class Pathway:
    kegg_id: str
    name: str
    p_value: str
    p_adjusted_value: str
    enrichment_score: str
    normalized_enrichment_score: str
    gene_set_size: str
    leading_edge_size: str
    core_enriched_genes: List[str]
    core_entrez_ids: List[str]


def pathway_to_dict(pathway: Pathway) -> Dict[str, Union[str, List[str]]]:
    return {
        "pathwayId": pathway.kegg_id,
        "pathwayName": pathway.name,
        "pValue": pathway.p_value,
        "pAdjusted": pathway.p_adjusted_value,
        "enrichmentScore": pathway.enrichment_score,
        "normalizedEnrichmentScore": pathway.normalized_enrichment_score,
        "geneSetSize": pathway.gene_set_size,
        "leadingEdgeSize": pathway.leading_edge_size,
        "coreEnrichedGenes": pathway.core_enriched_genes,
    }


KEGG_TABLE_PATH = "./res/KEGG/table.csv"


def parse_kegg_pathway_table() -> List[Pathway]:
    pathways: List[Pathway] = []

    with Path(KEGG_TABLE_PATH).open("r") as f:
        for row in csv.DictReader(f):
            gene_names = row["coreEnrichedGenes"].split(" ")
            pathways.append(
                Pathway(
                    kegg_id=row["ID"],
                    name=row["Description"],
                    p_value=f'{row["pvalue"]:.3}',
                    p_adjusted_value=f'{row["p.adjust"]:.3}',
                    enrichment_score=f'{row["enrichmentScore"]:.6}',
                    normalized_enrichment_score=f'{row["NES"]:.6}',
                    gene_set_size=row["setSize"],
                    leading_edge_size=str(len(gene_names)),
                    core_enriched_genes=gene_names,
                    core_entrez_ids=row["coreEntrezIDs"].split("/"),
                )
            )

    return pathways


def parse_contrast_csv(contrast_csv: LatchFile):
    contrast_csv_path = Path(contrast_csv).resolve()
    with open(contrast_csv_path, "r") as f:
        # Gene column name is an empty string in current differential gene
        # expression workflow's contrast csv outputs
        return {
            row[""]: [
                f"{row['log2FoldChange']:.4}",
                f"{row['pvalue']:.3}",
                f"{row['padj']:.3}",
            ]
            for row in csv.DictReader(f)
        }


def local_results_path(name: Optional[str] = None) -> Path:
    path = "./res"
    if name is not None:
        path = f"{path}/{name}"
    return Path(path)


def remote_results_path(
    output_location: LatchDir,
    report_name: str,
    name: Optional[str] = None,
) -> str:
    base_path = output_location.remote_path
    assert base_path is not None
    path = os.path.join(base_path, slugify(report_name, keep_spaces=True))
    if name is not None:
        path = os.path.join(path, slugify(name, keep_spaces=True))
    return path


TEMPLATE_REPORT_PATH = "./template.html"


def create_report(report_name: str, kegg_data: Optional[Dict[str, Any]]) -> Path:
    report_path = local_results_path("Report.html")
    with report_path.open("w") as freport:
        with Path(TEMPLATE_REPORT_PATH).open("r") as ftemplate:
            template = Template(ftemplate.read())

            data: Dict[str, Any] = {"reportName": report_name}
            if kegg_data is not None:
                data.update(kegg_data)
                data["showKEGG"] = True
            else:
                data["showKEGG"] = False

            html = template.render(data)
            freport.write(html)

    return report_path


def make_gene_annotation_view_rect(
    image_width: float,
    graphics: ElementTree,
) -> Dict[str, float]:

    scale = PATHVIEW_IMAGE_WIDTH / image_width
    attributes = ("x", "y", "width", "height")

    # Want type of "rect"
    if graphics.get("type") == "line":
        return

    x, y, width, height = [float(graphics.get(a)) * scale for a in attributes]

    # Shift x, y coordinates from center to top-left
    return {
        "x": x - width / 2,
        "y": y - height / 2,
        "width": width,
        "height": height,
    }


def parse_gene_groups(
    pathway: Pathway,
    contrast_genes: Set[str],
    pathview_path: Path,
    species: Species,
) -> List[Dict[str, Any]]:
    gene_groups = []
    entrez_id_to_gene = dict(zip(pathway.core_entrez_ids, pathway.core_enriched_genes))

    pathview_image_width = imagesize.get(str(pathview_path))[0]

    # Parse pathview XML to create hoverable gene annotations
    xml_path = Path(str(pathview_path).replace(".pathview.png", ".xml"))
    with xml_path.open("r") as f:
        xmldoc = ElementTree.parse(f)

    for entry in xmldoc.iter("entry"):
        if entry.get("type") != "gene":
            continue

        raw_names = entry.get("name")
        if raw_names is None:
            continue

        graphics = entry.find("graphics")
        if graphics is None:
            continue

        genes: Dict[str, bool] = {}
        for name in raw_names.split(" "):
            if not name.startswith(f"{species_2_kegg[species]}:"):
                continue
            entrez_id = name.removeprefix(f"{species_2_kegg[species]}:")
            if entrez_id not in entrez_id_to_gene:
                continue
            genes[entrez_id_to_gene[entrez_id]] = True

        graphics_raw_genes = graphics.get("name")
        if graphics_raw_genes is not None:
            graphics_genes = graphics_raw_genes.removesuffix("...").split(", ")
            for new_gene in graphics_genes:
                if new_gene in contrast_genes and new_gene not in genes:
                    genes[new_gene] = False

        if len(genes) == 0:
            continue

        view = make_gene_annotation_view_rect(pathview_image_width, graphics)
        if view:
            gene_groups.append(
                {
                    "view": view,
                    "core": any(genes.values()),
                    "genes": [{"name": k, "core": v} for k, v in genes.items()],
                }
            )

    return gene_groups


def run_rscript(
    contrast_csv: LatchFile, number_of_pathways: int, species: Species
) -> None:
    print("Running go pathway")

    returncode, stdout = run_and_capture_output(
        [
            "Rscript",
            "/root/go_pathway.r",
            str(Path(contrast_csv).resolve()),
            str(number_of_pathways),
            species_2_msig[species],
            species_2_go[species],
            species_2_kegg[species],
        ]
    )

    print("Finished go pathway")
    emit_messages(ERROR_PATTERN, "error", stdout)
    emit_messages(WARNING_PATTERN, "warning", stdout)

    if returncode != 0:
        raise RuntimeError(f"R script failed with exit code '{returncode}'")


@medium_task
def go_pathway(
    contrast_csv: LatchFile,
    report_name: str,
    species: Species,
    output_location: LatchDir,
    number_of_pathways: int,
) -> LatchDir:
    run_rscript(contrast_csv, number_of_pathways, species)

    try:
        pathways = parse_kegg_pathway_table()
        kegg_data = {}
    except FileNotFoundError:
        # Output a report without KEGG plots and without pathway viewers
        pathways = None
        kegg_data = None

    if kegg_data is not None:
        contrast_csv_data = parse_contrast_csv(contrast_csv)
        pathviews: List[Dict[str, str]] = []
        pathway_id_to_gene_groups: Dict[str, List[Dict[str, Any]]] = {}
        contrast_genes = set(contrast_csv_data.keys())

        assert pathways is not None
        for p in pathways:
            pathview_path = Path(f"./{p.kegg_id}.pathview.png")

            pathway_id_to_gene_groups[p.kegg_id] = parse_gene_groups(
                p, contrast_genes, pathview_path, species
            )

            relative_remote_path = f"Pathview/{slugify(p.name)}.png"
            path = local_results_path(relative_remote_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            pathview_path.rename(path)
            pathviews.append({"id": p.kegg_id, "path": f"./{relative_remote_path}"})

        mapping = pathway_to_genesets_mapping(set(x.kegg_id for x in pathways))
        kegg_data["pathwayIdToGeneSets"] = js_injectable(mapping)

        kegg_data["pathwayIdToGeneGroups"] = js_injectable(pathway_id_to_gene_groups)
        kegg_data["pathwayData"] = js_injectable([pathway_to_dict(x) for x in pathways])
        kegg_data["contrastData"] = js_injectable(contrast_csv_data)
        kegg_data["pathviews"] = pathviews

    create_report(report_name, kegg_data)

    path = remote_results_path(output_location, report_name)
    return LatchDir(str(local_results_path()), path)


metadata = LatchMetadata(
    display_name="Pathway Enrichment Analysis",
    documentation="https://www.latch.wiki/bulk-rna-seq-end-to-end#49afbcfd1f9d4ef381644d0e11e44bd2",
    wiki_url="https://www.latch.wiki/bulk-rna-seq-end-to-end#49afbcfd1f9d4ef381644d0e11e44bd2",
    video_tutorial="https://www.loom.com/share/c3848161fd6347cd8b560c5c904116d3",
    author=LatchAuthor(
        name="Latch Verified",
        github="https://github.com/latch-verified",
    ),
    repository="https://github.com/latch-verified/pathway",
    license="MIT",
    parameters={
        "contrast_csv": LatchParameter(
            display_name="Contrast Data",
            description="Contrast CSV from the Differential Gene Expression workflow",
            batch_table_column=True,
        ),
        "report_name": LatchParameter(
            display_name="Report Name",
            description="Name of your report",
            batch_table_column=True,
        ),
        "species": LatchParameter(
            display_name="Species",
            description="Which species to condition pathways on",
        ),
        "number_of_pathways": LatchParameter(
            display_name="Number of top pathways to display",
            description="Pathways are ranked by their enrichment score",
        ),
        "output_location": LatchParameter(
            display_name="Output Location",
            description="Your results folder",
            output=True,
        ),
    },
)


@workflow(metadata)
def gene_ontology_pathway_analysis(
    contrast_csv: LatchFile,
    report_name: str,
    species: Species = Species.HUMAN,
    number_of_pathways: int = 20,
    output_location: LatchDir = LatchDir("latch:///Pathway Analysis/"),
) -> LatchDir:
    """Use differential expression contrast data to calculate the gene ontology and pathways for the most significant genes.

    Under the hood, the workflow uses Gene Set Enrichment Analysis (GSEA) and
    searches for the top 20 pathways using mSigDB, KEGG, and GO pathway
    databases. For more information, see our [wiki]("https://www.latch.wiki/bulk-rna-seq-end-to-end#49afbcfd1f9d4ef381644d0e11e44bd2").
    """
    return go_pathway(
        contrast_csv=contrast_csv,
        report_name=report_name,
        species=species,
        number_of_pathways=number_of_pathways,
        output_location=output_location,
    )


if __name__ == "wf":
    LaunchPlan.create(
        "wf.__init__.gene_ontology_pathway_analysis.Abs vs Gingiva (Foote, 2019)",
        gene_ontology_pathway_analysis,
        default_inputs={
            "contrast_csv": LatchFile(
                "s3://latch-public/welcome/gsea_pathway/Abs vs Gingiva (Condition).csv"
            ),
            "report_name": "Test Data - Abs vs Gingiva (Foote, 2019)",
        },
    )

    LaunchPlan.create(
        "wf.__init__.gene_ontology_pathway_analysis.CoCl2 vs Control (Knyazev, 2021)",
        gene_ontology_pathway_analysis,
        default_inputs={
            "contrast_csv": LatchFile(
                "s3://latch-public/welcome/gsea_pathway/CoCl2 vs Control (Condition).csv"
            ),
            "report_name": "Test Data - CoCl2 vs Control (Knyazev, 2021)",
        },
    )
