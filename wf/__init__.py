import csv
import json
import os
import re
import subprocess
from pathlib import Path
from xml.etree import ElementTree

import imagesize
from jinja2 import Template
from latch import large_task, message, workflow
from latch.types import LatchAuthor, LatchDir, LatchFile, LatchMetadata, LatchParameter

_PATHVIEW_IMAGE_WIDTH = 800


error_pattern = re.compile("__LATCH_ERROR_START__(.*?)__LATCH_ERROR_END__", re.DOTALL)


_WARNING_PATTERN = re.compile(
    "__LATCH_WARNING_START__(.*?)__LATCH_WARNING_END__", re.DOTALL
)


def _capture_output(command: list[str]) -> tuple[int, str]:
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


def _slugify(value: str) -> str:
    return value.replace("/", "_").replace(" ", "_")


def _base_output_path(output_location: LatchDir, report_name: str) -> str:
    output = output_location.remote_path
    assert output is not None
    if output[-1] != "/":
        output += "/"
    output += report_name.replace("/", "_")
    return output


def _pathway_entrez_ids():
    path = Path("/root/tempres/genesets.txt").resolve()
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
    return dict(zip(pathway_ids, zip(entrez_ids, gene_names)))


def _make_injectable(data):
    # Replace quotes with escaped quotes to prevent premature string termination
    # when injecting this data into JavaScript code
    return json.dumps(data).replace('"', '\\"')


@large_task
def go_pathway(
    contrast_csv: LatchFile,
    report_name: str,
    output_location: LatchDir,
    number_of_pathways: int,
) -> tuple[LatchDir, list[LatchFile]]:
    contrast_csv_path = Path(contrast_csv).resolve()

    print("Running go pathway", flush=True)
    returncode, stdout = _capture_output(
        [
            "Rscript",
            "/root/go_pathway.r",
            str(contrast_csv_path),
            str(number_of_pathways),
        ]
    )
    print("Finished go pathway")

    for item in re.findall(error_pattern, stdout):
        message(
            "error",
            {"title": "Gene ontology and pathway analysis ERROR", "body": item},
        )
    for item in re.findall(_WARNING_PATTERN, stdout):
        message(
            "warning",
            {"title": "Gene ontology and pathway analysis warning", "body": item},
        )

    if returncode != 0:
        raise RuntimeError(f"R script failed with exit code '{returncode}'")

    contrast_data = {}

    with open(contrast_csv_path, "r") as f:
        rows = csv.DictReader(f)
        # Gene column name is an empty string (at least in current DESeq2 contrast csv outputs)
        contrast_data = {
            row[""]: [
                f"{row['log2FoldChange']:.4}",
                f"{row['pvalue']:.3}",
                f"{row['padj']:.3}",
            ]
            for row in rows
        }

    contrast_genes = set(contrast_data.keys())
    output = _base_output_path(output_location, report_name)
    pathviews: list[LatchFile] = []
    data = {
        "reportName": report_name,
        "pathwayData": [],
        "pathviews": [],
        "pathwayIdToGeneGroups": {},
        "contrastData": _make_injectable(contrast_data),
    }

    with Path("./res/KEGG/table.csv").open("r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_names = row["coreEnrichedGenes"].split(" ")
            core_entrez_ids = row["coreEntrezIDs"].split("/")
            entrez_id_to_gene = dict(zip(core_entrez_ids, gene_names))
            pathway_id = row["ID"]
            pathway_name = row["Description"]

            data["pathwayData"].append(
                {
                    "pathwayId": pathway_id,
                    "pathwayName": pathway_name,
                    "pValue": f'{row["pvalue"]:.3}',
                    "pAdjusted": f'{row["p.adjust"]:.3}',
                    "enrichmentScore": f'{row["enrichmentScore"]:.6}',
                    "normalizedEnrichmentScore": f'{row["NES"]:.6}',
                    "geneSetSize": row["setSize"],
                    "leadingEdgeSize": str(len(gene_names)),
                    "coreEnrichedGenes": row["coreEnrichedGenes"],
                }
            )

            image_path = f"./{pathway_id}.pathview.png"

            print("Saving pathview image", image_path)

            relative_remote_path = f"Pathview/{_slugify(pathway_name)}.png"
            absolute_remote_path = os.path.join(output, relative_remote_path)
            pathviews.append(LatchFile(image_path, absolute_remote_path))

            # Parse pathview XML to create hoverable gene annotations
            image_width = imagesize.get(image_path)[0]
            xml_path = image_path.replace(".pathview.png", ".xml")
            with open(xml_path, "r") as f:
                xmldoc = ElementTree.parse(f)

            gene_groups = []
            for entry in xmldoc.iter("entry"):
                if entry.get("type") != "gene":
                    continue

                raw_names = entry.get("name")
                if raw_names is None:
                    continue

                graphics = entry.find("graphics")
                if graphics is None:
                    continue

                genes: dict[str, bool] = {}
                for name in raw_names.split(" "):
                    if not name.startswith("hsa:"):
                        continue

                    entrez_id = name.removeprefix("hsa:")
                    if entrez_id not in entrez_id_to_gene:
                        continue

                    gene = entrez_id_to_gene[entrez_id]
                    genes[gene] = True

                if graphics.get("name") is not None:
                    all_genes = graphics.get("name").removesuffix("...").split(", ")
                    for new_gene in all_genes:
                        if new_gene in contrast_genes and new_gene not in genes:
                            genes[new_gene] = False

                if len(genes) == 0:
                    continue

                scale = _PATHVIEW_IMAGE_WIDTH / image_width
                x, y, width, height = [
                    float(graphics.get(attr)) * scale
                    for attr in ("x", "y", "width", "height")
                ]

                # Shift x, y coordinates from center to top-left
                gene_groups.append(
                    {
                        "view": {
                            "x": x - width / 2,
                            "y": y - height / 2,
                            "width": width,
                            "height": height,
                        },
                        "core": any(genes.values()),
                        "genes": [{"name": k, "core": v} for k, v in genes.items()],
                    }
                )

            data["pathviews"].append(
                {
                    "id": pathway_id,
                    "path": f"./{relative_remote_path}",
                }
            )

            data["pathwayIdToGeneGroups"][pathway_id] = gene_groups

        pathway_ids = set(x["pathwayId"] for x in data["pathwayData"])
        peids = {k: v for k, v in _pathway_entrez_ids().items() if k in pathway_ids}
        data["pathwayIdToGenes"] = _make_injectable(peids)

        data["pathwayIdToGeneGroups"] = _make_injectable(data["pathwayIdToGeneGroups"])
        data["pathwayData"] = _make_injectable(data["pathwayData"])

    print("Writing report HTML")
    with Path("./res/Report.html").open("w") as fw:
        with Path("./template.html").open("r") as f:
            template = Template(f.read())
            html = template.render(data)
            fw.write(html)

    print(f"Uploading results to {output}")

    local_results_directory = Path("./res").resolve()
    main_results = LatchDir(str(local_results_directory), output)

    return main_results, pathviews


metadata = LatchMetadata(
    display_name="Gene Ontology and Pathway Analysis",
    author=LatchAuthor(),
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
    number_of_pathways: int = 20,
    output_location: LatchDir = LatchDir("latch:///Pathway Analysis/"),
) -> tuple[LatchDir, list[LatchFile]]:
    """Use differential expression contrast data to calculate the gene ontology and pathways for the most significant genes

    Gene ontology
    https://www.pnas.org/doi/10.1073/pnas.0506580102
    """
    return go_pathway(
        contrast_csv=contrast_csv,
        report_name=report_name,
        number_of_pathways=number_of_pathways,
        output_location=output_location,
    )
