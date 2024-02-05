"""Merge FASTQs."""

from pathlib import Path

from resolwe.process import Data, DataField, ListField, Process, SchedulingClass


def get_label(data, warning):
    """Get relation partition label of data object."""
    label = None
    for relation in data.relations:
        if relation.category == "Replicate":
            label = next(
                p.label for p in relation.partitions if p.entity_id == data.entity_id
            )
        else:
            other_relation = next(
                p for p in relation.partitions if p.entity_id == data.entity_id
            )
            if other_relation:
                warning(
                    f"Sample {data.entity.name} has defined {relation.category} "
                    "relation. Samples will only be merged based on Replicate relations."
                )

    return label


def create_symlinks(paths):
    """Create symlinks with unique file names."""
    container_paths = [f"{n}_{Path(path).name}" for n, path in enumerate(paths)]
    for container_path, path in zip(container_paths, paths):
        Path(container_path).symlink_to(path)
    return container_paths


def group_paths(data_objects, warning, second_pair=False):
    """Group read paths grouped by relation labels."""
    labeled_paths = {}
    for data in data_objects:
        label = get_label(data=data, warning=warning)
        if second_pair:
            read_paths = [fastq.path for fastq in data.output.fastq2]
        else:
            read_paths = [fastq.path for fastq in data.output.fastq]
        if label in labeled_paths:
            labeled_paths[label].extend(read_paths)
        else:
            labeled_paths[label] = read_paths
    return labeled_paths.items()


class MergeFastqSingle(Process):
    """Merge single-end FASTQs into one sample.

    Samples are merged based on the defined replicate group relations
    and then uploaded as separate samples.
    """

    slug = "merge-fastq-single"
    name = "Merge FASTQ (single-end)"
    process_type = "data:mergereads:single"
    version = "2.2.2"
    category = "FASTQ processing"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"}
        },
        "relations": [{"type": "group"}],
    }
    data_name = "Merge FASTQ (single-end)"

    class Input:
        """Input fields to process MergeFastqSingle."""

        reads = ListField(
            DataField(data_type="reads:fastq:single:"),
            label="Select relations",
            description="Define and select replicate relations.",
            relation_type="group",
        )

    def run(self, inputs, outputs):
        """Run the analysis."""

        # Check if user has selected multiple read objects from the same sample
        data_by_sample = {}
        for data in inputs.reads:
            if data.entity_id in data_by_sample:
                self.warning(
                    "There are multiple read objects for "
                    f"{data.entity_name}. Using only the first one."
                )
                if int(data.id) < int(data_by_sample[data.entity_id].id):
                    data_by_sample[data.entity_id] = data
            else:
                data_by_sample[data.entity_id] = data

        reads = [*data_by_sample.values()]
        labeled_reads = group_paths(data_objects=reads, warning=self.warning)

        for label, paths in labeled_reads:
            if label is None:
                self.error(
                    "Missing replicate relations. Please make sure you have selected and defined "
                    "replicate sample relations."
                )
            symlinks = create_symlinks(paths=paths)

            self.run_process(
                slug="upload-fastq-single",
                inputs={
                    "src": symlinks,
                    "merge_lanes": True,
                },
            )

            merged_objects = Data.filter(entity__name=symlinks[0])
            # Sort by id and select the newest data object.
            merged_objects.sort(key=lambda x: x.id)
            merged_data = merged_objects[-1]
            merged_data.name = label
            merged_data.entity.name = label


class MergeFastqPaired(Process):
    """Merge paired-end FASTQs into one sample.

    Samples are merged based on the defined replicate group relations
    and then uploaded as separate samples.
    """

    slug = "merge-fastq-paired"
    name = "Merge FASTQ (paired-end)"
    process_type = "data:mergereads:paired"
    version = "2.2.2"
    category = "FASTQ processing"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"}
        },
        "relations": [{"type": "group"}],
    }
    data_name = "Merge FASTQ (paired-end)"

    class Input:
        """Input fields to process MergeFastqPaired."""

        reads = ListField(
            DataField(data_type="reads:fastq:paired:"),
            label="Select relations",
            description="Define and select Replicate relations.",
            relation_type="group",
        )

    def run(self, inputs, outputs):
        """Run the analysis."""
        # Check if user has selected multiple read objects from the same sample
        data_by_sample = {}
        for data in inputs.reads:
            if data.entity_id in data_by_sample:
                self.warning(
                    "There are multiple read objects for "
                    f"{data.entity_name}. Using only the first one."
                )
                if int(data.id) < int(data_by_sample[data.entity_id].id):
                    data_by_sample[data.entity_id] = data
            else:
                data_by_sample[data.entity_id] = data

        reads = [*data_by_sample.values()]
        labeled_reads = group_paths(data_objects=reads, warning=self.warning)
        labeled_reads_2 = group_paths(
            data_objects=reads, second_pair=True, warning=self.warning
        )
        for (label, paths), (_, paths_2) in zip(labeled_reads, labeled_reads_2):
            if label is None:
                self.error(
                    "Missing replicate relations. Please make sure you have selected and defined "
                    "replicate sample relations."
                )
            symlinks = create_symlinks(paths=paths)
            symlinks_2 = create_symlinks(paths=paths_2)
            self.run_process(
                slug="upload-fastq-paired",
                inputs={
                    "src1": symlinks,
                    "src2": symlinks_2,
                    "merge_lanes": True,
                },
            )

            merged_objects = Data.filter(entity__name=symlinks[0])
            # Sort by id and select the newest data object.
            merged_objects.sort(key=lambda x: x.id)
            merged_data = merged_objects[-1]
            merged_data.name = label
            merged_data.entity.name = label
