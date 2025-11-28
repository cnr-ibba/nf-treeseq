#! /usr/bin/env python
import json
import tskit
import tszip
import datetime
import argparse


# define metadata schema for individuals and populations
POPULATION_METADATA_SCHEMA = tskit.MetadataSchema({
    "codec": "json",
    "type": "object",
    "properties": {
        "breed": {"type": "string"}
    },
    "required": ["breed"]
})


INDIVIDUAL_METADATA_SCHEMA = tskit.MetadataSchema({
    "codec": "json",
    "type": "object",
    "properties": {
        "sample_id": {"type": "string"}
    },
    "required": ["sample_id"]
})


def main():
    parser = argparse.ArgumentParser(description="Annotate tree sequence with sample metadata.")
    parser.add_argument("--input_tsz", help="Input tree sequence file")
    parser.add_argument("--sample_file", help="File with sample ID and breed information")
    parser.add_argument("--output_tsz", help="Output annotated tree sequence file")
    parser.add_argument("--software_name", help="Software name for provenance record")
    parser.add_argument("--software_version", help="Software version for provenance record")
    args = parser.parse_args()

    # Read sample metadata from the provided file
    sample_info = []
    with open(args.sample_file, "r") as f:
        sample_info = []

        for line in f:
            if line.strip():
                breed, sample_id = line.strip().split()
                sample_info.append((sample_id, breed))

    print(f"Loaded metadata for {len(sample_info)} samples.")

    # Load the input tree sequence
    input_tsz = tszip.load(args.input_tsz)

    # create a copy of the table that can be modified
    tables = input_tsz.dump_tables()

    # now I need to determine how many distinct populations (breeds) there are
    breeds = set(breed for _, breed in sample_info)

    # apply population metadata
    tables.populations.metadata_schema = POPULATION_METADATA_SCHEMA

    breed_to_id = {}

    for breed in breeds:
        metadata = {"breed": breed}
        pop_id = tables.populations.add_row(metadata=metadata)
        breed_to_id[breed] = pop_id

    # apply individual metadata and set population
    tables.individuals.metadata_schema = INDIVIDUAL_METADATA_SCHEMA

    individual_to_id = {}

    for sample_id, _ in sample_info:
        if "sample_id" not in individual_to_id:
            metadata = {"sample_id": sample_id}
            ind_id = tables.individuals.add_row(metadata=metadata)
            individual_to_id[sample_id] = ind_id

    # now set the population for each individual
    # we are talking about diploid individuals here
    for i, (sample_id, breed) in enumerate(sample_info):
        pop_id = breed_to_id[breed]
        ind_id = individual_to_id[sample_id]

        # update both nodes for the diploid individual
        for j in range(2):
            node_id = i * 2 + j
            node = tables.nodes[node_id]

            # Replace the node with an updated one
            tables.nodes[node_id] = node.replace(
                population=pop_id,
                individual=ind_id
            )

    # Create provenance record
    provenance_record = {
        "software": {
            "name": args.software_name,
            "version": args.software_version
        },
        "parameters": {
            "input_file": "data/toInfer/threads/ts300I2k.vcf.gz",
            "metadata_added": True,
            "populations_added": len(breed_to_id),
            "individuals_added": len(individual_to_id)
        },
        "timestamp": datetime.datetime.now().isoformat(),
        "description": "TreeSequence generated with threads and metadata added for populations and individuals"
    }

    # Add provenance to tables
    tables.provenances.add_row(
        timestamp=provenance_record["timestamp"],
        record=json.dumps(provenance_record)
    )

    # write a tsz file as output
    output_tsz = tables.tree_sequence()

    print(f"Num of populations: {output_tsz.num_populations}")
    print(f"Num of individuals: {output_tsz.num_individuals}")
    print(f"Num of nodes: {output_tsz.num_nodes}")

    tszip.compress(output_tsz, args.output_tsz)


if __name__ == "__main__":
    main()
