import pandas as pd
import numpy as np
from utils.window import build_TSS_window
import argparse

def clean_data():
    df = pd.read_csv("/ceph/ibmi/it/thesis_data/neumannarthur/data/ecoli_data/software/test_results/MasterTable.tsv", sep='\t')

    filtered_df = df[(df["detected"] == 1) & (df["enriched"] == 1)].copy()
    data = filtered_df.drop_duplicates(subset=["SuperStrand", "SuperPos"])

    def _remove_greater_100(column):
        return pd.to_numeric(column.replace({">100": "101", "> 100": "101"}), errors="coerce")

    data_clean = data.copy()
    data_clean["enrichmentFactor"] = _remove_greater_100(data["enrichmentFactor"])
    data_clean["stepFactor"] = _remove_greater_100(data["stepFactor"])
    data_clean["stepHeight"] = pd.to_numeric(data_clean["stepHeight"])

    data_by_pos = (
        data_clean.sort_values(
            by=["contigPos", "enrichmentFactor", "stepFactor", "stepHeight"],
            ascending=[True, False, False, False]
        )
        .drop_duplicates(subset="contigPos", keep="first")
    )

    data_by_pos_condition = (
        data_clean.sort_values(
            by=["contigPos", "Genome", "enrichmentFactor", "stepFactor", "stepHeight"],
            ascending=[True, True, False, False, False]
        )
        .drop_duplicates(subset=["contigPos", "Genome"], keep="first")
    )

    data_by_pos.to_csv("/ceph/ibmi/it/thesis_data/neumannarthur/results/data_by_position.csv", sep="\t", index=False)
    data_by_pos_condition.to_csv("/ceph/ibmi/it/thesis_data/neumannarthur/results/data_by_pos_condition.csv", sep="\t", index=False)

    return data_by_pos, data_by_pos_condition



def main():
    parser = argparse.ArgumentParser(description="TSS Predictor")
    subparsers = parser.add_subparsers(dest="command", required=True)

    data_parser = subparsers.add_parser("data", help="Data-related commands. Contains subcommands for cleainng, viewing and windowing given datasets.")
    data_subparser = data_parser.add_subparsers(dest="data", required=True)

    data_subparser.add_parser("clean", help="Cleans the raw dataset and outputs as processed file")

    view_parser = data_subparser.add_parser("view", help="Cleans the raw dataset and outputs as processed file.")
    view_parser.add_argument(
        "--dataset",
        nargs="?",
        const="Master",
        choices=["Master", "Clean"],
        help="Which dataset to view (default Clean)",
    )
    view_parser.add_argument(
        "--rows",
        type=int,
        default=5000,
        help="[Optional] How many rows to see in the CLI preview (default 5000)",
        required=False
    )
    view_parser.add_argument(
        "--save",
        type=str,
        help="[Optional] Save preview to .tsv file on specified path",
        required=False
    )

    window = data_subparser.add_parser("window", help="Perform window-based analysis on the cleaned dataset")
    window.add_argument(
        "--size",
        type=int,
        help="[Optional] Choose how far window stretches in each direction of TSS (Default 1000)"
    )
    window.add_argument(
        "--conds",
        nargs="+",
        choices=["Control", "nov", "rif", "tet"],
        default="Control",
        help="[Optional] List of conditions to filter by (default 'Control')",
    )
    window.add_argument(
        "--strand",
        choices=["+", "-", "master", "all"],
        help="[Optional] Specify which strand to window on (default does +, - and both)",
        required=False
    )
    window.add_argument(
        "--preview",
        nargs="?",
        type=str,
        help="[Optional] Choose path to save preview to",
        required=False
    )
    window.add_argument(
        "--wiggle_f",
        type=str,
        help="[Optional] Define which forward wiggle file the enrichment factors should be extracted from",
        required=False
    )
    window.add_argument(
        "--wiggle_b",
        type=str,
        help="[Optional] Define which backward wiggle file the enrichment factors should be extracted from",
        required=False
    )
    window.add_argument(
        "--output_dir",
        type=str,
        help="[Optional] Change output directory",
        required=False
    )

    args = parser.parse_args()

    if args.command == "data":
        print("We have data")

        if args.data == "view":
            row_nums = args.rows if args.rows else 5000
            save_path = args.save if args.save else None

            if args.dataset == "Master":
                master_df = pd.read_csv("/ceph/ibmi/it/thesis_data/neumannarthur/data/ecoli_data/software/test_results/MasterTable.tsv", sep="\t", nrows=row_nums)
                print(master_df.head(row_nums))
                if save_path is not None:
                    print("Saving to", save_path)
                    master_df.to_csv(f"{str(save_path)}/master_table_preview.tsv", sep="\t", index=False)
                    print(f"Succesfully saved {row_nums} row view of master to {save_path}")

            elif args.dataset == "Clean":
                clean_df = pd.read_csv("/ceph/ibmi/it/thesis_data/neumannarthur/results/data_by_pos_condition.csv", sep="\t", nrows=row_nums)
                print(clean_df.head(row_nums))
                if save_path is not None:
                    print("Saving to", save_path)
                    clean_df.to_csv(f"{str(save_path)}/clean_preview.tsv", sep="\t", index=False)
                    print(f"Succesfully saved {row_nums} row view of master to {save_path}")

        elif args.data == "window":
            dataframe = pd.read_csv("/ceph/ibmi/it/thesis_data/neumannarthur/results/data_by_pos_condition.csv", sep="\t")
            conditions = args.conds if args.conds else None
            strand = args.strand if args.strand else None
            save=True if args.preview else True # TODO: Changed logic, need to remove this param!!!
            preview_path=args.preview if args.preview else "/ceph/ibmi/it/thesis_data/neumannarthur/results/previews"
            window_size=args.size if args.size else 1000
            wiggle_forward = args.wiggle_f if args.wiggle_f else "/ceph/ibmi/it/thesis_data/neumannarthur/data/ecoli_data/wiggles/rnaseq_files/SRR21871309_VCE_Non_enriched_ctrl_div_by_8609484.0_multi_by_1000000.0_forward.wig"
            wiggle_backward = args.wiggle_b if args.wiggle_b else "/ceph/ibmi/it/thesis_data/neumannarthur/data/ecoli_data/wiggles/rnaseq_files/SRR21871309_VCE_Non_enriched_ctrl_div_by_8609484.0_multi_by_1000000.0_reverse.wig"
            output_dir = args.output_dir if args.output_dir else "/ceph/ibmi/it/thesis_data/neumannarthur/results"
            build_TSS_window(dataframe, window_size, conditions, strand, save, 50, preview_path, output_dir, wiggle_forward, wiggle_backward)

            

        elif args.data == "clean":
            print("Cleaning data...")
            data_by_pos, data_by_pos_condition = clean_data() # Tupel with cleaned_data[0] = data_by_pos, cleaned_data[1] = data_by_pos_condition
            print("Data cleaned!")


if __name__ == "__main__":
    main()