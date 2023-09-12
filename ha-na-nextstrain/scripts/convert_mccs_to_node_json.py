import argparse
from augur.utils import write_json
import json
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--mccs", help="MCCs JSON file from TreeKnit")
    parser.add_argument("--field-name", default="MCC", help="name of the annotation field to use in the node data JSON")
    parser.add_argument("--min-size", type=int, default=1, help="minimum size of MCC to include in output")
    parser.add_argument("--output", help="node data JSON file")
    args = parser.parse_args()

    data = {}
    total_mccs = 0
    total_strains = 0
    skipped_mccs = 0
    skipped_strains = 0
    with open(args.mccs, "r", encoding="utf-8") as fh:
        mccs_json = json.load(fh)

    if not "MCC_dict" in mccs_json:
        print(f"Could not find expected `MCC_dict` key in the provided JSON: {args.mcc}", file=sys.stderr)
        sys.exit(1)

    mcc_keys = list(mccs_json["MCC_dict"].keys())
    print(f"Found the following keys in the MCC dictionary: {mcc_keys}")
    selected_mcc_key = mcc_keys[0]
    print(f"Selected the first MCC key: {selected_mcc_key}")

    mccs = mccs_json["MCC_dict"][selected_mcc_key]["mccs"]

    for mcc in mccs:
        strains = mcc
        total_strains += len(strains)

        if len(strains) >= args.min_size:
            print(f"MCC {total_mccs} has {len(strains)} strains")
            for strain in strains:
                data[strain] = {
                    args.field_name: f"MCC_{total_mccs}",
                }

            total_mccs += 1
        else:
            for strain in strains:
                data[strain] = {
                    args.field_name: "unassigned",
                }

            skipped_mccs += 1
            skipped_strains += len(strains)

    write_json({"nodes": data}, args.output)

    print(f"Inspected MCCs with {total_strains} total strains.")
    print(f"Annotated {total_mccs} MCCs, totalling {total_strains - skipped_strains} strains ({((total_strains - skipped_strains) / total_strains) * 100:.0f}%).")
    print(f"Skipped {skipped_mccs} MCCs smaller than {args.min_size}, totalling {skipped_strains} strains.")
