"""Inspect or regenerate the saved q0 result without repeating its IBPs.

In a notebook, import SavedQ0Result and pass the existing Symbolica API:
    saved = SavedQ0Result.load(api=symbolica)
    retained = saved.expression("retained")
    numerical = saved.expression("numerical", format="json")
The masters remain formal in every variant. See docs/five_loop_q0_coefficients.typ.
"""

import argparse
import json

from q0_coefficients import SavedQ0Result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("summary", "verify", "show"):
        command = commands.add_parser(name)
        command.add_argument(
            "--directory", help="saved artifacts; defaults to the bundled result"
        )
        if name == "show":
            command.add_argument(
                "--variant",
                choices=("retained", "symbolic", "numerical"),
                default="retained",
            )
            command.add_argument(
                "--format", choices=("native", "text", "json"), default="native"
            )
    command = commands.add_parser(
        "regenerate", help="postprocess a saved result; no IBP or master evaluation"
    )
    command.add_argument("--input", required=True)
    command.add_argument("--output", required=True, help="new artifact directory")
    command.add_argument(
        "--oracle",
        action="store_true",
        help="also cross-check with optional mpmath at fixed380digits",
    )
    args = parser.parse_args(argv)

    # One canonical Community host. Importing this example as a helper imports
    # no kernel; callers can supply their application's already loaded API.
    import symbolica as core

    try:
        if args.command == "regenerate":
            from q0_coefficients.generate import Q0Generator
            from symbolica.community.hepkit.vakint import Vakint

            directory = Q0Generator(core, Vakint).generate(
                args.input, args.output, validate_oracle=args.oracle
            )
            saved = SavedQ0Result.load(directory, api=core)
        else:
            saved = SavedQ0Result.load(args.directory, api=core)
        if args.command == "show":
            print(saved.expression(args.variant, format=args.format))
        elif args.command == "verify":
            print(json.dumps(saved.verify(), indent=2))
        else:
            print(json.dumps(saved.summary(), indent=2))
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(2, f"error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
