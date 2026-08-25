let
  system = "x86_64-linux";
  # Nix CI fetches this file without cloning the repository, so keep this
  # generated graph synchronized with nix/ci-workspace-graph.json.
  workspaceGraph = builtins.fromJSON ''
    {
      "normal_dependencies": {
        "clinnet": [],
        "feynkit": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "feynkit-ufo"
        ],
        "feynkit-cff": [
          "feynkit-graph"
        ],
        "feynkit-generator": [
          "feynkit-graph",
          "feynkit-model",
          "linnet"
        ],
        "feynkit-graph": [
          "feynkit-kinematics",
          "feynkit-model",
          "linnet"
        ],
        "feynkit-kinematics": [
          "linnet"
        ],
        "feynkit-model": [],
        "feynkit-py": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "feynkit-ufo"
        ],
        "feynkit-ufo": [
          "feynkit-model"
        ],
        "gammaloop-api": [
          "feynkit-ufo",
          "gammaloop-tracing-filter",
          "gammaloop-workspace-hack",
          "gammalooprs",
          "idenso",
          "linnet",
          "spenso",
          "symbolica-utils"
        ],
        "gammaloop-integration-tests": [
          "feynkit-generator",
          "feynkit-model",
          "gammaloop-api",
          "gammaloop-workspace-hack",
          "gammalooprs",
          "spenso",
          "symbolica-utils",
          "vakint"
        ],
        "gammaloop-tracing-filter": [
          "gammaloop-tracing-filter-macros",
          "gammaloop-workspace-hack",
          "spenso"
        ],
        "gammaloop-tracing-filter-macros": [],
        "gammaloop-workspace-hack": [],
        "gammalooprs": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-tracing-filter",
          "gammaloop-workspace-hack",
          "idenso",
          "linnet",
          "spenso",
          "spenso-hep-lib",
          "symbolica-utils",
          "vakint"
        ],
        "idenso": [
          "linnet",
          "spenso",
          "spenso-macros",
          "symbolica-utils"
        ],
        "kurvst": [],
        "linnest": [
          "linnet"
        ],
        "linnet": [
          "gammaloop-workspace-hack",
          "symbolica-utils"
        ],
        "linnet-py": [
          "linnet"
        ],
        "spenso": [
          "gammaloop-workspace-hack",
          "linnet",
          "spenso-macros",
          "symbolica-utils"
        ],
        "spenso-hep-lib": [
          "gammaloop-workspace-hack",
          "idenso",
          "spenso",
          "spenso-macros"
        ],
        "spenso-macros": [],
        "spynso3": [
          "idenso",
          "spenso",
          "spenso-hep-lib",
          "spenso-macros"
        ],
        "symbolica-utils": [],
        "vakint": []
      },
      "normal_dependency_features": {
        "clinnet": {},
        "feynkit": {
          "feynkit-cff": [],
          "feynkit-generator": [],
          "feynkit-graph": [],
          "feynkit-kinematics": [],
          "feynkit-model": [],
          "feynkit-ufo": []
        },
        "feynkit-cff": {
          "feynkit-graph": []
        },
        "feynkit-generator": {
          "feynkit-graph": [],
          "feynkit-model": [],
          "linnet": [
            "symbolica"
          ]
        },
        "feynkit-graph": {
          "feynkit-kinematics": [],
          "feynkit-model": [],
          "linnet": [
            "serde"
          ]
        },
        "feynkit-kinematics": {
          "linnet": []
        },
        "feynkit-model": {},
        "feynkit-py": {
          "feynkit-cff": [],
          "feynkit-generator": [],
          "feynkit-graph": [],
          "feynkit-kinematics": [],
          "feynkit-model": [],
          "feynkit-ufo": []
        },
        "feynkit-ufo": {
          "feynkit-model": []
        },
        "gammaloop-api": {
          "feynkit-ufo": [],
          "gammaloop-tracing-filter": [
            "clap"
          ],
          "gammaloop-workspace-hack": [],
          "gammalooprs": [],
          "idenso": [],
          "linnet": [
            "bincode",
            "serde",
            "symbolica"
          ],
          "spenso": [
            "shadowing"
          ],
          "symbolica-utils": []
        },
        "gammaloop-integration-tests": {
          "feynkit-generator": [],
          "feynkit-model": [],
          "gammaloop-api": [],
          "gammaloop-workspace-hack": [],
          "gammalooprs": [],
          "spenso": [
            "shadowing"
          ],
          "symbolica-utils": [],
          "vakint": []
        },
        "gammaloop-tracing-filter": {
          "gammaloop-tracing-filter-macros": [],
          "gammaloop-workspace-hack": [],
          "spenso": [
            "shadowing"
          ]
        },
        "gammaloop-tracing-filter-macros": {},
        "gammaloop-workspace-hack": {},
        "gammalooprs": {
          "feynkit-cff": [],
          "feynkit-generator": [],
          "feynkit-graph": [],
          "feynkit-kinematics": [],
          "feynkit-model": [],
          "gammaloop-tracing-filter": [
            "symbolica"
          ],
          "gammaloop-workspace-hack": [],
          "idenso": [
            "bincode"
          ],
          "linnet": [
            "bincode",
            "serde",
            "symbolica"
          ],
          "spenso": [
            "shadowing"
          ],
          "spenso-hep-lib": [],
          "symbolica-utils": [],
          "vakint": []
        },
        "idenso": {
          "linnet": [],
          "spenso": [
            "shadowing"
          ],
          "spenso-macros": [],
          "symbolica-utils": []
        },
        "kurvst": {},
        "linnest": {
          "linnet": [
            "drawing",
            "rkyv",
            "serde"
          ]
        },
        "linnet": {
          "gammaloop-workspace-hack": [],
          "symbolica-utils": []
        },
        "linnet-py": {
          "linnet": []
        },
        "spenso": {
          "gammaloop-workspace-hack": [],
          "linnet": [
            "bincode",
            "serde"
          ],
          "spenso-macros": [],
          "symbolica-utils": []
        },
        "spenso-hep-lib": {
          "gammaloop-workspace-hack": [],
          "idenso": [
            "reference-cases"
          ],
          "spenso": [
            "shadowing"
          ],
          "spenso-macros": []
        },
        "spenso-macros": {},
        "spynso3": {
          "idenso": [],
          "spenso": [
            "python",
            "shadowing"
          ],
          "spenso-hep-lib": [],
          "spenso-macros": []
        },
        "symbolica-utils": {},
        "vakint": {}
      },
      "package_dirs": {
        "clinnet": "crates/clinnet",
        "feynkit": "crates/feynkit",
        "feynkit-cff": "crates/feynkit-cff",
        "feynkit-generator": "crates/feynkit-generator",
        "feynkit-graph": "crates/feynkit-graph",
        "feynkit-kinematics": "crates/feynkit-kinematics",
        "feynkit-model": "crates/feynkit-model",
        "feynkit-py": "crates/feynkit-py",
        "feynkit-ufo": "crates/feynkit-ufo",
        "gammaloop-api": "crates/gammaloop-api",
        "gammaloop-integration-tests": "tests",
        "gammaloop-tracing-filter": "crates/gammaloop-tracing-filter",
        "gammaloop-tracing-filter-macros": "crates/gammaloop-tracing-filter-macros",
        "gammaloop-workspace-hack": "crates/gammaloop-workspace-hack",
        "gammalooprs": "crates/gammalooprs",
        "idenso": "crates/idenso",
        "kurvst": "crates/kurvst",
        "linnest": "crates/linnest",
        "linnet": "crates/linnet",
        "linnet-py": "crates/linnet-py",
        "spenso": "crates/spenso",
        "spenso-hep-lib": "crates/spenso-hep-lib",
        "spenso-macros": "crates/spenso-macros",
        "spynso3": "crates/spynso3",
        "symbolica-utils": "crates/symbolica-utils",
        "vakint": "crates/vakint"
      },
      "packages": [
        "clinnet",
        "feynkit",
        "feynkit-cff",
        "feynkit-generator",
        "feynkit-graph",
        "feynkit-kinematics",
        "feynkit-model",
        "feynkit-py",
        "feynkit-ufo",
        "gammaloop-api",
        "gammaloop-integration-tests",
        "gammaloop-tracing-filter",
        "gammaloop-tracing-filter-macros",
        "gammaloop-workspace-hack",
        "gammalooprs",
        "idenso",
        "kurvst",
        "linnest",
        "linnet",
        "linnet-py",
        "spenso",
        "spenso-hep-lib",
        "spenso-macros",
        "spynso3",
        "symbolica-utils",
        "vakint"
      ],
      "resolved_normal_dependencies": {
        "clinnet": [],
        "feynkit": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-workspace-hack",
          "linnet",
          "symbolica-utils"
        ],
        "feynkit-cff": [
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-workspace-hack",
          "linnet"
        ],
        "feynkit-generator": [
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-workspace-hack",
          "linnet",
          "symbolica-utils"
        ],
        "feynkit-graph": [
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-workspace-hack",
          "linnet"
        ],
        "feynkit-kinematics": [
          "gammaloop-workspace-hack",
          "linnet"
        ],
        "feynkit-model": [],
        "feynkit-py": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "feynkit-ufo",
          "gammaloop-workspace-hack",
          "linnet",
          "symbolica-utils"
        ],
        "feynkit-ufo": [
          "feynkit-model"
        ],
        "gammaloop-api": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-tracing-filter",
          "gammaloop-tracing-filter-macros",
          "gammaloop-workspace-hack",
          "gammalooprs",
          "idenso",
          "linnet",
          "spenso",
          "spenso-hep-lib",
          "spenso-macros",
          "symbolica-utils",
          "vakint"
        ],
        "gammaloop-integration-tests": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-api",
          "gammaloop-tracing-filter",
          "gammaloop-tracing-filter-macros",
          "gammaloop-workspace-hack",
          "gammalooprs",
          "idenso",
          "linnet",
          "spenso",
          "spenso-hep-lib",
          "spenso-macros",
          "symbolica-utils",
          "vakint"
        ],
        "gammaloop-tracing-filter": [
          "gammaloop-tracing-filter-macros",
          "gammaloop-workspace-hack",
          "linnet",
          "spenso",
          "spenso-macros",
          "symbolica-utils"
        ],
        "gammaloop-tracing-filter-macros": [],
        "gammaloop-workspace-hack": [],
        "gammalooprs": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-tracing-filter",
          "gammaloop-tracing-filter-macros",
          "gammaloop-workspace-hack",
          "idenso",
          "linnet",
          "spenso",
          "spenso-hep-lib",
          "spenso-macros",
          "symbolica-utils",
          "vakint"
        ],
        "idenso": [
          "gammaloop-workspace-hack",
          "linnet",
          "spenso",
          "spenso-macros",
          "symbolica-utils"
        ],
        "kurvst": [],
        "linnest": [
          "gammaloop-workspace-hack",
          "linnet"
        ],
        "linnet": [
          "gammaloop-workspace-hack",
          "symbolica-utils"
        ],
        "linnet-py": [
          "gammaloop-workspace-hack",
          "linnet"
        ],
        "spenso": [
          "gammaloop-workspace-hack",
          "linnet",
          "spenso-macros",
          "symbolica-utils"
        ],
        "spenso-hep-lib": [
          "gammaloop-workspace-hack",
          "idenso",
          "linnet",
          "spenso",
          "spenso-macros",
          "symbolica-utils"
        ],
        "spenso-macros": [],
        "spynso3": [
          "gammaloop-workspace-hack",
          "idenso",
          "linnet",
          "spenso",
          "spenso-hep-lib",
          "spenso-macros",
          "symbolica-utils"
        ],
        "symbolica-utils": [],
        "vakint": []
      },
      "resolved_test_dependencies": {
        "clinnet": [],
        "feynkit": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-workspace-hack",
          "linnet",
          "symbolica-utils"
        ],
        "feynkit-cff": [
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-workspace-hack",
          "linnet"
        ],
        "feynkit-generator": [
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-workspace-hack",
          "linnet",
          "symbolica-utils"
        ],
        "feynkit-graph": [
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-workspace-hack",
          "linnet"
        ],
        "feynkit-kinematics": [
          "gammaloop-workspace-hack",
          "linnet"
        ],
        "feynkit-model": [],
        "feynkit-py": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "feynkit-ufo",
          "gammaloop-workspace-hack",
          "idenso",
          "linnet",
          "spenso",
          "spenso-hep-lib",
          "spenso-macros",
          "spynso3",
          "symbolica-utils"
        ],
        "feynkit-ufo": [
          "feynkit-model"
        ],
        "gammaloop-api": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-tracing-filter",
          "gammaloop-tracing-filter-macros",
          "gammaloop-workspace-hack",
          "gammalooprs",
          "idenso",
          "linnet",
          "spenso",
          "spenso-hep-lib",
          "spenso-macros",
          "symbolica-utils",
          "vakint"
        ],
        "gammaloop-integration-tests": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-api",
          "gammaloop-tracing-filter",
          "gammaloop-tracing-filter-macros",
          "gammaloop-workspace-hack",
          "gammalooprs",
          "idenso",
          "linnet",
          "spenso",
          "spenso-hep-lib",
          "spenso-macros",
          "symbolica-utils",
          "vakint"
        ],
        "gammaloop-tracing-filter": [
          "gammaloop-tracing-filter-macros",
          "gammaloop-workspace-hack",
          "linnet",
          "spenso",
          "spenso-macros",
          "symbolica-utils"
        ],
        "gammaloop-tracing-filter-macros": [],
        "gammaloop-workspace-hack": [],
        "gammalooprs": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-tracing-filter",
          "gammaloop-tracing-filter-macros",
          "gammaloop-workspace-hack",
          "idenso",
          "linnet",
          "spenso",
          "spenso-hep-lib",
          "spenso-macros",
          "symbolica-utils",
          "vakint"
        ],
        "idenso": [
          "gammaloop-workspace-hack",
          "linnet",
          "spenso",
          "spenso-macros",
          "symbolica-utils"
        ],
        "kurvst": [],
        "linnest": [
          "gammaloop-workspace-hack",
          "linnet"
        ],
        "linnet": [
          "gammaloop-workspace-hack",
          "linnest",
          "symbolica-utils"
        ],
        "linnet-py": [
          "gammaloop-workspace-hack",
          "linnet"
        ],
        "spenso": [
          "gammaloop-workspace-hack",
          "linnet",
          "spenso-macros",
          "symbolica-utils"
        ],
        "spenso-hep-lib": [
          "gammaloop-workspace-hack",
          "idenso",
          "linnet",
          "spenso",
          "spenso-macros",
          "symbolica-utils"
        ],
        "spenso-macros": [
          "gammaloop-workspace-hack",
          "linnet",
          "spenso",
          "symbolica-utils"
        ],
        "spynso3": [
          "gammaloop-workspace-hack",
          "idenso",
          "linnet",
          "spenso",
          "spenso-hep-lib",
          "spenso-macros",
          "symbolica-utils"
        ],
        "symbolica-utils": [],
        "vakint": []
      },
      "symbolica_normal_packages": [
        "feynkit",
        "feynkit-cff",
        "feynkit-generator",
        "feynkit-graph",
        "feynkit-kinematics",
        "feynkit-py",
        "gammaloop-api",
        "gammaloop-integration-tests",
        "gammaloop-tracing-filter",
        "gammaloop-workspace-hack",
        "gammalooprs",
        "idenso",
        "linnest",
        "linnet",
        "linnet-py",
        "spenso",
        "spenso-hep-lib",
        "spenso-macros",
        "spynso3",
        "symbolica-utils",
        "vakint"
      ],
      "symbolica_test_packages": [
        "feynkit",
        "feynkit-cff",
        "feynkit-generator",
        "feynkit-graph",
        "feynkit-kinematics",
        "feynkit-py",
        "gammaloop-api",
        "gammaloop-integration-tests",
        "gammaloop-tracing-filter",
        "gammaloop-workspace-hack",
        "gammalooprs",
        "idenso",
        "linnest",
        "linnet",
        "linnet-py",
        "spenso",
        "spenso-hep-lib",
        "spenso-macros",
        "spynso3",
        "symbolica-utils",
        "vakint"
      ],
      "test_dependencies": {
        "clinnet": [],
        "feynkit": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "feynkit-ufo"
        ],
        "feynkit-cff": [
          "feynkit-graph"
        ],
        "feynkit-generator": [
          "feynkit-graph",
          "feynkit-model",
          "linnet"
        ],
        "feynkit-graph": [
          "feynkit-kinematics",
          "feynkit-model",
          "linnet"
        ],
        "feynkit-kinematics": [
          "linnet"
        ],
        "feynkit-model": [],
        "feynkit-py": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "feynkit-ufo",
          "spynso3"
        ],
        "feynkit-ufo": [
          "feynkit-model"
        ],
        "gammaloop-api": [
          "feynkit-ufo",
          "gammaloop-tracing-filter",
          "gammaloop-workspace-hack",
          "gammalooprs",
          "idenso",
          "linnet",
          "spenso",
          "symbolica-utils"
        ],
        "gammaloop-integration-tests": [
          "feynkit-generator",
          "feynkit-model",
          "gammaloop-api",
          "gammaloop-workspace-hack",
          "gammalooprs",
          "spenso",
          "symbolica-utils",
          "vakint"
        ],
        "gammaloop-tracing-filter": [
          "gammaloop-tracing-filter-macros",
          "gammaloop-workspace-hack",
          "spenso"
        ],
        "gammaloop-tracing-filter-macros": [],
        "gammaloop-workspace-hack": [],
        "gammalooprs": [
          "feynkit-cff",
          "feynkit-generator",
          "feynkit-graph",
          "feynkit-kinematics",
          "feynkit-model",
          "gammaloop-tracing-filter",
          "gammaloop-workspace-hack",
          "idenso",
          "linnet",
          "spenso",
          "spenso-hep-lib",
          "symbolica-utils",
          "vakint"
        ],
        "idenso": [
          "linnet",
          "spenso",
          "spenso-macros",
          "symbolica-utils"
        ],
        "kurvst": [],
        "linnest": [
          "linnet"
        ],
        "linnet": [
          "gammaloop-workspace-hack",
          "linnest",
          "symbolica-utils"
        ],
        "linnet-py": [
          "linnet"
        ],
        "spenso": [
          "gammaloop-workspace-hack",
          "linnet",
          "spenso-macros",
          "symbolica-utils"
        ],
        "spenso-hep-lib": [
          "gammaloop-workspace-hack",
          "idenso",
          "spenso",
          "spenso-macros"
        ],
        "spenso-macros": [
          "linnet",
          "spenso"
        ],
        "spynso3": [
          "idenso",
          "spenso",
          "spenso-hep-lib",
          "spenso-macros"
        ],
        "symbolica-utils": [],
        "vakint": []
      },
      "test_dependency_features": {
        "clinnet": {},
        "feynkit": {
          "feynkit-cff": [],
          "feynkit-generator": [],
          "feynkit-graph": [],
          "feynkit-kinematics": [],
          "feynkit-model": [],
          "feynkit-ufo": []
        },
        "feynkit-cff": {
          "feynkit-graph": []
        },
        "feynkit-generator": {
          "feynkit-graph": [],
          "feynkit-model": [],
          "linnet": [
            "symbolica"
          ]
        },
        "feynkit-graph": {
          "feynkit-kinematics": [],
          "feynkit-model": [],
          "linnet": [
            "serde"
          ]
        },
        "feynkit-kinematics": {
          "linnet": []
        },
        "feynkit-model": {},
        "feynkit-py": {
          "feynkit-cff": [],
          "feynkit-generator": [],
          "feynkit-graph": [],
          "feynkit-kinematics": [],
          "feynkit-model": [],
          "feynkit-ufo": [],
          "spynso3": []
        },
        "feynkit-ufo": {
          "feynkit-model": []
        },
        "gammaloop-api": {
          "feynkit-ufo": [],
          "gammaloop-tracing-filter": [
            "clap"
          ],
          "gammaloop-workspace-hack": [],
          "gammalooprs": [],
          "idenso": [],
          "linnet": [
            "bincode",
            "serde",
            "symbolica"
          ],
          "spenso": [
            "shadowing"
          ],
          "symbolica-utils": []
        },
        "gammaloop-integration-tests": {
          "feynkit-generator": [],
          "feynkit-model": [],
          "gammaloop-api": [],
          "gammaloop-workspace-hack": [],
          "gammalooprs": [],
          "spenso": [
            "shadowing"
          ],
          "symbolica-utils": [],
          "vakint": []
        },
        "gammaloop-tracing-filter": {
          "gammaloop-tracing-filter-macros": [],
          "gammaloop-workspace-hack": [],
          "spenso": [
            "shadowing"
          ]
        },
        "gammaloop-tracing-filter-macros": {},
        "gammaloop-workspace-hack": {},
        "gammalooprs": {
          "feynkit-cff": [],
          "feynkit-generator": [],
          "feynkit-graph": [],
          "feynkit-kinematics": [],
          "feynkit-model": [],
          "gammaloop-tracing-filter": [
            "symbolica"
          ],
          "gammaloop-workspace-hack": [],
          "idenso": [
            "bincode"
          ],
          "linnet": [
            "bincode",
            "serde",
            "symbolica"
          ],
          "spenso": [
            "shadowing"
          ],
          "spenso-hep-lib": [],
          "symbolica-utils": [],
          "vakint": []
        },
        "idenso": {
          "linnet": [],
          "spenso": [
            "shadowing"
          ],
          "spenso-macros": [],
          "symbolica-utils": []
        },
        "kurvst": {},
        "linnest": {
          "linnet": [
            "drawing",
            "rkyv",
            "serde"
          ]
        },
        "linnet": {
          "gammaloop-workspace-hack": [],
          "linnest": [],
          "symbolica-utils": []
        },
        "linnet-py": {
          "linnet": []
        },
        "spenso": {
          "gammaloop-workspace-hack": [],
          "linnet": [
            "bincode",
            "serde"
          ],
          "spenso-macros": [],
          "symbolica-utils": []
        },
        "spenso-hep-lib": {
          "gammaloop-workspace-hack": [],
          "idenso": [
            "reference-cases"
          ],
          "spenso": [
            "shadowing"
          ],
          "spenso-macros": []
        },
        "spenso-macros": {
          "linnet": [],
          "spenso": []
        },
        "spynso3": {
          "idenso": [],
          "spenso": [
            "python",
            "shadowing"
          ],
          "spenso-hep-lib": [],
          "spenso-macros": []
        },
        "symbolica-utils": {},
        "vakint": {}
      },
      "workspace_root": "."
    }
  '';
  unique = values:
    builtins.attrNames (builtins.listToAttrs (map (value: {
        name = value;
        value = true;
      })
      values));
  workspacePackages = workspaceGraph.packages;
  cratePackageDepsAttr = package: "packages.${system}.crate-deps-${package}";
  cratePackageAttr = package: "packages.${system}.crate-${package}";
  crateTestDependencyAttr = representative: "packages.${system}.crate-test-dependencies-${representative}";
  crateTestBinaryAttr = package: "packages.${system}.crate-test-binaries-${package}";
  nextestContextualTestDependencyAttr = target: package:
    "packages.${system}.crate-test-dependencies-${target}-${package}";
  nextestContextualTestBinaryAttr = target: package:
    "packages.${system}.crate-test-binaries-${target}-${package}";
  workspaceHackPackage = "gammaloop-workspace-hack";
  # clinnet is binary-only: the flake exposes crate-clinnet, but no
  # crate-deps-clinnet artifact.
  workspacePackagesWithDependencyArtifacts = builtins.filter (package: package != "clinnet") workspacePackages;
  nonWorkspaceHackPackages = builtins.filter (package: package != workspaceHackPackage) workspacePackages;
  workspaceHackCacheAttr = cratePackageDepsAttr workspaceHackPackage;
  gammaloopApiPackageArtifactsAttr = "packages.${system}.gammaloopApiPackageArtifacts";
  workspacePackageGraphAttr = package: cratePackageAttr package;
  mergeDependencySets = sets: let
    attrs = unique (builtins.concatLists (map builtins.attrNames sets));
  in
    builtins.listToAttrs (map (attr: {
        name = attr;
        value = unique (builtins.concatLists (map (set: set.${attr} or []) sets));
      })
      attrs);
  workspaceDependencyNamesFor = package:
    workspaceGraph.resolved_normal_dependencies.${package} or [];
  workspaceTestDependencyNamesFor = package:
    workspaceGraph.test_dependencies.${package} or (workspaceGraph.normal_dependencies.${package} or []);
  workspaceDependencyClosureFor = dependencyNamesFor: package:
    unique (map (entry: entry.key) (builtins.genericClosure {
      startSet = [{key = package;}];
      operator = entry: map (dependency: {key = dependency;}) (dependencyNamesFor entry.key);
    }));
  workspaceTestDependencyClosureFor =
    workspaceDependencyClosureFor workspaceTestDependencyNamesFor;
  workspaceTestComponentMembersFor = package:
    unique (builtins.filter (
        other:
          builtins.elem other (workspaceTestDependencyClosureFor package)
          && builtins.elem package (workspaceTestDependencyClosureFor other)
      )
      workspacePackages);
  workspaceTestComponentRepresentativeFor = package:
    builtins.head (workspaceTestComponentMembersFor package);
  workspaceTestComponentRepresentatives =
    unique (map workspaceTestComponentRepresentativeFor workspacePackages);
  workspaceTestDependencyComponentRepresentatives =
    builtins.filter (representative: representative != workspaceHackPackage) workspaceTestComponentRepresentatives;
  workspaceTestComponentMembers =
    builtins.listToAttrs (map (representative: {
        name = representative;
        value = workspaceTestComponentMembersFor representative;
      })
      workspaceTestComponentRepresentatives);
  workspaceTestComponentDependencyRepresentativesFor = representative:
    unique (builtins.filter (dependencyRepresentative: dependencyRepresentative != representative) (
      map workspaceTestComponentRepresentativeFor (
        builtins.concatLists (map workspaceTestDependencyNamesFor workspaceTestComponentMembers.${representative})
      )
    ));
  workspaceCratePackageDependencies = builtins.listToAttrs (
    builtins.filter (entry: entry.value != []) (map (package: {
        name = workspacePackageGraphAttr package;
        value = map workspacePackageGraphAttr (workspaceDependencyNamesFor package);
      })
      workspacePackages)
  );
  workspaceCratePackageDependencyEdges = builtins.concatLists (map (dependent:
      map (dependency: {
        inherit dependency dependent;
      })
      (workspaceCratePackageDependencies.${dependent} or []))
    (builtins.attrNames workspaceCratePackageDependencies));
  workspaceCratePackageCacheDependencies = builtins.listToAttrs (map (package: {
      name = workspacePackageGraphAttr package;
      value = [(cratePackageDepsAttr package)];
    })
    workspacePackagesWithDependencyArtifacts);
  workspaceCratePackageCacheArtifactDependencies = builtins.listToAttrs (
    builtins.filter (entry: entry.value != []) (map (package: {
        name = cratePackageDepsAttr package;
        value =
          (
            if package == workspaceHackPackage
            then []
            else ["packages.${system}.cargoArtifacts"]
          )
          ++ map cratePackageDepsAttr (workspaceDependencyNamesFor package)
          ++ (
            if package != workspaceHackPackage && (workspaceDependencyNamesFor package) == [] && builtins.elem package workspaceGraph.symbolica_normal_packages
            then [workspaceHackCacheAttr]
            else []
          );
      })
      workspacePackagesWithDependencyArtifacts)
  );
  workspaceCratePackageCacheArtifactDependencyEdges = builtins.concatLists (map (dependent:
      map (dependency: {
        inherit dependency dependent;
      })
      (workspaceCratePackageCacheArtifactDependencies.${dependent} or []))
    (builtins.attrNames workspaceCratePackageCacheArtifactDependencies));
  workspaceTestDependencyArtifactDependencies = builtins.listToAttrs (map (representative: {
      name = crateTestDependencyAttr representative;
      value =
        [
          "packages.${system}.cargoArtifacts"
          workspaceHackCacheAttr
        ]
        ++ map crateTestDependencyAttr (
          builtins.filter (
            dependencyRepresentative: dependencyRepresentative != workspaceHackPackage
          )
          (workspaceTestComponentDependencyRepresentativesFor representative)
        );
    })
    workspaceTestDependencyComponentRepresentatives);
  workspaceTestBinaryArtifactDependencies = builtins.listToAttrs (map (package: {
      name = crateTestBinaryAttr package;
      value = [(crateTestDependencyAttr (workspaceTestComponentRepresentativeFor package))];
    })
    nonWorkspaceHackPackages);
  nextestPackageGroups = {
    core = [
      "gammaloop-api"
      "gammaloop-tracing-filter"
      "gammaloop-tracing-filter-macros"
      "gammalooprs"
    ];
    integration = ["gammaloop-integration-tests"];
    feynkit = [
      "feynkit"
      "feynkit-cff"
      "feynkit-generator"
      "feynkit-graph"
      "feynkit-kinematics"
      "feynkit-model"
      "feynkit-py"
      "feynkit-ufo"
    ];
    "python-api" = ["gammaloop-integration-tests"];
    clinnet = ["clinnet"];
    linnet = [
      "kurvst"
      "linnet"
      "linnet-py"
      "linnest"
    ];
    spenso = [
      "idenso"
      "spenso"
      "spenso-hep-lib"
      "spenso-macros"
      "symbolica-utils"
    ];
    vakint = ["vakint"];
  };
  nextestArchiveAttr = target: "checks.${system}.gammaloop-nextest-binaries-${target}";
  nextestPackageArtifactAttrFor = target: package:
    if target == "python-api"
    then nextestContextualTestBinaryAttr target package
    else crateTestBinaryAttr package;
  nextestArchiveDependenciesFor = target:
    ["packages.${system}.cargoArtifacts"]
    ++ unique (map (nextestPackageArtifactAttrFor target) nextestPackageGroups.${target})
    ++ (
      if target == "python-api"
      then ["packages.${system}.gammaloop-python-module"]
      else []
    );
  # The Hakari workspace-hack deps artifact is the root for the
  # Symbolica-containing cache DAG. Higher-level crate cache jobs reach it
  # through their Guppy-resolved workspace cache dependencies.
  nextestBinaryChecks = map nextestArchiveAttr (builtins.attrNames nextestPackageGroups);
  dependencies = mergeDependencySets [
    workspaceCratePackageCacheArtifactDependencies
    workspaceCratePackageCacheDependencies
    workspaceCratePackageDependencies
    workspaceTestDependencyArtifactDependencies
    workspaceTestBinaryArtifactDependencies
    {
      "packages.${system}.gammaloop" = [
        gammaloopApiPackageArtifactsAttr
        "checks.${system}.gammaloop-fmt"
      ];
      "checks.${system}.gammaloop" = ["packages.${system}.gammaloop"];
      "packages.${system}.default" = ["packages.${system}.gammaloop"];
      "packages.${system}.gammaloop-python-module" =
        workspaceCratePackageDependencies.${cratePackageAttr "gammaloop-api"} or [];
      ${gammaloopApiPackageArtifactsAttr} = [(cratePackageAttr "gammaloop-api")];
      "packages.${system}.cargoArtifacts" = [workspaceHackCacheAttr];
      ${nextestContextualTestDependencyAttr "python-api" "gammaloop-integration-tests"} =
        workspaceTestDependencyArtifactDependencies.${crateTestDependencyAttr (workspaceTestComponentRepresentativeFor "gammaloop-integration-tests")}
        ++ ["packages.${system}.gammaloop-python-module"];
      ${nextestContextualTestBinaryAttr "python-api" "gammaloop-integration-tests"} = [
        (nextestContextualTestDependencyAttr "python-api" "gammaloop-integration-tests")
        "packages.${system}.gammaloop-python-module"
      ];
      "checks.${system}.gammaloop-check" = ["packages.${system}.cargoArtifacts"];
      "checks.${system}.gammaloop-clippy" = ["packages.${system}.cargoArtifacts"];
      "checks.${system}.gammaloop-doc" = ["packages.${system}.cargoArtifacts"];
      "checks.${system}.gammaloop-doctest" = ["packages.${system}.cargoArtifacts"];
      "packages.${system}.workspaceBuildArtifacts" = ["packages.${system}.cargoArtifacts"];
      "checks.${system}.gammaloop-nextest-binaries-core" = nextestArchiveDependenciesFor "core";
      "checks.${system}.gammaloop-nextest-binaries-clinnet" = nextestArchiveDependenciesFor "clinnet";
      "checks.${system}.gammaloop-nextest-binaries-feynkit" = nextestArchiveDependenciesFor "feynkit";
      "checks.${system}.gammaloop-nextest-binaries-integration" = nextestArchiveDependenciesFor "integration";
      "checks.${system}.gammaloop-nextest-binaries-python-api" = nextestArchiveDependenciesFor "python-api";
      "checks.${system}.gammaloop-nextest-binaries-linnet" = nextestArchiveDependenciesFor "linnet";
      "checks.${system}.gammaloop-nextest-binaries-spenso" = nextestArchiveDependenciesFor "spenso";
      "checks.${system}.gammaloop-nextest-binaries-vakint" = nextestArchiveDependenciesFor "vakint";
      "checks.${system}.gammaloop-nextest-binaries" = nextestBinaryChecks;
      "packages.${system}.linnest-wasm" = ["packages.${system}.linnestWasmCargoArtifacts"];
      "checks.${system}.linnest-wasm" = ["packages.${system}.linnest-wasm"];
      "packages.${system}.gammaloop-llvm-coverage" = ["packages.${system}.gammaloop"];
      "packages.${system}.nix-ci-check-gammaloop-doctest" = ["packages.${system}.cargoArtifacts"];
      "packages.${system}.nix-ci-check-gammaloop-nextest" =
        nextestBinaryChecks
        ++ ["packages.${system}.gammaloop-python-module"];
      "packages.${system}.nix-ci-check-gammaloop-nextest-clinnet" = ["checks.${system}.gammaloop-nextest-binaries-clinnet"];
      "packages.${system}.nix-ci-check-gammaloop-nextest-core" = ["checks.${system}.gammaloop-nextest-binaries-core"];
      "packages.${system}.nix-ci-check-gammaloop-nextest-feynkit" = ["checks.${system}.gammaloop-nextest-binaries-feynkit"];
      "packages.${system}.nix-ci-check-gammaloop-nextest-integration" = [
        "checks.${system}.gammaloop-nextest-binaries-integration"
      ];
      "packages.${system}.nix-ci-check-gammaloop-nextest-python-api" = [
        "checks.${system}.gammaloop-nextest-binaries-python-api"
        "packages.${system}.gammaloop-python-module"
      ];
      "packages.${system}.nix-ci-check-gammaloop-nextest-linnet" = ["checks.${system}.gammaloop-nextest-binaries-linnet"];
      "packages.${system}.nix-ci-check-gammaloop-nextest-spenso" = ["checks.${system}.gammaloop-nextest-binaries-spenso"];
      "packages.${system}.nix-ci-check-gammaloop-nextest-vakint" = ["checks.${system}.gammaloop-nextest-binaries-vakint"];
    }
  ];
  missingWorkspaceCratePackageEdges =
    builtins.filter (
      edge: !(builtins.elem edge.dependency (dependencies.${edge.dependent} or []))
    )
    workspaceCratePackageDependencyEdges;
  missingWorkspaceCratePackageCacheArtifactEdges =
    builtins.filter (
      edge: !(builtins.elem edge.dependency (dependencies.${edge.dependent} or []))
    )
    workspaceCratePackageCacheArtifactDependencyEdges;
  reciprocalWorkspaceCratePackageEdges =
    builtins.filter (
      edge: builtins.elem edge.dependent (workspaceCratePackageDependencies.${edge.dependency} or [])
    )
    workspaceCratePackageDependencyEdges;
  formatDependencyEdge = edge: "${edge.dependency} -> ${edge.dependent}";
  selfDependencies =
    builtins.filter (
      attr: builtins.elem attr (dependencies.${attr} or [])
    )
    (builtins.attrNames dependencies);
  validatedDependencies =
    assert missingWorkspaceCratePackageEdges == []
    || builtins.throw "manual NixCI dependency graph is missing workspace crate package dependency edges: ${builtins.concatStringsSep ", " (map formatDependencyEdge missingWorkspaceCratePackageEdges)}";
    assert missingWorkspaceCratePackageCacheArtifactEdges == []
    || builtins.throw "manual NixCI dependency graph is missing workspace crate package cache dependency edges: ${builtins.concatStringsSep ", " (map formatDependencyEdge missingWorkspaceCratePackageCacheArtifactEdges)}";
    assert reciprocalWorkspaceCratePackageEdges == []
    || builtins.throw "manual NixCI dependency graph contains reciprocal workspace crate package dependency edges: ${builtins.concatStringsSep ", " (map formatDependencyEdge reciprocalWorkspaceCratePackageEdges)}";
    assert selfDependencies == []
    || builtins.throw "manual NixCI dependency graph contains self dependencies: ${builtins.concatStringsSep ", " selfDependencies}";
      dependencies;
  doNotBuild = unique (
    [
      "checks.${system}.gammaloop"
      "checks.${system}.gammaloop-doctest"
      "checks.${system}.gammaloop-nextest"
      "checks.${system}.gammaloop-nextest-binaries"
      "checks.${system}.gammaloop-nextest-clinnet"
      "checks.${system}.gammaloop-nextest-core"
      "checks.${system}.gammaloop-nextest-feynkit"
      "checks.${system}.gammaloop-nextest-integration"
      "checks.${system}.gammaloop-nextest-python-api"
      "checks.${system}.gammaloop-nextest-linnet"
      "checks.${system}.gammaloop-nextest-spenso"
      "checks.${system}.gammaloop-nextest-vakint"
      "packages.${system}.default"
      "packages.${system}.crane-ci-prebuild"
      "packages.${system}.workspaceBuildArtifacts"
      "packages.${system}.gammaloop-llvm-coverage"
      "packages.${system}.nix-ci-check-gammaloop-nextest"
    ]
    ++ [
      (crateTestDependencyAttr "spynso3")
      (crateTestBinaryAttr workspaceHackPackage)
      (crateTestBinaryAttr "spynso3")
      (workspacePackageGraphAttr workspaceHackPackage)
    ]
    ++ map cratePackageDepsAttr (
      builtins.filter (
        package: package != workspaceHackPackage && package != "gammalooprs"
      )
      workspacePackagesWithDependencyArtifacts
    )
    ++ map cratePackageAttr nonWorkspaceHackPackages
  );
  # NixCI only schedules jobs it actually builds, so a dependency edge that
  # references a doNotBuild job is rejected as pointing at a non-existent job.
  # The manual graph above is constructed over the full crate/artifact DAG
  # (which keeps the drift and cycle asserts meaningful); here hidden paths are
  # contracted to their nearest built dependency so their ordering is retained.
  doNotBuildSet = builtins.listToAttrs (map (job: {
      name = job;
      value = true;
    })
    doNotBuild);
  isBuiltJob = job: !(doNotBuildSet ? ${job});
  builtDependencyFrontierFor = deps: dependent:
    unique (map (entry: entry.key) (builtins.filter (
        entry: isBuiltJob entry.key
      ) (builtins.genericClosure {
        startSet = map (dependency: {key = dependency;}) (deps.${dependent} or []);
        operator = entry:
          if isBuiltJob entry.key
          then []
          else map (dependency: {key = dependency;}) (deps.${entry.key} or []);
      })));
  buildableDependencies = deps:
    builtins.listToAttrs (
      builtins.filter (entry: entry.value != []) (map (dependent: {
          name = dependent;
          value = builtDependencyFrontierFor deps dependent;
        })
        (builtins.filter isBuiltJob (builtins.attrNames deps)))
    );
  projectedDependencies = buildableDependencies validatedDependencies;
  projectedDependencyClosureFor = dependent:
    map (entry: entry.key) (builtins.genericClosure {
      startSet = map (dependency: {key = dependency;}) (projectedDependencies.${dependent} or []);
      operator = entry:
        map (dependency: {key = dependency;}) (projectedDependencies.${entry.key} or []);
    });
  projectedDependencyCycles = builtins.filter (
    attr: builtins.elem attr (projectedDependencyClosureFor attr)
  ) (builtins.attrNames projectedDependencies);
  validatedProjectedDependencies =
    assert projectedDependencyCycles == []
    || builtins.throw "projected NixCI dependency graph contains cycles through: ${builtins.concatStringsSep ", " projectedDependencyCycles}";
      projectedDependencies;
in {
  systems = [system];
  inherit doNotBuild;
  fail-fast = false;
  fail-on-dangling-dependencies = true;
  # Keep dependency discovery manual. With generated Rust outputs,
  # automatic discovery asks NixCI to compute derivation paths for many
  # package/check attrs during `show`, including attrs listed in doNotBuild.
  # The manual graph below uses the Hakari workspace-hack cache artifact as the
  # root for Symbolica-containing cache jobs and orders nextest archive jobs
  # after the package-local test-binary artifacts that the archives reuse. The
  # exported artifact attrs are revision-scoped symlink barriers around stable
  # Cargo artifacts, so a memoized top-level result cannot release consumers
  # without one worker realizing and publishing the underlying closure. The
  # graph is constructed over the full crate/artifact DAG so the drift and
  # cycle asserts stay meaningful, then hidden paths are contracted to their
  # nearest built producer because NixCI rejects edges to jobs it does not
  # build. Ordinary crate package attrs are not CI roots, so test-binary
  # generation can start before unrelated final package outputs.
  # See https://nix-ci.com/documentation/automatic-dependency-discovery
  # and https://nix-ci.com/documentation/manually-specified-dependencies
  dependency-discovery.enable = false;
  dependencies = validatedProjectedDependencies;
  test = {
    gammaloop-doctest = {
      package = "packages.${system}.nix-ci-check-gammaloop-doctest";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest-core = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest-core";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest-clinnet = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest-clinnet";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest-feynkit = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest-feynkit";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest-integration = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest-integration";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest-python-api = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest-python-api";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest-linnet = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest-linnet";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest-spenso = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest-spenso";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest-vakint = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest-vakint";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };
  };
  deploy = {
    ci-passed = {
      package = "packages.${system}.nix-ci-passed";
      system = system;
      branches = "all";
    };
  };
}
