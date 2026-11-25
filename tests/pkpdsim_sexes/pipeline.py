"""The pipeline of the toy model per sex, `pipeline` of its model object.

Every step writes below `{results}`, so `run_pipeline -r DIR` moves the whole
run; the processes of a step are `{cores}`.
"""

from pkpdsim.pipeline import Pipeline, Step

PIPELINE = Pipeline(
    description="The experiments and the local sensitivity of the toy model per sex.",
    steps=(
        Step(
            "experiments",
            "run_simulation",
            ("--action", "list_experiments"),
            description="list the experiments",
        ),
        Step(
            "sensitivity",
            "run_sensitivity",
            ("local", "--workers", "{cores}", "-o", "{results}/sensitivity"),
            description="the local sensitivity per sex",
        ),
    ),
)
