# Licenses

A model is published under two licenses: the model and its documentation under CC BY-SA 4.0, the source code which creates and simulates it under MIT. The license texts, the copyright argument, the functions `write_license` and `license_text` and the files they write are described in [Licenses](https://matthiaskoenig.github.io/sbmlsim/licenses/) of sbmlsim ([`sbmlsim.licenses`](https://matthiaskoenig.github.io/sbmlsim/api/licenses/)), and the licenses travel with the model in its [COMBINE archive](https://combinearchive.org). This page keeps the archive of the models of the PBPK library.

## The COMBINE archive of a model

A model of the PBPK library writes the texts and the archive with one call, `pkpdsim.pbpk.create_archive(definition, model_output_dir, model=..., copyright="2024-2026 Matthias König", readme=...)` (the licenses are written next to the models, the archive `<model>_model.omex` has the licenses, the README and the models as entries, the flat model as master). The factory of another model writes the texts next to its models and adds them to the entries of the archive, here with [pymetadata](https://github.com/matthiaskoenig/pymetadata):

```python
# not executed: a fragment of the factory of a model
from pymetadata.omex import EntryFormat, ManifestEntry
from sbmlsim.licenses import License, write_license

entries = {
    license.value: {
        "path": write_license(results_dir, license, copyright=copyright),
        "entry": ManifestEntry(
            location=f"./{license.filename}", format=EntryFormat.TXT, master=False
        ),
    }
    for license, copyright in [
        (License.CC_BY_SA_4_0, None),
        (License.MIT, "2024-2025 Matthias König"),
    ]
}
```
